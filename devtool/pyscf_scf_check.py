#!/usr/bin/env python3
"""Inspect ProteinDF SCF convergence and orbital energies against PySCF references.

This tool compares converged ProteinDF density matrices, orbital energies,
and SCF stationarity conditions against PySCF. It maps atomic orbital (AO)
orderings between ProteinDF and PySCF, evaluates the stationary condition
commutator [F, P*S], and compares variational energy and orbital energies.
"""

import argparse
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# Import helper functions from pyscf_compare
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pyscf_compare


def load_symmetric_matrix(path):
    """Load a packed symmetric matrix from ProteinDF binary file (.mat)."""
    with open(path, "rb") as f:
        data = f.read()

    filesize = len(data)
    endian_list = ["=", "<", ">"]
    header_struct_list = ["bii", "iii", "bll", "ill"]

    found = False
    for endian in endian_list:
        for hs in header_struct_list:
            fmt = endian + hs
            sz = struct.calcsize(fmt)
            if sz > filesize:
                continue
            h = struct.unpack(fmt, data[:sz])
            mtype, row, col = h[0], h[1], h[2]
            expected = sz + 8 * row * (row + 1) // 2
            if expected == filesize:
                found = True
                header_fmt = fmt
                header_sz = sz
                break
        if found:
            break

    if not found:
        raise ValueError(f"Cannot parse symmetric matrix file: {path} (size: {filesize})")

    vals = struct.unpack(f"{endian}{row*(row+1)//2}d", data[header_sz:])
    mat = np.zeros((row, col), dtype=np.float64)
    idx = 0
    for r in range(row):
        for c in range(r + 1):
            val = vals[idx]
            mat[r, c] = val
            mat[c, r] = val
            idx += 1

    return mat


def load_vector(path):
    """Load a vector from ProteinDF binary file (.vtr)."""
    with open(path, "rb") as f:
        data = f.read()

    filesize = len(data)
    endian_list = ["=", "<", ">"]
    header_struct_list = ["i", "l", "q"]

    found = False
    for endian in endian_list:
        for hs in header_struct_list:
            fmt = endian + hs
            sz = struct.calcsize(fmt)
            if sz > filesize:
                continue
            h = struct.unpack(fmt, data[:sz])
            n = h[0]
            expected = sz + 8 * n
            if expected == filesize:
                found = True
                header_fmt = fmt
                header_sz = sz
                break
        if found:
            break

    if not found:
        raise ValueError(f"Cannot parse vector file: {path} (size: {filesize})")

    vals = struct.unpack(f"{endian}{n}d", data[header_sz:])
    return np.array(vals, dtype=np.float64)


def build_ao_mapping(parsed_input, basis_file):
    """Build index permutation mapping from ProteinDF AO ordering to PySCF AO ordering."""
    shell_map = {
        0: [0],                      # s: s
        1: [0, 1, 2],                # p: px, py, pz
        2: [0, 3, 1, 4, 2],          # d: dxy, dxz, dyz, dx2-y2, dz2 -> dxy, dyz, dz2, dxz, dx2-y2
    }

    p2p = []
    ao_offset = 0
    for sym in parsed_input["atoms"]:
        bname = parsed_input["basis_map"][sym]
        basis_list = pyscf_compare.load_basis2(basis_file, bname)
        for shell in basis_list:
            l_val = shell[0]
            if l_val not in shell_map:
                raise NotImplementedError(f"Angular momentum l={l_val} not supported in AO mapping")
            mapping = shell_map[l_val]
            for local_idx in mapping:
                p2p.append(ao_offset + local_idx)
            ao_offset += len(mapping)
    return p2p


def permute_matrix(mat, p2p):
    """Permute rows and columns of a matrix according to mapping p2p (PDF to PySCF)."""
    nao = len(p2p)
    perm = np.zeros((nao, nao), dtype=np.float64)
    for i, j in enumerate(p2p):
        perm[j, i] = 1.0
    return perm @ mat @ perm.T


def find_pdf_home(top_dir):
    """Locate ProteinDF installation directory."""
    env_pdf = os.environ.get("PDF_HOME", "")
    if env_pdf and os.path.isfile(os.path.join(env_pdf, "bin", "PDF.x")):
        return env_pdf

    candidates = [
        os.path.join(top_dir, "build-regress", "install"),
        os.path.join(top_dir, "build-check", "install"),
        os.path.join(top_dir, "install"),
    ]
    for c in candidates:
        if os.path.isfile(os.path.join(c, "bin", "PDF.x")):
            return c

    # Search in PATH
    pdf_x = shutil.which("PDF.x")
    if pdf_x:
        return str(Path(pdf_x).resolve().parent.parent)

    return None


def modify_userinput_for_tight_conv(content, threshold=1e-8, threshold_energy=1e-10, cut_value=1e-16):
    """Update convergence criteria in fl_Userinput content."""
    updates = {
        r"convergence/threshold\s*=.*": f"convergence/threshold\t= {threshold}",
        r"convergence/threshold-energy\s*=.*": f"convergence/threshold-energy\t= {threshold_energy}",
        r"cut-value\s*=.*": f"cut-value\t= {cut_value}",
        r"CDAM_tau\s*=.*": "CDAM_tau = 1.0E-16",
        r"CD_epsilon\s*=.*": "CD_epsilon = 1.0E-10",
        r"max-iteration\s*=.*": "max-iteration\t= 200",
    }
    new_content = content
    for pat, repl in updates.items():
        if re.search(pat, new_content):
            new_content = re.sub(pat, repl, new_content)
        else:
            # If not found, insert under SCF section if possible
            if "convergence/threshold" in pat:
                new_content = re.sub(r"(>>>>SCF\n)", rf"\1\t{repl}\n", new_content)

    return new_content


def run_proteindf(entry_src_dir, work_dir, pdf_home, tight=False, verbose=False):
    """Run ProteinDF in work_dir and return calculation summary."""
    userinput_src = os.path.join(entry_src_dir, "fl_Userinput")
    with open(userinput_src, "r", encoding="utf-8") as f:
        content = f.read()

    if tight:
        content = modify_userinput_for_tight_conv(content)

    userinput_dst = os.path.join(work_dir, "fl_Userinput")
    with open(userinput_dst, "w", encoding="utf-8") as f:
        f.write(content)

    fl_work = os.path.join(work_dir, "fl_Work")
    os.makedirs(fl_work, exist_ok=True)

    env = os.environ.copy()
    env["PDF_HOME"] = pdf_home
    env["PATH"] = f"{os.path.join(pdf_home, 'bin')}:{env.get('PATH', '')}"
    if "OMP_NUM_THREADS" not in env:
        env["OMP_NUM_THREADS"] = str(os.cpu_count() or 4)

    pdf_bin = os.path.join(pdf_home, "bin", "PDF.x")
    log_file = os.path.join(work_dir, "pdf.log")

    with open(log_file, "w", encoding="utf-8") as out:
        proc = subprocess.run([pdf_bin], cwd=work_dir, env=env, stdout=out, stderr=subprocess.STDOUT)

    if proc.returncode != 0:
        with open(log_file, "r", encoding="utf-8") as f:
            tail = "".join(f.readlines()[-30:])
        raise RuntimeError(f"ProteinDF failed with return code {proc.returncode}:\n{tail}")

    # Parse log for iteration and total energy
    converged_iter = None
    total_energy = None
    with open(log_file, "r", encoding="utf-8") as f:
        for line in f:
            m_te = re.search(r"^\s*(\d+)\s*th TE\s*=\s*([-\d.]+)", line)
            if m_te:
                converged_iter = int(m_te.group(1))
                total_energy = float(m_te.group(2))

    if converged_iter is None:
        raise RuntimeError(f"Could not determine converged iteration from {log_file}")

    return {
        "iteration": converged_iter,
        "total_energy": total_energy,
        "log_file": log_file,
    }


def extract_proteindf_data(work_dir, is_uks, iteration):
    """Read converged matrices and vectors from fl_Work directory."""
    fl_work = os.path.join(work_dir, "fl_Work")
    spq_path = os.path.join(fl_work, "Spq.mat")
    if not os.path.isfile(spq_path):
        raise FileNotFoundError(f"Overlap matrix not found: {spq_path}")
    pdf_S = load_symmetric_matrix(spq_path)

    if is_uks:
        pa_path = os.path.join(fl_work, f"Ppq.uks_alpha{iteration}.mat")
        pb_path = os.path.join(fl_work, f"Ppq.uks_beta{iteration}.mat")
        ea_path = os.path.join(fl_work, f"eigenvalues.uks_alpha{iteration}.vtr")
        eb_path = os.path.join(fl_work, f"eigenvalues.uks_beta{iteration}.vtr")

        p_alpha = load_symmetric_matrix(pa_path)
        p_beta = load_symmetric_matrix(pb_path)
        mo_energy_a = load_vector(ea_path)
        mo_energy_b = load_vector(eb_path)

        return {
            "S": pdf_S,
            "P_alpha": p_alpha,
            "P_beta": p_beta,
            "mo_energy_alpha": mo_energy_a,
            "mo_energy_beta": mo_energy_b,
        }
    else:
        p_path = os.path.join(fl_work, f"Ppq.rks{iteration}.mat")
        e_path = os.path.join(fl_work, f"eigenvalues.rks{iteration}.vtr")

        p_tot = load_symmetric_matrix(p_path)
        mo_energy = load_vector(e_path)

        return {
            "S": pdf_S,
            "P_tot": p_tot,
            "mo_energy": mo_energy,
        }


def check_scf_entry(entry_dir, basis_file, pdf_home, tight=False,
                    grid_levels=(5, 7), keep_work=False, external_work_dir=None):
    """Run SCF check on a single entry and return structured comparison results."""
    from pyscf import gto, scf, dft

    userinput_path = os.path.join(entry_dir, "fl_Userinput")
    parsed_input = pyscf_compare.parse_fl_userinput(userinput_path)
    is_uks = parsed_input["is_uks"]

    # 1. Run ProteinDF or use existing work directory
    work_dir_obj = None
    if external_work_dir:
        work_dir = external_work_dir
        # Determine iteration from fl_Work
        fl_work = os.path.join(work_dir, "fl_Work")
        itr_candidates = []
        pattern = re.compile(r"Ppq\.(rks|uks_alpha)(\d+)\.mat")
        for f in os.listdir(fl_work):
            m = pattern.match(f)
            if m:
                itr_candidates.append(int(m.group(2)))
        iteration = max(itr_candidates) if itr_candidates else 1
        pdf_calc = {"iteration": iteration, "total_energy": None}
    else:
        if keep_work:
            work_dir = tempfile.mkdtemp(prefix=f"pdf-check-{os.path.basename(entry_dir)}-")
        else:
            work_dir_obj = tempfile.TemporaryDirectory(prefix=f"pdf-check-{os.path.basename(entry_dir)}-")
            work_dir = work_dir_obj.name

        pdf_calc = run_proteindf(entry_dir, work_dir, pdf_home, tight=tight)
        iteration = pdf_calc["iteration"]

    try:
        # 2. Extract ProteinDF results
        pdf_data = extract_proteindf_data(work_dir, is_uks, iteration)

        # 3. Build AO mapping and map density matrices
        p2p = build_ao_mapping(parsed_input, basis_file)
        if is_uks:
            dm_pdf_a = permute_matrix(pdf_data["P_alpha"], p2p)
            dm_pdf_b = permute_matrix(pdf_data["P_beta"], p2p)
            dm_pdf = np.array([dm_pdf_a, dm_pdf_b])
        else:
            dm_pdf_tot = permute_matrix(pdf_data["P_tot"], p2p)
            dm_pdf = dm_pdf_tot

        # 4. Set up PySCF Molecule
        pyscf_basis = {}
        for sym, bname in parsed_input["basis_map"].items():
            pyscf_basis[sym] = pyscf_compare.load_basis2(basis_file, bname)

        mol = gto.M(
            atom=parsed_input["atom_str"],
            basis=pyscf_basis,
            unit=parsed_input["unit"],
            charge=parsed_input["charge"],
            spin=parsed_input["spin"],
            cart=False,
            verbose=0,
        )

        xc_func = pyscf_compare.resolve_xc_functional(parsed_input["xc"])
        is_hf = (xc_func is None)

        results = {
            "entry": os.path.basename(entry_dir),
            "method": "UKS" if (is_uks and not is_hf) else ("UHF" if is_uks else ("RKS" if not is_hf else "RHF")),
            "xc": parsed_input["xc"] or "HF",
            "is_uks": is_uks,
            "is_hf": is_hf,
            "tight": tight,
            "pdf_iter": iteration,
            "pdf_energy": pdf_calc["total_energy"],
            "grids": [],
        }

        # Select grid levels
        lvls = [None] if is_hf else grid_levels

        for g_lvl in lvls:
            if is_hf:
                mf = scf.UHF(mol) if is_uks else scf.RHF(mol)
            else:
                mf = dft.UKS(mol) if is_uks else dft.RKS(mol)
                mf.xc = xc_func
                mf.grids.level = g_lvl

            mf.conv_tol = 1e-12
            mf.max_cycle = 150
            e_pyscf = mf.kernel()
            dm_pyscf = mf.make_rdm1()
            pyscf_S = mol.intor("int1e_ovlp")

            # Evaluate Fock and energy for ProteinDF density
            h1e = mf.get_hcore()
            vhf_pdf = mf.get_veff(mol, dm_pdf)
            fock_pdf = mf.get_fock(h1e, vhf=vhf_pdf, dm=dm_pdf)
            e_pdf_in_pyscf = mf.energy_tot(dm=dm_pdf, vhf=vhf_pdf)
            fock_pyscf = mf.get_fock(dm=dm_pyscf)

            delta_e = e_pdf_in_pyscf - e_pyscf

            # Stationarity condition: [F_sigma, P_sigma * S]
            residuals = {}
            if is_uks:
                spins = [("alpha", 0), ("beta", 1)]
                for sname, sidx in spins:
                    F_p = fock_pdf[sidx]
                    P_p = dm_pdf[sidx]
                    res_p = F_p @ P_p @ pyscf_S - pyscf_S @ P_p @ F_p

                    F_y = fock_pyscf[sidx]
                    P_y = dm_pyscf[sidx]
                    res_y = F_y @ P_y @ pyscf_S - pyscf_S @ P_y @ F_y

                    residuals[sname] = {
                        "pdf_max_abs": float(np.max(np.abs(res_p))),
                        "pdf_frob": float(np.linalg.norm(res_p)),
                        "pyscf_max_abs": float(np.max(np.abs(res_y))),
                        "pyscf_frob": float(np.linalg.norm(res_y)),
                    }
            else:
                # For RKS/RHF, evaluate with P_sigma = P / 2
                P_p_sigma = dm_pdf / 2.0
                F_p = fock_pdf
                res_p = F_p @ P_p_sigma @ pyscf_S - pyscf_S @ P_p_sigma @ F_p

                P_y_sigma = dm_pyscf / 2.0
                F_y = fock_pyscf
                res_y = F_y @ P_y_sigma @ pyscf_S - pyscf_S @ P_y_sigma @ F_y

                residuals["tot"] = {
                    "pdf_max_abs": float(np.max(np.abs(res_p))),
                    "pdf_frob": float(np.linalg.norm(res_p)),
                    "pyscf_max_abs": float(np.max(np.abs(res_y))),
                    "pyscf_frob": float(np.linalg.norm(res_y)),
                }

            # Orbital energies around HOMO
            mo_data = {}
            if is_uks:
                n_a = parsed_input["n_alpha"]
                n_b = parsed_input["n_beta"]
                mo_data["alpha"] = {
                    "pdf_occ": pdf_data["mo_energy_alpha"][max(0, n_a - 3):n_a].tolist(),
                    "pyscf_occ": mf.mo_energy[0][max(0, n_a - 3):n_a].tolist(),
                    "pdf_vir": pdf_data["mo_energy_alpha"][n_a:n_a + 2].tolist(),
                    "pyscf_vir": mf.mo_energy[0][n_a:n_a + 2].tolist(),
                }
                mo_data["beta"] = {
                    "pdf_occ": pdf_data["mo_energy_beta"][max(0, n_b - 3):n_b].tolist(),
                    "pyscf_occ": mf.mo_energy[1][max(0, n_b - 3):n_b].tolist(),
                    "pdf_vir": pdf_data["mo_energy_beta"][n_b:n_b + 2].tolist(),
                    "pyscf_vir": mf.mo_energy[1][n_b:n_b + 2].tolist(),
                }
            else:
                n_occ = parsed_input["n_alpha"]
                mo_data["tot"] = {
                    "pdf_occ": pdf_data["mo_energy"][max(0, n_occ - 3):n_occ].tolist(),
                    "pyscf_occ": mf.mo_energy[max(0, n_occ - 3):n_occ].tolist(),
                    "pdf_vir": pdf_data["mo_energy"][n_occ:n_occ + 2].tolist(),
                    "pyscf_vir": mf.mo_energy[n_occ:n_occ + 2].tolist(),
                }

            results["grids"].append({
                "grid_level": g_lvl,
                "e_pyscf": float(e_pyscf),
                "e_pdf_in_pyscf": float(e_pdf_in_pyscf),
                "delta_e": float(delta_e),
                "residuals": residuals,
                "mo_energies": mo_data,
            })

        return results

    finally:
        if work_dir_obj:
            work_dir_obj.cleanup()


def format_report(results_list, mode="full"):
    """Format structured results into report tables."""
    lines = []
    lines.append("### (b) SCF Stationarity Condition & (c) Variational Energy Difference")
    lines.append("")
    lines.append("| Entry | Conv | Grid | E[P_PDF] (PySCF) | E[P_PySCF] | Delta E (a.u.) | Spin | PDF res max | PDF res frob | PySCF res max | PySCF res frob |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")

    for r in results_list:
        entry = r["entry"]
        conv_label = "tight" if r["tight"] else "default"
        for g in r["grids"]:
            glvl = str(g["grid_level"]) if g["grid_level"] is not None else "-"
            e_pdf = f"{g['e_pdf_in_pyscf']:.10f}"
            e_pyscf = f"{g['e_pyscf']:.10f}"
            de = f"{g['delta_e']:+.2e}"

            first = True
            for spin, res in g["residuals"].items():
                p_max = f"{res['pdf_max_abs']:.2e}"
                p_fr = f"{res['pdf_frob']:.2e}"
                y_max = f"{res['pyscf_max_abs']:.2e}"
                y_fr = f"{res['pyscf_frob']:.2e}"
                if first:
                    lines.append(f"| {entry} | {conv_label} | {glvl} | {e_pdf} | {e_pyscf} | {de} | {spin} | {p_max} | {p_fr} | {y_max} | {y_fr} |")
                    first = False
                else:
                    lines.append(f"| | | | | | | {spin} | {p_max} | {p_fr} | {y_max} | {y_fr} |")

    lines.append("")
    lines.append("### (a) Orbital Energies around HOMO (Occ top 3, Vir low 2)")
    lines.append("")
    for r in results_list:
        entry = r["entry"]
        conv_label = "tight" if r["tight"] else "default"
        lines.append(f"#### {entry} ({conv_label})")
        # Use first grid level
        g = r["grids"][0]
        for spin, mo in g["mo_energies"].items():
            lines.append(f"- **Spin: {spin}**")
            p_occ = [f"{x:.6f}" for x in mo["pdf_occ"]]
            y_occ = [f"{x:.6f}" for x in mo["pyscf_occ"]]
            p_vir = [f"{x:.6f}" for x in mo["pdf_vir"]]
            y_vir = [f"{x:.6f}" for x in mo["pyscf_vir"]]
            lines.append(f"  - PDF   occ: `[{', '.join(p_occ)}]`, vir: `[{', '.join(p_vir)}]`")
            lines.append(f"  - PySCF occ: `[{', '.join(y_occ)}]`, vir: `[{', '.join(y_vir)}]`")
            diff_occ = [f"{(po - yo):+.2e}" for po, yo in zip(mo["pdf_occ"], mo["pyscf_occ"])]
            diff_vir = [f"{(pv - yv):+.2e}" for pv, yv in zip(mo["pdf_vir"], mo["pyscf_vir"])]
            lines.append(f"  - Diff  occ: `[{', '.join(diff_occ)}]`, vir: `[{', '.join(diff_vir)}]`")
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect ProteinDF SCF convergence and orbital energies against PySCF references."
    )
    parser.add_argument("--entry", help="Single test entry name (e.g. O2_UBLYP)")
    parser.add_argument("--entries", help="Comma-separated test entry names")
    parser.add_argument("--suite", default="serial_dev", help="Test suite name (default: serial_dev)")
    parser.add_argument("--test-dir", help="Path to ProteinDF_test directory")
    parser.add_argument("--basis-file", help="Path to data/basis2 file")
    parser.add_argument("--pdf-home", help="Path to ProteinDF install directory")
    parser.add_argument("--tight", action="store_true", help="Run with tight convergence settings")
    parser.add_argument("--both-conv", action="store_true", help="Run with both default and tight convergence")
    parser.add_argument("--grid-levels", default="5,7", help="Comma-separated PySCF grid levels (default: 5,7)")
    parser.add_argument("--keep-work", action="store_true", help="Keep temporary work directories")
    parser.add_argument("--work-dir", help="Use existing work directory (skip ProteinDF calculation)")

    args = parser.parse_args()

    top_dir, common_dir = pyscf_compare.find_git_paths()
    conf = pyscf_compare.parse_conf_file(os.path.join(common_dir, "regress.conf"))

    test_dir = args.test_dir or os.environ.get("PROTEINDF_TEST_DIR") or conf.get("PROTEINDF_TEST_DIR")
    if not test_dir or not os.path.isdir(test_dir):
        sys.exit(f"ERROR: ProteinDF_test directory not found: {test_dir}")

    basis_file = args.basis_file or os.path.join(top_dir, "data", "basis2")
    if not os.path.isfile(basis_file):
        sys.exit(f"ERROR: Basis file not found: {basis_file}")

    pdf_home = args.pdf_home or find_pdf_home(top_dir)
    if not args.work_dir and (not pdf_home or not os.path.isfile(os.path.join(pdf_home, "bin", "PDF.x"))):
        sys.exit(f"ERROR: PDF.x not found in: {pdf_home}")

    grid_levels = [int(x.strip()) for x in args.grid_levels.split(",") if x.strip()]

    target_entries = []
    if args.entry:
        target_entries.append(args.entry)
    elif args.entries:
        target_entries.extend([x.strip() for x in args.entries.split(",") if x.strip()])
    else:
        target_entries = ["O2_UBLYP", "O2_USVWN5", "O2_UHF", "N2_RBLYP", "N2_UBLYP", "O2_UB3LYP"]

    suite_dir = os.path.join(test_dir, args.suite)
    all_results = []

    conv_modes = [False, True] if args.both_conv else [args.tight]

    for entry_name in target_entries:
        entry_dir = os.path.join(suite_dir, entry_name)
        if not os.path.isdir(entry_dir):
            print(f"WARNING: Entry not found: {entry_dir}", file=sys.stderr)
            continue

        for tight in conv_modes:
            mode_str = "tight" if tight else "default"
            print(f"--- Running {entry_name} ({mode_str} conv) ---", file=sys.stderr)
            res = check_scf_entry(
                entry_dir,
                basis_file,
                pdf_home,
                tight=tight,
                grid_levels=grid_levels,
                keep_work=args.keep_work,
                external_work_dir=args.work_dir,
            )
            all_results.append(res)

    print(format_report(all_results))


if __name__ == "__main__":
    main()
