#!/usr/bin/env python3
"""Compare ProteinDF total energies with PySCF reference calculations.

This script parses ProteinDF input files (fl_Userinput), extracts molecular
geometry, orbital basis sets from data/basis2, method, and XC functional,
converts the basis sets to PySCF format, and executes PySCF SCF calculations.
The resulting SCF total energy is compared against the ProteinDF database
(pdfresults_std.db or a user-specified database).

Database reading is handled via Python's standard sqlite3 module directly,
so the script does not require external tools or older Python 2 environments.
PySCF calculations require a Python environment with pyscf installed
(specified via PYSCF_PYTHON or regress.conf).
"""

import argparse
import os
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

# Periodic table symbol to atomic number mapping
ATOMIC_NUMBERS = {
    "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "NE": 10, "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15,
    "S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21, "TI": 22,
    "V": 23, "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29,
    "ZN": 30, "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36,
}


def find_git_paths():
    """Find repository top directory and git common directory."""
    try:
        top = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    except Exception:
        top = str(Path(__file__).resolve().parent.parent)

    try:
        common = subprocess.check_output(
            ["git", "rev-parse", "--git-common-dir"], text=True
        ).strip()
    except Exception:
        common = os.path.join(top, ".git")

    return top, common


def parse_conf_file(conf_path):
    """Parse a simple KEY=VALUE bash config file."""
    config = {}
    if not os.path.isfile(conf_path):
        return config
    with open(conf_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                v = v.strip().strip("\"'")
                v = os.path.expandvars(v)
                config[k] = v
    return config


def load_basis2(basis_path, basis_name):
    """Load an orbital basis set from ProteinDF data/basis2 and convert to PySCF format.

    Basis format in basis2:
      <basis_name>
      <s_count> <p_count> <d_count> <f_count> <g_count> ...
      for each shell:
        <num_pgto>
        <exponent> <coefficient>
        ...
    """
    with open(basis_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    idx = 0
    found = False
    while idx < len(lines):
        line = lines[idx].strip()
        if line == basis_name:
            found = True
            break
        idx += 1

    if not found:
        raise ValueError(f"Basis set '{basis_name}' not found in {basis_path}")

    idx += 1
    counts_line = ""
    while idx < len(lines):
        counts_line = lines[idx].strip()
        if counts_line and not counts_line.startswith("#"):
            break
        idx += 1

    counts = [int(x) for x in counts_line.split()[:5]]
    shell_l = [0, 1, 2, 3, 4]  # s, p, d, f, g
    idx += 1

    basis_list = []
    for l_val, num_cgto in zip(shell_l, counts):
        for _ in range(num_cgto):
            while idx < len(lines):
                line = lines[idx].strip()
                if line and not line.startswith("#"):
                    break
                idx += 1
            num_pgto = int(lines[idx].split()[0])
            idx += 1

            pgtos = []
            for _ in range(num_pgto):
                while idx < len(lines):
                    line = lines[idx].strip()
                    if line and not line.startswith("#"):
                        break
                    idx += 1
                parts = lines[idx].split()
                exp_val = float(parts[0])
                c_val = float(parts[1]) if len(parts) > 1 else 1.0
                pgtos.append([exp_val, c_val])
                idx += 1
            basis_list.append([l_val] + pgtos)

    return basis_list


def parse_fl_userinput(input_path):
    """Parse ProteinDF fl_Userinput to extract calculation settings."""
    with open(input_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Geometry unit
    m_unit = re.search(r"geometry/cartesian/unit\s*=\s*(\w+)", content, re.IGNORECASE)
    unit = m_unit.group(1).lower() if m_unit else "angstrom"

    # Geometry coordinates
    m_geom = re.search(r"geometry/cartesian/input\s*=\s*\{([^}]+)\}", content, re.IGNORECASE)
    if not m_geom:
        raise ValueError(f"No geometry/cartesian/input found in {input_path}")

    atom_lines = []
    atoms = []
    for line in m_geom.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            sym = parts[0]
            x, y, z = parts[1], parts[2], parts[3]
            atom_lines.append(f"{sym} {x} {y} {z}")
            atoms.append(sym)

    # Orbital basis sets
    m_basis = re.search(r"basis-set/orbital\s*=\s*\{([^}]+)\}", content, re.IGNORECASE)
    if not m_basis:
        raise ValueError(f"No basis-set/orbital found in {input_path}")

    basis_map = {}
    for line in m_basis.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m_b = re.search(r"(\w+)\s*=\s*[\"']([^\"']+)[\"']", line)
        if m_b:
            sym, bname = m_b.group(1), m_b.group(2)
            basis_map[sym] = bname

    # Method
    m_method = re.search(r"^\s*method\s*=\s*(\w+)", content, re.MULTILINE | re.IGNORECASE)
    method_str = m_method.group(1).lower() if m_method else "rks"

    # XC potential
    m_xc = re.search(r"xc-potential\s*=\s*([^\s\n\r]+)", content, re.IGNORECASE)
    xc_str = m_xc.group(1).strip() if m_xc else ""

    # Electrons
    m_rks_e = re.search(r"method/rks/electron-number\s*=\s*(\d+)", content, re.IGNORECASE)
    m_alpha = re.search(r"method/uks/alpha_electrons\s*=\s*(\d+)", content, re.IGNORECASE)
    m_beta = re.search(r"method/uks/beta_electrons\s*=\s*(\d+)", content, re.IGNORECASE)

    # Total nuclear charge
    total_z = 0
    for sym in atoms:
        clean_sym = re.sub(r"[^A-Za-z]", "", sym).upper()
        if clean_sym in ATOMIC_NUMBERS:
            total_z += ATOMIC_NUMBERS[clean_sym]
        else:
            raise ValueError(f"Unknown element symbol '{sym}'")

    is_uks = (method_str == "uks")
    if is_uks:
        if not m_alpha or not m_beta:
            raise ValueError("UKS requires method/uks/alpha_electrons and beta_electrons")
        n_alpha = int(m_alpha.group(1))
        n_beta = int(m_beta.group(1))
        tot_elec = n_alpha + n_beta
        charge = total_z - tot_elec
        spin = n_alpha - n_beta  # 2S
    else:
        tot_elec = int(m_rks_e.group(1)) if m_rks_e else total_z
        charge = total_z - tot_elec
        spin = 0
        n_alpha = tot_elec // 2
        n_beta = tot_elec // 2

    return {
        "unit": unit,
        "atom_str": "; ".join(atom_lines),
        "atoms": atoms,
        "basis_map": basis_map,
        "method": method_str,
        "xc": xc_str,
        "charge": charge,
        "spin": spin,
        "n_alpha": n_alpha,
        "n_beta": n_beta,
        "is_uks": is_uks,
    }


def read_proteindf_db(db_path):
    """Read conditions and converged total energy from ProteinDF SQLite DB."""
    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")

    con = sqlite3.connect(db_path)
    cur = con.cursor()

    cond_row = cur.execute(
        "SELECT num_of_atoms, num_of_AOs, method, xc_functional, iterations FROM conditions"
    ).fetchone()
    if not cond_row:
        raise ValueError(f"No condition record in {db_path}")

    num_atoms, num_aos, method, xc_func, iters = cond_row

    energy_row = cur.execute(
        "SELECT energy FROM total_energies WHERE iteration = ?", (iters,)
    ).fetchone()
    if not energy_row:
        energy_row = cur.execute(
            "SELECT energy FROM total_energies ORDER BY iteration DESC LIMIT 1"
        ).fetchone()

    con.close()

    total_energy = energy_row[0] if energy_row else None
    return {
        "num_atoms": num_atoms,
        "num_aos": num_aos,
        "method": method,
        "xc_functional": xc_func,
        "iterations": iters,
        "total_energy": total_energy,
    }


def resolve_xc_functional(xc_name, vwn_variant=None, custom_xc=None):
    """Resolve the XC functional string for PySCF / libxc.

    - SVWN / SVWN5:
        Default: 'LDA_X, LDA_C_VWN' (VWN5)
        If vwn_variant is specified: 'LDA_X, LDA_C_<vwn_variant>'
    - BLYP:
        Default: 'B88, LYP'
    - B3LYP:
        Default: 'B3LYP' (libxc HYB_GGA_XC_B3LYP / VWN_RPA)
        If vwn_variant is overridden:
          '0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*LDA_C_<vwn_variant>'
    """
    if custom_xc:
        return custom_xc

    xc_upper = xc_name.upper()

    if xc_upper in ("HF", ""):
        return None

    if xc_upper in ("SVWN", "SVWN5"):
        if not vwn_variant:
            return "LDA_X, LDA_C_VWN"
        vwn_clean = vwn_variant.upper()
        if not vwn_clean.startswith("LDA_C_"):
            if vwn_clean in ("VWN", "VWN5"):
                vwn_clean = "LDA_C_VWN"
            else:
                vwn_clean = f"LDA_C_{vwn_clean}"
        return f"LDA_X, {vwn_clean}"

    if xc_upper == "BLYP":
        return "B88, LYP"

    if xc_upper == "B3LYP":
        if not vwn_variant or vwn_variant.upper() in ("DEFAULT", "VWN_RPA", "RPA"):
            return "B3LYP"
        vwn_clean = vwn_variant.upper()
        if not vwn_clean.startswith("LDA_C_"):
            if vwn_clean in ("VWN", "VWN5"):
                vwn_clean = "LDA_C_VWN"
            else:
                vwn_clean = f"LDA_C_{vwn_clean}"
        return f"0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*{vwn_clean}"

    return xc_name


def run_pyscf_calc(parsed_input, basis_file, xc_override=None, vwn_variant=None, grid_level=3):
    """Run PySCF calculation and return results."""
    from pyscf import gto, scf, dft

    pyscf_basis = {}
    for sym, bname in parsed_input["basis_map"].items():
        pyscf_basis[sym] = load_basis2(basis_file, bname)

    mol = gto.M(
        atom=parsed_input["atom_str"],
        basis=pyscf_basis,
        unit=parsed_input["unit"],
        charge=parsed_input["charge"],
        spin=parsed_input["spin"],
        cart=False,
        verbose=0,
    )

    xc_func = resolve_xc_functional(
        parsed_input["xc"], vwn_variant=vwn_variant, custom_xc=xc_override
    )
    is_hf = (xc_func is None)
    is_uks = parsed_input["is_uks"]

    if is_hf:
        calc = scf.UHF(mol) if is_uks else scf.RHF(mol)
        method_label = "UHF" if is_uks else "RHF"
    else:
        calc = dft.UKS(mol) if is_uks else dft.RKS(mol)
        calc.xc = xc_func
        calc.grids.level = grid_level
        method_label = "UKS" if is_uks else "RKS"

    calc.conv_tol = 1e-11
    calc.max_cycle = 150
    energy = calc.kernel()

    s2_val = None
    mult_val = None
    if is_uks:
        s2_val, mult_val = calc.spin_square()

    return {
        "nao": mol.nao,
        "energy": energy,
        "converged": calc.converged,
        "method_label": method_label,
        "xc_used": xc_func if not is_hf else "HF",
        "s2": s2_val,
        "mult": mult_val,
        "grid_level": grid_level if not is_hf else None,
    }


def compare_entry(entry_dir, basis_file, db_override=None, xc_override=None,
                  vwn_variant=None, grid_levels=(3,)):
    """Compare a single test entry between ProteinDF and PySCF across grid levels."""
    userinput_path = os.path.join(entry_dir, "fl_Userinput")
    db_path = db_override or os.path.join(entry_dir, "pdfresults_std.db")

    parsed_input = parse_fl_userinput(userinput_path)
    pdf_res = read_proteindf_db(db_path)

    results = []
    for g_lvl in grid_levels:
        pyscf_res = run_pyscf_calc(
            parsed_input,
            basis_file,
            xc_override=xc_override,
            vwn_variant=vwn_variant,
            grid_level=g_lvl,
        )

        diff = pyscf_res["energy"] - pdf_res["total_energy"]
        results.append({
            "entry": os.path.basename(entry_dir),
            "method": pyscf_res["method_label"],
            "xc_input": parsed_input["xc"] or "HF",
            "xc_used": pyscf_res["xc_used"],
            "grid_level": pyscf_res["grid_level"],
            "pdf_nao": pdf_res["num_aos"],
            "pyscf_nao": pyscf_res["nao"],
            "nao_match": (pdf_res["num_aos"] == pyscf_res["nao"]),
            "pdf_energy": pdf_res["total_energy"],
            "pyscf_energy": pyscf_res["energy"],
            "diff": diff,
            "converged": pyscf_res["converged"],
            "spin": parsed_input["spin"],
            "s2": pyscf_res["s2"],
            "mult": pyscf_res["mult"],
            "n_alpha": parsed_input["n_alpha"],
            "n_beta": parsed_input["n_beta"],
        })

    return results


def check_vwn_variants(suite_dir, basis_file, output_format="table"):
    """Evaluate candidate VWN variants for Ne_SVWN5 and Ne_B3LYP."""
    print("=== VWN Variant Analysis for Ne_SVWN5 and Ne_B3LYP ===")

    # 1. Ne_SVWN5
    ne_svwn5_dir = os.path.join(suite_dir, "Ne_SVWN5")
    parsed_svwn5 = parse_fl_userinput(os.path.join(ne_svwn5_dir, "fl_Userinput"))
    pdf_svwn5 = read_proteindf_db(os.path.join(ne_svwn5_dir, "pdfresults_std.db"))

    svwn_candidates = [
        ("LDA_C_VWN (VWN5)", "LDA_X, LDA_C_VWN"),
        ("LDA_C_VWN_RPA", "LDA_X, LDA_C_VWN_RPA"),
        ("LDA_C_VWN_1", "LDA_X, LDA_C_VWN_1"),
        ("LDA_C_VWN_2", "LDA_X, LDA_C_VWN_2"),
        ("LDA_C_VWN_3", "LDA_X, LDA_C_VWN_3"),
        ("LDA_C_VWN_4", "LDA_X, LDA_C_VWN_4"),
    ]

    svwn_results = []
    for label, xc_str in svwn_candidates:
        res = run_pyscf_calc(parsed_svwn5, basis_file, xc_override=xc_str, grid_level=7)
        diff = res["energy"] - pdf_svwn5["total_energy"]
        svwn_results.append({
            "candidate": label,
            "xc_str": xc_str,
            "pdf_energy": pdf_svwn5["total_energy"],
            "pyscf_energy": res["energy"],
            "diff": diff,
        })

    # 2. Ne_B3LYP
    ne_b3lyp_dir = os.path.join(suite_dir, "Ne_B3LYP")
    parsed_b3lyp = parse_fl_userinput(os.path.join(ne_b3lyp_dir, "fl_Userinput"))
    pdf_b3lyp = read_proteindf_db(os.path.join(ne_b3lyp_dir, "pdfresults_std.db"))

    b3lyp_candidates = [
        ("B3LYP (default, VWN_RPA)", "B3LYP"),
        ("B3LYP with VWN_RPA", "0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*LDA_C_VWN_RPA"),
        ("B3LYP5 (VWN5)", "B3LYP5"),
        ("B3LYP with VWN (VWN5)", "0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*LDA_C_VWN"),
        ("B3LYP with VWN_3", "0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*LDA_C_VWN_3"),
        ("B3LYP with VWN_1", "0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*LDA_C_VWN_1"),
        ("B3LYP with VWN_2", "0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*LDA_C_VWN_2"),
    ]

    b3lyp_results = []
    for label, xc_str in b3lyp_candidates:
        res = run_pyscf_calc(parsed_b3lyp, basis_file, xc_override=xc_str, grid_level=7)
        diff = res["energy"] - pdf_b3lyp["total_energy"]
        b3lyp_results.append({
            "candidate": label,
            "xc_str": xc_str,
            "pdf_energy": pdf_b3lyp["total_energy"],
            "pyscf_energy": res["energy"],
            "diff": diff,
        })

    if output_format == "markdown":
        print("#### Ne_SVWN5 (PDF: {:.10f} Eh)".format(pdf_svwn5["total_energy"]))
        print("| Candidate | PySCF XC string | PySCF Energy (Eh) | Diff (Eh) |")
        print("|:---|:---|:---|:---|")
        for r in svwn_results:
            print(f"| {r['candidate']} | `{r['xc_str']}` | {r['pyscf_energy']:.10f} | {r['diff']:+.10e} |")

        print("\n#### Ne_B3LYP (PDF: {:.10f} Eh)".format(pdf_b3lyp["total_energy"]))
        print("| Candidate | PySCF XC string | PySCF Energy (Eh) | Diff (Eh) |")
        print("|:---|:---|:---|:---|")
        for r in b3lyp_results:
            print(f"| {r['candidate']} | `{r['xc_str']}` | {r['pyscf_energy']:.10f} | {r['diff']:+.10e} |")
    else:
        print("\n--- Ne_SVWN5 ---")
        for r in svwn_results:
            print(f"{r['candidate']:24s}: PySCF={r['pyscf_energy']:.10f}  Diff={r['diff']:+.10e}")
        print("\n--- Ne_B3LYP ---")
        for r in b3lyp_results:
            print(f"{r['candidate']:30s}: PySCF={r['pyscf_energy']:.10f}  Diff={r['diff']:+.10e}")


def print_results_table(results_list, output_format="table", show_spin=False):
    """Print results in plain table or markdown format."""
    if output_format == "markdown":
        if show_spin:
            print("| Entry | Method | 2S (α, β) | <S^2> | Mult (2S+1) | PySCF Energy (Eh) | PDF Energy (Eh) | Diff (Eh) |")
            print("|:---|:---:|:---:|:---:|:---:|:---|:---|:---|")
            for r in results_list:
                spin_str = f"{r['spin']} ({r['n_alpha']}, {r['n_beta']})"
                s2_str = f"{r['s2']:.4f}" if r["s2"] is not None else "-"
                mult_str = f"{r['mult']:.4f}" if r["mult"] is not None else "-"
                print(f"| {r['entry']} | {r['method']} | {spin_str} | {s2_str} | {mult_str} | {r['pyscf_energy']:.10f} | {r['pdf_energy']:.10f} | {r['diff']:+.10e} |")
        else:
            print("| Entry | Method | XC (PySCF) | Grid | AO (PDF/PySCF) | PDF Energy (Eh) | PySCF Energy (Eh) | Diff (Eh) | Conv |")
            print("|:---|:---|:---|:---:|:---:|:---|:---|:---|:---:|")
            for r in results_list:
                grid_str = str(r["grid_level"]) if r["grid_level"] is not None else "-"
                ao_str = f"{r['pdf_nao']}/{r['pyscf_nao']}"
                conv_str = "YES" if r["converged"] else "NO"
                print(f"| {r['entry']} | {r['method']} | {r['xc_used']} | {grid_str} | {ao_str} | {r['pdf_energy']:.10f} | {r['pyscf_energy']:.10f} | {r['diff']:+.10e} | {conv_str} |")
    else:
        header = f"{'Entry':16s} {'Method':6s} {'XC (PySCF)':22s} {'Grid':4s} {'AO(PDF/Py)':10s} {'PDF Energy (Eh)':18s} {'PySCF Energy (Eh)':18s} {'Diff (Eh)':16s} {'Conv':4s}"
        print(header)
        print("-" * len(header))
        for r in results_list:
            grid_str = str(r["grid_level"]) if r["grid_level"] is not None else "-"
            ao_str = f"{r['pdf_nao']}/{r['pyscf_nao']}"
            conv_str = "YES" if r["converged"] else "NO"
            print(f"{r['entry']:16s} {r['method']:6s} {r['xc_used'][:22]:22s} {grid_str:4s} {ao_str:10s} {r['pdf_energy']:18.10f} {r['pyscf_energy']:18.10f} {r['diff']:+16.10e} {conv_str:4s}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare ProteinDF SCF total energy against PySCF reference calculations."
    )
    parser.add_argument("--suite", default="serial_dev", help="Test suite name (default: serial_dev)")
    parser.add_argument("--entry", help="Single test entry name (e.g. Ne_B3LYP)")
    parser.add_argument("--entries", help="Comma-separated test entry names")
    parser.add_argument("--all", action="store_true", help="Run all entries in suite")
    parser.add_argument("--test-dir", help="Path to ProteinDF_test directory")
    parser.add_argument("--basis-file", help="Path to data/basis2 file")
    parser.add_argument("--db", help="Path to specific calculation DB (default: pdfresults_std.db in entry)")
    parser.add_argument("--xc", help="Override XC functional string for PySCF")
    parser.add_argument("--vwn", default=None, help="VWN variant for SVWN/B3LYP (options: VWN, VWN_RPA, VWN_3, etc.)")
    parser.add_argument("--grid-level", default="3", help="Comma-separated PySCF grid levels (default: 3)")
    parser.add_argument("--check-vwn", action="store_true", help="Run analysis of candidate VWN variants on Ne_SVWN5 and Ne_B3LYP")
    parser.add_argument("--show-spin", action="store_true", help="Display spin multiplet and <S^2> details")
    parser.add_argument("--format", choices=["table", "markdown"], default="table", help="Output format")
    args = parser.parse_args()

    top_dir, common_dir = find_git_paths()
    conf_path = os.path.join(common_dir, "regress.conf")
    conf = parse_conf_file(conf_path)

    test_dir = args.test_dir or os.environ.get("PROTEINDF_TEST_DIR") or conf.get("PROTEINDF_TEST_DIR")
    if not test_dir or not os.path.isdir(test_dir):
        sys.exit(f"Error: ProteinDF_test directory not found ({test_dir}). Please set PROTEINDF_TEST_DIR.")

    basis_file = args.basis_file or os.path.join(top_dir, "data", "basis2")
    if not os.path.isfile(basis_file):
        sys.exit(f"Error: Basis set file not found: {basis_file}")

    suite_dir = os.path.join(test_dir, args.suite)
    if not os.path.isdir(suite_dir):
        sys.exit(f"Error: Suite directory not found: {suite_dir}")

    if args.check_vwn:
        check_vwn_variants(suite_dir, basis_file, output_format=args.format)
        return

    # Determine entries to run
    if args.entry:
        entries = [args.entry]
    elif args.entries:
        entries = [e.strip() for e in args.entries.split(",") if e.strip()]
    elif args.all:
        entries = sorted([
            d for d in os.listdir(suite_dir)
            if os.path.isdir(os.path.join(suite_dir, d))
            and os.path.isfile(os.path.join(suite_dir, d, "fl_Userinput"))
        ])
    else:
        entries = ["Ne_HF"]

    grid_levels = [int(x.strip()) for x in args.grid_level.split(",") if x.strip()]

    all_results = []
    for entry in entries:
        entry_dir = os.path.join(suite_dir, entry)
        if not os.path.isdir(entry_dir):
            print(f"Warning: Entry directory not found: {entry_dir}, skipping.", file=sys.stderr)
            continue
        try:
            res = compare_entry(
                entry_dir,
                basis_file,
                db_override=args.db,
                xc_override=args.xc,
                vwn_variant=args.vwn,
                grid_levels=grid_levels,
            )
            all_results.extend(res)
        except Exception as e:
            print(f"Error processing entry {entry}: {e}", file=sys.stderr)

    if all_results:
        print_results_table(all_results, output_format=args.format, show_spin=args.show_spin)


if __name__ == "__main__":
    main()
