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
import math
import os
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

try:
    import h5py

    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

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
    mm_charges = []
    for line in m_geom.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            sym = parts[0]
            if sym.upper() == "X":
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                chg = float(parts[4]) if len(parts) >= 5 else 0.0
                mm_charges.append((x, y, z, chg))
            else:
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
        "mm_charges": mm_charges,
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


def read_proteindf_h5(h5_path):
    """Read conditions and converged total energy from ProteinDF reference.h5."""
    if not HAS_H5PY:
        raise RuntimeError("h5py is not installed")
    if not os.path.isfile(h5_path):
        raise FileNotFoundError(f"HDF5 file not found: {h5_path}")

    with h5py.File(h5_path, "r") as f:
        # Total energy: last value of TEs
        if "TEs" in f and len(f["TEs"]) > 0:
            total_energy = float(f["TEs"][-1])
        else:
            total_energy = None

        def _get_attr(attr_name):
            val = f.attrs.get(attr_name)
            if val is None and "control" in f:
                val = f["control"].attrs.get(attr_name)
            if val is not None:
                if isinstance(val, (bytes, bytearray)):
                    val = val.decode("utf-8")
                return val
            return None

        num_atoms = _get_attr("num_of_atoms")
        if num_atoms is None and "molecule/xyzc" in f:
            num_atoms = len(f["molecule/xyzc"])
        if num_atoms is not None:
            num_atoms = int(num_atoms)

        num_aos = _get_attr("num_of_AOs")
        if num_aos is not None:
            num_aos = int(num_aos)

        method = _get_attr("method")
        if method is not None:
            method = str(method)

        xc_func = _get_attr("xc_functional")
        if xc_func is not None:
            xc_func = str(xc_func)

        iters = _get_attr("iterations")
        if iters is None and "TEs" in f:
            iters = len(f["TEs"])
        if iters is not None:
            iters = int(iters)

    return {
        "num_atoms": num_atoms,
        "num_aos": num_aos,
        "method": method,
        "xc_functional": xc_func,
        "iterations": iters,
        "total_energy": total_energy,
    }


def read_proteindf_reference(entry_dir, ref_mode="auto", db_override=None):
    """Read ProteinDF reference data from HDF5 or SQLite DB.

    ref_mode can be 'auto', 'h5', or 'db'.
    If db_override is given, SQLite DB is used directly.
    Returns a tuple of (data_dict, ref_format_used).
    """
    if db_override:
        return read_proteindf_db(db_override), "db"

    h5_path = os.path.join(entry_dir, "reference.h5")
    db_path = os.path.join(entry_dir, "pdfresults_std.db")

    if ref_mode == "auto":
        if HAS_H5PY and os.path.isfile(h5_path):
            return read_proteindf_h5(h5_path), "h5"
        elif os.path.isfile(db_path):
            return read_proteindf_db(db_path), "db"
        else:
            entry_name = os.path.basename(entry_dir)
            raise FileNotFoundError(
                f"Neither reference.h5 nor pdfresults_std.db found in {entry_name}"
            )
    elif ref_mode == "h5":
        if not HAS_H5PY:
            raise RuntimeError("h5py is not installed but --ref h5 was requested")
        if not os.path.isfile(h5_path):
            entry_name = os.path.basename(entry_dir)
            raise FileNotFoundError(f"reference.h5 not found in {entry_name}")
        return read_proteindf_h5(h5_path), "h5"
    elif ref_mode == "db":
        if not os.path.isfile(db_path):
            entry_name = os.path.basename(entry_dir)
            raise FileNotFoundError(f"pdfresults_std.db not found in {entry_name}")
        return read_proteindf_db(db_path), "db"
    else:
        raise ValueError(f"Unknown ref_mode: {ref_mode}")


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

    mm_list = parsed_input.get("mm_charges", [])
    mm_self_energy = 0.0
    if mm_list:
        from pyscf import qmmm

        coords = [(x, y, z) for x, y, z, q in mm_list]
        charges = [q for x, y, z, q in mm_list]
        calc = qmmm.mm_charge(calc, coords, charges, unit=parsed_input["unit"])
        if len(mm_list) > 1:
            to_bohr = (
                1.0 / 0.529177210903
                if parsed_input["unit"].lower().startswith("ang")
                else 1.0
            )
            for i in range(len(mm_list)):
                for j in range(i + 1, len(mm_list)):
                    dx = (coords[i][0] - coords[j][0]) * to_bohr
                    dy = (coords[i][1] - coords[j][1]) * to_bohr
                    dz = (coords[i][2] - coords[j][2]) * to_bohr
                    r = math.sqrt(dx * dx + dy * dy + dz * dz)
                    if r > 1e-12:
                        mm_self_energy += (charges[i] * charges[j]) / r

    energy = calc.kernel()
    if mm_self_energy != 0.0:
        energy += mm_self_energy

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


def compare_entry(entry_dir, basis_file, db_override=None, ref_mode="auto",
                  xc_override=None, vwn_variant=None, grid_levels=(3,),
                  threshold=1e-4):
    """Compare a single test entry between ProteinDF and PySCF across grid levels."""
    userinput_path = os.path.join(entry_dir, "fl_Userinput")

    try:
        parsed_input = parse_fl_userinput(userinput_path)
    except Exception as e:
        return [{
            "entry": os.path.basename(entry_dir),
            "ref_format": "-",
            "method": "-",
            "xc_input": "-",
            "xc_used": "-",
            "grid_level": None,
            "pdf_nao": None,
            "pyscf_nao": None,
            "nao_match": False,
            "pdf_energy": None,
            "pyscf_energy": None,
            "diff": None,
            "converged": False,
            "spin": None,
            "s2": None,
            "mult": None,
            "n_alpha": None,
            "n_beta": None,
            "status": f"SKIP: {e}",
        }]

    # Reference data reading (HDF5 or DB)
    try:
        pdf_res, ref_format = read_proteindf_reference(
            entry_dir, ref_mode=ref_mode, db_override=db_override
        )
    except Exception as e:
        return [{
            "entry": os.path.basename(entry_dir),
            "ref_format": "-",
            "method": (parsed_input.get("method") or "-").upper(),
            "xc_input": parsed_input.get("xc") or "HF",
            "xc_used": "-",
            "grid_level": None,
            "pdf_nao": None,
            "pyscf_nao": None,
            "nao_match": False,
            "pdf_energy": None,
            "pyscf_energy": None,
            "diff": None,
            "converged": False,
            "spin": parsed_input.get("spin"),
            "s2": None,
            "mult": None,
            "n_alpha": parsed_input.get("n_alpha"),
            "n_beta": parsed_input.get("n_beta"),
            "status": f"SKIP: {e}",
        }]

    # Check for XC fitting (~ in XC functional)
    xc_in = parsed_input.get("xc") or ""
    xc_ref = pdf_res.get("xc_functional") or ""
    if "~" in xc_in or "~" in xc_ref:
        return [{
            "entry": os.path.basename(entry_dir),
            "ref_format": ref_format,
            "method": (parsed_input.get("method") or pdf_res.get("method") or "-").upper(),
            "xc_input": xc_in or "HF",
            "xc_used": "-",
            "grid_level": None,
            "pdf_nao": pdf_res.get("num_aos"),
            "pyscf_nao": None,
            "nao_match": False,
            "pdf_energy": pdf_res.get("total_energy"),
            "pyscf_energy": None,
            "diff": None,
            "converged": False,
            "spin": parsed_input.get("spin"),
            "s2": None,
            "mult": None,
            "n_alpha": parsed_input.get("n_alpha"),
            "n_beta": parsed_input.get("n_beta"),
            "status": "SKIP: XC fitting (~) not supported in PySCF",
        }]

    if parsed_input.get("method") == "roks":
        return [{
            "entry": os.path.basename(entry_dir),
            "ref_format": ref_format,
            "method": "ROKS",
            "xc_input": xc_in or "HF",
            "xc_used": "-",
            "grid_level": None,
            "pdf_nao": pdf_res.get("num_aos"),
            "pyscf_nao": None,
            "nao_match": False,
            "pdf_energy": pdf_res.get("total_energy"),
            "pyscf_energy": None,
            "diff": None,
            "converged": False,
            "spin": parsed_input.get("spin"),
            "s2": None,
            "mult": None,
            "n_alpha": parsed_input.get("n_alpha"),
            "n_beta": parsed_input.get("n_beta"),
            "status": "SKIP: ROKS not supported in PySCF comparison",
        }]

    results = []
    for g_lvl in grid_levels:
        pyscf_res = run_pyscf_calc(
            parsed_input,
            basis_file,
            xc_override=xc_override,
            vwn_variant=vwn_variant,
            grid_level=g_lvl,
        )

        diff = (
            pyscf_res["energy"] - pdf_res["total_energy"]
            if (pyscf_res["energy"] is not None and pdf_res["total_energy"] is not None)
            else None
        )
        if not pyscf_res["converged"]:
            status = "NOT_CONV"
        elif diff is not None and abs(diff) > threshold:
            status = "FAIL"
        else:
            status = "OK"

        results.append({
            "entry": os.path.basename(entry_dir),
            "ref_format": ref_format,
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
            "status": status,
        })

    return results


def check_vwn_variants(suite_dir, basis_file, output_format="table"):
    """Evaluate candidate VWN variants for Ne_SVWN5 and Ne_B3LYP."""
    print("=== VWN Variant Analysis for Ne_SVWN5 and Ne_B3LYP ===")

    # 1. Ne_SVWN5
    ne_svwn5_dir = os.path.join(suite_dir, "Ne_SVWN5")
    parsed_svwn5 = parse_fl_userinput(os.path.join(ne_svwn5_dir, "fl_Userinput"))
    pdf_svwn5, _ = read_proteindf_reference(ne_svwn5_dir, ref_mode="auto")

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
    pdf_b3lyp, _ = read_proteindf_reference(ne_b3lyp_dir, ref_mode="auto")

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
            print(
                "| Entry | Ref | Method | 2S (α, β) | <S^2> | Mult (2S+1) | "
                "PySCF Energy (Eh) | PDF Energy (Eh) | Diff (Eh) | Status |"
            )
            print("|:---|:---:|:---:|:---:|:---:|:---:|:---|:---|:---|:---|")
            for r in results_list:
                ref_str = r.get("ref_format", "-")
                spin_str = (
                    f"{r['spin']} ({r['n_alpha']}, {r['n_beta']})"
                    if r.get("spin") is not None and r.get("n_alpha") is not None
                    else "-"
                )
                s2_str = f"{r['s2']:.4f}" if r.get("s2") is not None else "-"
                mult_str = f"{r['mult']:.4f}" if r.get("mult") is not None else "-"
                pyscf_e_str = (
                    f"{r['pyscf_energy']:.10f}"
                    if r.get("pyscf_energy") is not None
                    else "-"
                )
                pdf_e_str = (
                    f"{r['pdf_energy']:.10f}"
                    if r.get("pdf_energy") is not None
                    else "-"
                )
                diff_str = (
                    f"{r['diff']:+.10e}"
                    if r.get("diff") is not None
                    else "-"
                )
                status_str = r.get("status", "OK")
                print(
                    f"| {r['entry']} | {ref_str} | {r['method']} | {spin_str} | "
                    f"{s2_str} | {mult_str} | {pyscf_e_str} | {pdf_e_str} | {diff_str} | {status_str} |"
                )
        else:
            print(
                "| Entry | Ref | Method | XC (PySCF) | Grid | AO (PDF/PySCF) | "
                "PDF Energy (Eh) | PySCF Energy (Eh) | Diff (Eh) | Status |"
            )
            print("|:---|:---:|:---|:---|:---:|:---:|:---|:---|:---|:---|")
            for r in results_list:
                ref_str = r.get("ref_format", "-")
                grid_str = (
                    str(r["grid_level"])
                    if r.get("grid_level") is not None
                    else "-"
                )
                if r.get("pdf_nao") is not None and r.get("pyscf_nao") is not None:
                    ao_str = f"{r['pdf_nao']}/{r['pyscf_nao']}"
                elif r.get("pdf_nao") is not None:
                    ao_str = f"{r['pdf_nao']}/-"
                else:
                    ao_str = "-"
                pdf_e_str = (
                    f"{r['pdf_energy']:.10f}"
                    if r.get("pdf_energy") is not None
                    else "-"
                )
                pyscf_e_str = (
                    f"{r['pyscf_energy']:.10f}"
                    if r.get("pyscf_energy") is not None
                    else "-"
                )
                diff_str = (
                    f"{r['diff']:+.10e}"
                    if r.get("diff") is not None
                    else "-"
                )
                status_str = r.get("status", "OK")
                print(
                    f"| {r['entry']} | {ref_str} | {r['method']} | {r['xc_used']} | "
                    f"{grid_str} | {ao_str} | {pdf_e_str} | {pyscf_e_str} | {diff_str} | {status_str} |"
                )
    else:
        if show_spin:
            header = (
                f"{'Entry':25s} {'Ref':5s} {'Method':6s} {'2S (α, β)':12s} {'<S^2>':7s} "
                f"{'Mult':7s} {'PySCF Energy (Eh)':18s} {'PDF Energy (Eh)':18s} {'Diff (Eh)':16s} {'Status':s}"
            )
            print(header)
            print("-" * len(header))
            for r in results_list:
                ref_str = r.get("ref_format", "-")
                spin_str = (
                    f"{r['spin']} ({r['n_alpha']}, {r['n_beta']})"
                    if r.get("spin") is not None and r.get("n_alpha") is not None
                    else "-"
                )
                s2_str = f"{r['s2']:.4f}" if r.get("s2") is not None else "-"
                mult_str = f"{r['mult']:.4f}" if r.get("mult") is not None else "-"
                pyscf_e_str = (
                    f"{r['pyscf_energy']:18.10f}"
                    if r.get("pyscf_energy") is not None
                    else f"{'-':>18s}"
                )
                pdf_e_str = (
                    f"{r['pdf_energy']:18.10f}"
                    if r.get("pdf_energy") is not None
                    else f"{'-':>18s}"
                )
                diff_str = (
                    f"{r['diff']:+16.10e}"
                    if r.get("diff") is not None
                    else f"{'-':>16s}"
                )
                status_str = r.get("status", "OK")
                print(
                    f"{r['entry']:25s} {ref_str:5s} {r['method']:6s} {spin_str:12s} "
                    f"{s2_str:7s} {mult_str:7s} {pyscf_e_str} {pdf_e_str} {diff_str} {status_str}"
                )
        else:
            header = (
                f"{'Entry':25s} {'Ref':5s} {'Method':6s} {'XC (PySCF)':22s} {'Grid':4s} "
                f"{'AO(PDF/Py)':10s} {'PDF Energy (Eh)':18s} {'PySCF Energy (Eh)':18s} {'Diff (Eh)':16s} {'Status':s}"
            )
            print(header)
            print("-" * len(header))
            for r in results_list:
                ref_str = r.get("ref_format", "-")
                grid_str = (
                    str(r["grid_level"])
                    if r.get("grid_level") is not None
                    else "-"
                )
                if r.get("pdf_nao") is not None and r.get("pyscf_nao") is not None:
                    ao_str = f"{r['pdf_nao']}/{r['pyscf_nao']}"
                elif r.get("pdf_nao") is not None:
                    ao_str = f"{r['pdf_nao']}/-"
                else:
                    ao_str = "-"
                pdf_e_str = (
                    f"{r['pdf_energy']:18.10f}"
                    if r.get("pdf_energy") is not None
                    else f"{'-':>18s}"
                )
                pyscf_e_str = (
                    f"{r['pyscf_energy']:18.10f}"
                    if r.get("pyscf_energy") is not None
                    else f"{'-':>18s}"
                )
                diff_str = (
                    f"{r['diff']:+16.10e}"
                    if r.get("diff") is not None
                    else f"{'-':>16s}"
                )
                status_str = r.get("status", "OK")
                print(
                    f"{r['entry']:25s} {ref_str:5s} {r['method']:6s} {r['xc_used'][:22]:22s} "
                    f"{grid_str:4s} {ao_str:10s} {pdf_e_str} {pyscf_e_str} {diff_str} {status_str}"
                )


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
    parser.add_argument(
        "--ref",
        choices=["auto", "h5", "db"],
        default="auto",
        help="Reference data format (auto: prefer h5 if available, else db; h5: require h5; db: require db)",
    )
    parser.add_argument("--db", help="Path to specific calculation DB (default: pdfresults_std.db in entry)")
    parser.add_argument("--xc", help="Override XC functional string for PySCF")
    parser.add_argument("--vwn", default=None, help="VWN variant for SVWN/B3LYP (options: VWN, VWN_RPA, VWN_3, etc.)")
    parser.add_argument("--grid-level", default="3", help="Comma-separated PySCF grid levels (default: 3)")
    parser.add_argument(
        "--threshold",
        type=float,
        default=1e-4,
        help="Energy difference threshold in Hartree for FAIL status (default: 1e-4)",
    )
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
                ref_mode=args.ref,
                xc_override=args.xc,
                vwn_variant=args.vwn,
                grid_levels=grid_levels,
                threshold=args.threshold,
            )
            all_results.extend(res)
        except Exception as e:
            print(f"Error processing entry {entry}: {e}", file=sys.stderr)

    if all_results:
        print_results_table(all_results, output_format=args.format, show_spin=args.show_spin)

    has_failure = any(r.get("status") in ("FAIL", "NOT_CONV") for r in all_results)
    if has_failure:
        sys.exit(1)


if __name__ == "__main__":
    main()
