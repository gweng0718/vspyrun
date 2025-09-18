# incar.py — debugged & lightly optimized
# - FIX: added missing colon on def incar_sys_rel
# - ADD: helper functions join_incar_blocks(), write_incar(), bool_flag() for composition
# No behavioral changes to your existing block generators.

import os
import itertools
import shutil
import subprocess
import numpy as np

# -------------------------------
# Utilities to compose INCAR text
# -------------------------------
import textwrap

def _strip_blank_edges(s: str) -> str:
    """Normalize block indentation and ensure trailing newline."""
    s = textwrap.dedent(s).strip("\n")
    return s + "\n"

def join_incar_blocks(*blocks: str) -> str:
    """Join multiple INCAR sections with a single blank line between them."""
    cleaned = [_strip_blank_edges(b) for b in blocks if b and b.strip()]
    return ("\n".join(cleaned) + ("\n" if cleaned else ""))

def write_incar(path: str, *blocks: str) -> None:
    """Write merged INCAR at `path` from the given blocks."""
    text = join_incar_blocks(*blocks)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

def bool_flag(value: bool | str) -> str:
    """Convert bool or string to VASP boolean (.TRUE./.FALSE.)."""
    if isinstance(value, str):
        v = value.strip().upper()
        if v in {".TRUE.", "TRUE", "T"}:  return ".TRUE."
        if v in {".FALSE.", "FALSE", "F"}: return ".FALSE."
        raise ValueError(f"Unrecognized boolean string: {value!r}")
    return ".TRUE." if value else ".FALSE."


def append_unique_line(filepath: str, line: str) -> bool:
    """Append `line` to `filepath` if it's not already present (ignoring whitespace).
       Returns True if appended, False if skipped."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            existing = [L.strip() for L in f]
    except FileNotFoundError:
        existing = []

    if line.strip() in existing:
        return False  # already there

    needs_nl = False
    if os.path.exists(filepath):
        with open(filepath, "rb") as f:
            data = f.read()
        needs_nl = (len(data) > 0 and not data.endswith(b"\n"))

    with open(filepath, "a", encoding="utf-8") as f:
        if needs_nl:
            f.write("\n")
        f.write(line.rstrip("\n") + "\n")
    return True

def incar_sys_scf(*, icharg="2", lcharg=".TRUE."):
    content = f"""
#system description
SYSTEM = 'scf'
NWRITE = 2
ISTART = 0
ICHARG = {icharg}
LWAVE= .FALSE.
LCHARG = {lcharg}
LREAL = AUTO
LORBIT = 0
ISYM = 2
"""
    return content

def incar_sys_bs(*, lwave=".FALSE."):
    content = f"""
#system description
SYSTEM = 'bands'
NWRITE = 2
ISTART = 0
ICHARG = 11
LWAVE= {lwave}
LCHARG = .FALSE.
LREAL = AUTO
LORBIT = 11
ISYM = 0
"""

    return content

def incar_sys_rel():
    content = f"""
#system description
SYSTEM = 'relax'
NWRITE = 2
ISTART = 0
ICHARG = 2
LWAVE= .FALSE.
LCHARG = .FALSE.
LREAL = AUTO
LORBIT = 0
ISYM = 2
"""
    return content

def incar_elec(*, algo="Normal", amin="0.1", encut="500", prec="Accurate", sigma="0.1", ispin="1", xcf="PE"):
    content = f"""
#electron scf
ALGO = {algo}
AMIN = {amin}
ENCUT = {encut}
PREC = {prec}
ISMEAR = 0
SIGMA = {sigma}
ISPIN = {ispin}
NELM = 500
EDIFF = 1E-7
GGA = {xcf}
LASPH = .TRUE.
"""
    return content

def incar_soc(*, lsoc=".TRUE.", magmom="100*0.0"):
    content = f"""
#spin-orbit
LSORBIT = {lsoc}
MAGMOM = {magmom}
"""
    return content

def incar_ecorr(*, ivdw="20", ldipol=".FALSE.", idipol="-1"):
    content = f"""
#energy-correction
IVDW = {ivdw}
LDIPOL = {ldipol}
IDIPOL = {idipol}
"""
    return content

def incar_geo(*, ibrion="-1", ediffg="5E-3", nsw="0", isif="2"):
    content = f"""
#ionic relaxation and optimization
IBRION = {ibrion}
EDIFFG = {ediffg}
NSW = {nsw}
ISIF = {isif}
NFREE = 1
POTIM = 0.5
"""
    return content

def incar_para(*, kpar="4", npar="8", ncore="8"):
    content = f"""
#parallelization related
KPAR = {kpar}
NPAR = {npar}
NCORE = {ncore}
NSIM = 4
"""
    return content

def incar_solv(*, lsol=".TRUE.", eb_k="78.4", tau="0.0", lambda_k="3.0", sig_k="0.4", lrhoion=".FALSE."):
    content = f"""
#VASPSOL
LSOL    = {lsol}
EB_K    = {eb_k}
TAU     = {tau}
LRHOION = {lrhoion}
LAMBDA_D_K = {lambda_k}
SIGMA_K    = {sig_k}
"""
    return content
