#!/usr/bin/env python3
"""
Minimal SCF input generator using your vspyrun package.
"""

import os
import sys

# --- CONFIG ------------------------------------------------------------
POTCAR_DIR = "/Users/gweng/codes/vasp/potpaw_PBE"

# INCAR (electronic settings) — keep as strings to match your function signatures
# electron scf
ALGO  = "Normal"
AMIN  = "0.01"
ENCUT = "500"
PREC  = "Normal"
SIGMA = "0.1"
ISPIN = "1"
XCF   = "PE"        # FIX: was undefined (you had GGA="PE"); vspyrun.incar_elec expects xcf
# spin-orbit
LSORBIT = ".FALSE."
# energy-correction
IVDW   = "20"
LDIPOL = ".FALSE."
IDIPOL = "-1"
# ionic relaxation and optimization
IBRION = "1"
NSW    = "100"
EDIFFG = "5E-3"
ISIF   = 2
# parallelization (GPU)
KPAR   = "4"
NPAR   = "8"
NCORE  = "1"
# POSCAR ops
vacuum = 15.0
n_top  = 8
n_bot  = 8
# KPOINTS
K_CENTER = "Gamma"   # "Gamma" or "Monkhorst-Pack"
K_MESH   = (8, 8, 8) # k-point mesh for SCF
# ----------------------------------------------------------------------

# Import from your package
from vspyrun import incar_sys_rel, incar_elec, incar_soc, incar_ecorr, incar_geo, incar_para, join_incar_blocks
from vspyrun import kpts_scf
from vspyrun import generate_potcar
from vspyrun import add_vacuum, center_atoms, write_pos_cart, apply_sd, read_pos
from vspyrun import write_file

def main(inpos_file, run_dir, pot_dir):

    # (1) POSCAR: read, add vacuum, center, apply selective dynamics, write out
    atoms = read_pos(inpos_file)
    atoms = add_vacuum(atoms, vacuum)
    atoms = center_atoms(atoms)
    atoms = apply_sd(atoms, n_bot, n_top)
    outpos_file = os.path.join(run_dir, "POSCAR")
    write_pos_cart(atoms, outpos_file)
    nat = len(atoms)
    ndim = 3*nat
    MAGMOM = f"{ndim}*0.0"

    # (2) Preparing INCAR
    incar_text = join_incar_blocks(
        incar_sys_rel(),
        incar_elec(algo=ALGO, amin=AMIN, encut=ENCUT, prec=PREC, sigma=SIGMA, ispin=ISPIN, xcf=XCF),
        incar_soc(lsoc=LSORBIT,magmom=MAGMOM),                 # FIX: add missing commas
        incar_ecorr(ivdw=IVDW, ldipol=LDIPOL, idipol=IDIPOL),
        incar_geo(ibrion=IBRION, nsw=NSW, ediffg=EDIFFG, isif=ISIF),
        incar_para(kpar=KPAR, npar=NPAR, ncore=NCORE),
    )
    incar_text = "\n".join([ln for ln in incar_text.splitlines() if ln.strip()]) + "\n"
    incar_file = os.path.join(run_dir, "INCAR")
    write_file(incar_file, incar_text)

    # (3) KPOINTS: uniform mesh
    nx, ny, nz = K_MESH
    kpoints_text = kpts_scf(kcen=K_CENTER, nx=str(nx), ny=str(ny), nz=str(nz)).lstrip("\n")
    kpts_file = os.path.join(run_dir, "KPOINTS")
    write_file(kpts_file, kpoints_text)

    # (4) POTCAR: assembled from element folders under POTCAR_DIR
    outpot_file = os.path.join(run_dir, "POTCAR")
    generate_potcar(inpos=outpos_file, potcar_dir=pot_dir, outpot=outpot_file)

    print(f"VASP inputs generated in {run_dir}")

if __name__ == "__main__":

    pot_dir = "/Users/gweng/codes/vasp/potpaw_PBE"
    if not os.path.isdir(pot_dir):
        print(f"ERROR: POTCAR_DIR does not exist: '{pot_dir}'. Set the correct path.")
        sys.exit(1)

    #prepare vasp_dir
    base_dir = os.getcwd()
    file_dir = os.path.join(base_dir, "files")
    inpos_file = os.path.join(file_dir, "input.vasp")
    if not os.path.isfile(inpos_file):
        print(f"ERROR: INPUT POSCAR not found at '{inpos_file}'. Set input path at top of this script.")
        sys.exit(1)

    run_dir = os.path.join(base_dir, "test")
    os.makedirs(run_dir, exist_ok=True)
    main(inpos_file, run_dir, pot_dir)

    #submit jobs
    from vspyrun import prep_slurm
    template = os.path.join(file_dir, "slurm.job")
    jname = "test"
    outslurm = os.path.join(run_dir, "slrum.job")
    prep_slurm(template, jname, outslurm)
    os.chdir(run_dir)
    #subprocess.run(['sbatch', 'slurm.job'])
    print(f"Successfully submit vasp job in {run_dir}!", flush=True)
