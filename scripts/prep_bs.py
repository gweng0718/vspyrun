#!/usr/bin/env python3
"""
Minimal SCF input generator using your vspyrun package.
"""

import os
import sys
import shutil
import subprocess

# --- CONFIG ------------------------------------------------------------
# INCAR (electronic settings) — keep as strings to match your function signatures
#Control
LWAVE = ".FALSE."
# electron scf
ALGO  = "Normal"
AMIN  = "0.01"
ENCUT = "500"
PREC  = "Accurate"
SIGMA = "0.1"
ISPIN = "1"
XCF   = "PE"        # FIX: was undefined (you had GGA="PE"); vspyrun.incar_elec expects xcf
# spin-orbit
LSORBIT = ".TRUE."
# energy-correction
IVDW = "20"
# ionic relaxation and optimization
IBRION = "-1"
NSW    = "0"
# parallelization
KPAR  = "8"
NPAR  = "4"
NCORE = "1"
# POSCAR ops
vacuum = 15.0
n_top  = 8
n_bot  = 8
# ----------------------------------------------------------------------

# Import from your package
from vspyrun import incar_sys_bs, incar_elec, incar_soc, incar_ecorr, incar_geo, incar_para, join_incar_blocks
from vspyrun import generate_potcar
from vspyrun import add_vacuum, center_atoms, write_pos_cart, apply_sd, read_pos
from vspyrun import write_file

def main(inpos_file, inkpts_file, chgcar_in, run_dir, pot_dir):

    # (1) POSCAR: read, add vacuum, center, apply selective dynamics, write out
    atoms = read_pos(inpos_file)
    outpos_file = os.path.join(run_dir, "POSCAR")
    shutil.copy(inpos_file,outpos_file)
    nat = len(atoms)
    ndim = 3*nat
    MAGMOM = f"{ndim}*0.0"

    # (2) Preparing INCAR
    incar_text = join_incar_blocks(
        incar_sys_bs(lwave=LWAVE),
        incar_elec(algo=ALGO, amin=AMIN, encut=ENCUT, prec=PREC, sigma=SIGMA, ispin=ISPIN, xcf=XCF),
        incar_soc(lsoc=LSORBIT, magmom=MAGMOM),                 # FIX: add missing commas
        incar_geo(ibrion=IBRION, nsw=NSW),
        incar_para(kpar=KPAR, npar=NPAR, ncore=NCORE),
    )
    incar_text = "\n".join([ln for ln in incar_text.splitlines() if ln.strip()]) + "\n"
    incar_file = os.path.join(run_dir, "INCAR")
    write_file(incar_file, incar_text)

    # (3) KPOINTS: copy from preexisting files
    kpts_file = os.path.join(run_dir, "KPOINTS")
    shutil.copy(inkpts_file, kpts_file)

    # (4) POTCAR: assembled from element folders under POTCAR_DIR
    outpot_file = os.path.join(run_dir, "POTCAR")
    generate_potcar(inpos=outpos_file, potcar_dir=pot_dir, outpot=outpot_file)

    # (5) CHGCAR: create a link to CHGCAR file from SCF
    chgcar_out = os.path.join(run_dir, "CHGCAR")
    os.symlink(chgcar_in, chgcar_out)

    print(f"VASP inputs generated in {run_dir}.")

if __name__ == "__main__":

    pot_dir = ""
    if not os.path.isdir(pot_dir):
        print(f"ERROR: POTCAR_DIR does not exist: '{pot_dir}'. Set the correct path.")
        sys.exit(1)

    #prepare vasp_dir
    base_dir = os.getcwd()
    scf_dir = ""
    file_dir = os.path.join(base_dir, "files")
    for j in range(1,4):

        runj_dir = os.path.join(base_dir, f"run{j}")
        os.makedirs(runj_dir, exist_ok=True)
        scfj_dir = os.path.join(scf_dir, f"run{j}")
        inpos_file = os.path.join(scfj_dir, "POSCAR")
        if not os.path.isfile(inpos_file):
            print(f"ERROR: INPUT POSCAR not found at '{inpos_file}'. Set input path at top of this script.")
            continue

        chgcar_in = os.path.join(scfj_dir, "CHGCAR")
        if not os.path.isfile(chgcar_in):
            print(f"{chgcar_in} not found!")
            continue

        for i in range(1,4):

            runi_dir = os.path.join(runj_dir, f"run{i}")
            os.makedirs(runi_dir, exist_ok=True)
            kpts_file = os.path.join(file_dir, f"K{i}")
            if not os.path.isfile(kpts_file):
                print(f"ERROR: INPUT KPOINTS not found at '{kpts_file}'. Set input path at top of this script.")
                continue

            main(inpos_file, kpts_file, chgcar_in, runi_dir, pot_dir)

            #submit jobs
            from vspyrun import prep_slurm
            template = os.path.join(file_dir, "slurm.job")
            jname = f"AZcr{j}K{i}"
            outslurm = os.path.join(runi_dir, "slurm.job")
            prep_slurm(template, jname, outslurm)
            os.chdir(runi_dir)
            subprocess.run(['sbatch', 'slurm.job'])
            print(f"Successfully submit vasp job in {runi_dir}!", flush=True)
