#!/usr/bin/env python3
"""
Minimal SCF input generator using your vspyrun package.
"""

import os
import sys
import subprocess

# --- CONFIG ------------------------------------------------------------
# INCAR (electronic settings) — keep as strings to match your function signatures
#system
SYSTEM = "AZCO2"
NWRITE = "1"
ISTART = "0"
ICHARG = "2"
LWAVE = ".FALSE."
LCHARG = ".TRUE."
LREAL = "AUTO"
LORBIT = "0"
ISYM = "2"
# electron scf
ALGO  = "Normal"
AMIN  = "0.1"
ENCUT = "400"
NELM  = "200"
PREC  = "Normal"
SIGMA = "0.1"
ISPIN = "1"
XCF   = "PE"        # FIX: was undefined (you had GGA="PE"); vspyrun.incar_elec expects xcf
LKILLSCF = ".TRUE."
# spin-orbit
LSORBIT = ".FALSE."
# energy-correction
IVDW = "20"
LDIPOL = ".FALSE."
IDIPOL = "-1"
# ionic relaxation and optimization
IBRION = "-1"
NSW    = "0"
EDIFFG = "1E-3"
# parallelization
KPAR  = "4"
NPAR  = "16"
NCORE = "1"
# POSCAR ops
vacuum = 15.0
n_top  = 29
n_bot  = 0
# KPOINTS
K_CENTER = "Gamma"   # "Gamma" or "Monkhorst-Pack"
K_MESH   = (4, 4, 1) # k-point mesh for SCF
# ----------------------------------------------------------------------

# Import from your package
from vspyrun import incar_sys, incar_elec, incar_soc, incar_ecorr, incar_geo, incar_para, incar_solv, get_vale, join_incar_blocks
from vspyrun import kpts_scf
from vspyrun import generate_potcar
from vspyrun import add_vacuum, center_atoms, write_pos_cart, apply_sd, read_pos, build_sc, keep_topn, shift_min_z, build_mir_sym
from vspyrun import write_file

def main(inpos_file, run_dir, pot_dir, charge):

    # (1) POSCAR: read, add vacuum, center, apply selective dynamics, write out
    atoms = read_pos(inpos_file)
    atoms = keep_topn(atoms, 92)
    atoms = shift_min_z(atoms, 1.0)
    atoms = add_vacuum(atoms, vacuum)
    atoms = apply_sd(atoms, n_bot, n_top)
    outpos_file = os.path.join(run_dir, "POSCAR")
    write_pos_cart(atoms, outpos_file)
    nat = len(atoms)
    ndim = 3*nat
    MAGMOM = f"{ndim}*0.0"

    # (2) KPOINTS: uniform mesh
    nx, ny, nz = K_MESH
    kpoints_text = kpts_scf(kcen=K_CENTER, nx=str(nx), ny=str(ny), nz=str(nz)).lstrip("\n")
    kpts_file = os.path.join(run_dir, "KPOINTS")
    write_file(kpts_file, kpoints_text)

    # (3) POTCAR: assembled from element folders under POTCAR_DIR
    outpot_file = os.path.join(run_dir, "POTCAR")
    generate_potcar(inpos=outpos_file, potcar_dir=pot_dir, outpot=outpot_file)

    # (4) Preparing INCAR
    #print(f"NELECT for VASPsol: {NELECT}")
    incar_text = join_incar_blocks(
        incar_sys(system=SYSTEM,nwrite=NWRITE,istart=ISTART,icharg=ICHARG,lwave=LWAVE,lcharg=LCHARG,lreal=LREAL,lorbit=LORBIT,isym=ISYM),    
        incar_elec(algo=ALGO, amin=AMIN, encut=ENCUT, prec=PREC, sigma=SIGMA, ispin=ISPIN, xcf=XCF, lkillscf=LKILLSCF, nelm=NELM),
        incar_soc(lsoc=LSORBIT, magmom=MAGMOM),                 # FIX: add missing commas
        incar_ecorr(ivdw=IVDW, ldipol=LDIPOL, idipol=IDIPOL),
        incar_geo(ibrion=IBRION, nsw=NSW, ediffg=EDIFFG),
        incar_para(kpar=KPAR, npar=NPAR, ncore=NCORE)
    )
    incar_text = "\n".join([ln for ln in incar_text.splitlines() if ln.strip()]) + "\n"
    incar_file = os.path.join(run_dir, "INCAR")
    write_file(incar_file, incar_text)

    print(f"VASP inputs generated in {run_dir}")

if __name__ == "__main__":

    pot_dir = "/lustre/grweng/codes/pseudo/potpaw_PBE"
    if not os.path.isdir(pot_dir):
        print(f"ERROR: POTCAR_DIR does not exist: '{POTCAR_DIR}'. Set the correct path.")
        sys.exit(1)

    #prepare vasp_dir
    base_dir = os.getcwd()
    file_dir = os.path.join(base_dir, "files")
    struc_dir = os.path.join(base_dir, "input_struc")
    for i in range(1,4):
        runi_dir = os.path.join(base_dir,f"run{i}")
        inpos_file = os.path.join(struc_dir, f"r{i}.vasp")
        if not os.path.isfile(inpos_file):
            print(f"ERROR: INPUT POSCAR not found at '{inpos_file}'. Set input path at top of this script.")
            continue
        os.makedirs(runi_dir, exist_ok=True)
        ncharge = 0 
        main(inpos_file, runi_dir, pot_dir, ncharge)
        #submit jobs
        from vspyrun import prep_slurm
        template = os.path.join(file_dir, "slurm.job")
        jname = f"AZCO1{i}"
        outslurm = os.path.join(runi_dir, "slurm.job")
        prep_slurm(template, jname, outslurm)
        os.chdir(runi_dir)
        #subprocess.run(['sbatch', 'slurm.job'])
        print(f"Successfully submit vasp job in {runi_dir}!", flush=True)
