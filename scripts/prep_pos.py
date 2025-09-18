#!/usr/bin/env python3
"""
Minimal SCF input generator using your vspyrun package.
"""

import os
import sys

# POSCAR ops
vacuum = 15.0
n_top  = 8
n_bot  = 8
# ----------------------------------------------------------------------

# Import from your package
from vspyrun import add_vacuum, center_atoms, write_pos_cart, apply_sd, read_pos, keep_topn

def main(inpos_file, run_dir):

    # (1) POSCAR: read, add vacuum, center, apply selective dynamics, write out
    atoms = read_pos(inpos_file)
    atoms = keep_topn(atoms, N)
    outpos_file = os.path.join(run_dir, "POSCAR")
    write_pos_cart(atoms, outpos_file)

    print(f"POSCAR generated in {run_dir}")

if __name__ == "__main__":

    #prepare vasp_dir
    base_dir = os.getcwd()
    file_dir = os.path.join(base_dir, "files")
    inpos_file = os.path.join(file_dir, "input.vasp")
    if not os.path.isfile(inpos_file):
        print(f"ERROR: INPUT POSCAR not found at '{inpos_file}'. Set input path at top of this script.")
        sys.exit(1)

    run_dir = os.path.join(base_dir, "test")
    os.makedirs(run_dir, exist_ok=True)
    main(inpos_file, run_dir)
