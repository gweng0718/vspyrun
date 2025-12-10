import re
from collections import OrderedDict
from ase.io import read as ase_read  # if you're using ASE


def parse_poscar_with_ase(poscar_path):
    atoms = ase_read(poscar_path, format="vasp")
    symbols = atoms.get_chemical_symbols()

    counts_ordered = OrderedDict()
    for s in symbols:
        counts_ordered[s] = counts_ordered.get(s, 0) + 1

    elements = list(counts_ordered.keys())
    counts = list(counts_ordered.values())
    return elements, counts


def parse_potcar(potcar_path):
    zvals = []
    with open(potcar_path, "r") as f:
        for line in f:
            if "ZVAL" in line:
                m = re.search(r"ZVAL\s*=\s*([-0-9.]+)", line)
                if m:
                    zvals.append(float(m.group(1)))

    if not zvals:
        raise RuntimeError(f"No ZVAL entries found in POTCAR: {potcar_path}")

    return zvals


def get_vale(poscar="POSCAR", potcar="POTCAR"):
    """
    Compute total valence electrons from POSCAR + POTCAR.
    """
    elements, counts = parse_poscar_with_ase(poscar)
    zvals = parse_potcar(potcar)

    if len(elements) != len(zvals):
        raise ValueError(
            f"Number of elements in POSCAR ({len(elements)}) "
            f"does not match number of ZVAL entries in POTCAR ({len(zvals)})."
        )

    total_val = sum(c * z for c, z in zip(counts, zvals))
    return total_val
