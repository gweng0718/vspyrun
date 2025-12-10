import sys
import numpy as np
import os
from ase.io import read, write
from ase import Atoms, Atom
from ase.constraints import FixAtoms

def read_pos(pos_file):
    atoms = read(pos_file, format="vasp")
    return atoms

def keep_topn(atoms: Atoms, n: int, *, use_scaled: bool = False) -> Atoms:
    """
    Return a new Atoms object containing only the top N atoms by z.

    Parameters
    ----------
    atoms : ase.Atoms
        Input structure.
    n : int
        Number of atoms to keep (highest z). If n <= 0, returns empty Atoms.
        If n >= len(atoms), returns a copy of the original.
    use_scaled : bool
        If True, rank by fractional z (useful if the slab wraps across the cell).
        Otherwise, rank by Cartesian z.

    Returns
    -------
    ase.Atoms
        Subset with the top N atoms, preserving original order.
    """
    total = len(atoms)
    if n <= 0:
        return Atoms(cell=atoms.get_cell(), pbc=atoms.get_pbc())
    if n >= total:
        sub = atoms.copy()
        sub.set_cell(atoms.get_cell())
        sub.set_pbc(atoms.get_pbc())
        return sub

    z = atoms.get_scaled_positions()[:, 2] if use_scaled else atoms.positions[:, 2]
    order = np.argsort(z)              # low -> high
    top_idx = order[-n:]               # indices of top N (unsorted)
    # preserve the original input order, like your rem_at_z does
    top_idx = np.sort(top_idx)

    sub = atoms[top_idx]               # slicing keeps arrays/constraints per-atom
    sub.set_cell(atoms.get_cell())     # ensure same cell/pbc
    sub.set_pbc(atoms.get_pbc())
    return sub

def rem_at_z(atoms, z_mirror, tol=1e-6):
    """
    Returns a new Atoms object that includes only those atoms with
    z >= z_mirror (within a small tolerance).
    """
    filtered_atoms = Atoms(cell=atoms.get_cell(), pbc=atoms.get_pbc())

    for atom in atoms:
        if atom.position[2] + tol >= z_mirror:
            # Keep if z >= z_mirror - tol
            filtered_atoms.append(atom)

    return filtered_atoms

def shift_min_z(atoms, z):
    z_coords = atoms.positions[:, 2]
    shift = np.min(z_coords)
    atoms.positions[:, 2] -= shift
    atoms.positions[:, 2] += z
  
    return atoms

def keep_at_el(atoms,el_list):
    """
    Returns a new Atoms object that includes only the atoms in the el_list.
    """
    mask = [sym in el_list for sym in atoms.get_chemical_symbols()]
    new_atoms = atoms[mask]
    
    return new_atoms

def add_ads(atoms1, atoms2, dz=0):
    """     
    Add adsorbates on the surface at a certain distance dz.
    """ 
    z_cor1 = atoms1.positions[:, 2]
    z_cor2 = atoms2.positions[:, 2]
    z_max = np.max(z_cor1)
    z_min = np.min(z_cor2)
    z_shift = z_max - z_min + dz
    atoms2.positions[:,2] += z_shift
    atoms_all = atoms1 + atoms2

    return atoms_all
    
def add_vacuum(atoms, vacuum_space):
    """
    Add vacuum_space along the c-axis
    """
    z_coords = atoms.positions[:, 2]
    z_min = np.min(z_coords)
    z_max = np.max(z_coords)
    thickness = z_max - z_min

    cell = atoms.get_cell().copy()
    new_c = thickness + vacuum_space
    cell[2, 2] = new_c
    atoms.set_cell(cell, scale_atoms=False)

    return atoms

def center_atoms(atoms):
    """
    shift the atoms to the center of the cell
    """
    cell = atoms.get_cell().copy()
    midpoint = 0.5 * cell[2, 2]
    z_coords = atoms.positions[:, 2]
    z_min = np.min(z_coords)
    z_max = np.max(z_coords)
    thickness = z_max - z_min
    shift = midpoint - (z_min + 0.5 * thickness)
    atoms.positions[:, 2] += shift

    return atoms

def build_mir_sym(atoms, z_mir, tol=0.1):
    """
    Reflect (mirror) the slab across z = z_mir.
    - If an atom's z is within 'tol' of z_mirror, skip it in the mirrored copy
      (so it is not duplicated).
    - Otherwise, z_new = 2*z_mir - z_old.
    x, y remain unchanged.
    """
    mir_atoms = Atoms(cell=atoms.get_cell(), pbc=atoms.get_pbc())

    for atom in atoms:
        z_old = atom.position[2]
        # Skip atoms that lie exactly on the mirror plane
        if abs(z_old - z_mir) <= tol:
            continue

        x_new = atom.position[0]
        y_new = atom.position[1]
        z_new = 2.0 * z_mir - z_old

        mir_atoms.append(
            Atom(symbol=atom.symbol, position=(x_new, y_new, z_new))
        )

    mir_atoms = atoms + mir_atoms

    return mir_atoms

def build_inver_sym(atoms, inv_center, tol_overlap=1e-3):
    """
    Apply inversion symmetry to all atoms, but drop any inverted atom whose site
    already exists in the original structure (within tol Å, using minimum-image).
    Returns only the non-duplicate inverted atoms (so you can do: atoms + returned).
    """
    # Prepare original geometry (for duplicate checks)
    cell = atoms.get_cell().array
    inv_cell = np.linalg.inv(cell)
    orig_cart = atoms.get_positions()
    orig_frac = orig_cart @ inv_cell  # fractional coords of originals

    new_atoms = atoms.copy()
    to_delete = []

    for i in range(len(new_atoms)):
        new_pos = 2 * inv_center - new_atoms[i].position

        # Compare new_pos to all original atoms using minimum-image convention
        f_new = new_pos @ inv_cell
        df = orig_frac - f_new                    # fractional deltas
        df -= np.round(df)                        # wrap to [-0.5, 0.5)
        dr = df @ cell                            # back to Cartesian
        dmin = float(np.min(np.linalg.norm(dr, axis=1)))

        if dmin < tol_overlap:
            # A preexisting atom already occupies this site -> discard inverted atom
            to_delete.append(i)
        else:
            # Keep this inverted atom
            new_atoms[i].position = new_pos

    # Remove duplicates (delete from the end to keep indices valid)
    for idx in sorted(to_delete, reverse=True):
        del new_atoms[idx]

    return new_atoms+atoms

def rep_lat_para(atoms1,atoms2):
    """
    substitue the lattice parameters for atoms2
    """
    sd_flags = atoms2.arrays.get("selective_dynamics", None)
    cell1 = atoms1.get_cell().copy()
    cell2 = atoms2.get_cell().copy()
    cell2[0,:] = cell1[0,:] #substituting a
    cell2[1,:] = cell1[1,:] #substituting b
    atoms2.set_cell(cell2, scale_atoms=True)
    atoms2.pbc = True 

    if sd_flags is not None and "selective_dynamics" not in atoms2.arrays:
        atoms2.set_array("selective_dynamics", sd_flags)

    return atoms2

def apply_sd(atoms_like, n_bot: int, n_top: int, *, use_scaled: bool = False):
    """
    Make bottom n_bot and top n_top atoms relax (T T T),
    freeze all others (F F F). Works on ase.Atoms or (Atoms, ...) tuples.

    Parameters
    ----------
    atoms_like : ase.Atoms or tuple/list whose first element is ase.Atoms
    n_bot, n_top : ints (>= 0)
    use_scaled : if True, rank by fractional z (helpful if slab wraps in z)

    Returns
    -------
    ase.Atoms  (same object, mutated in-place)
    """
    # Accept (Atoms, ...) tuple as well as bare Atoms
    atoms = atoms_like[0] if isinstance(atoms_like, (tuple, list)) and hasattr(atoms_like[0], "positions") else atoms_like
    if not hasattr(atoms, "positions"):
        raise TypeError("apply_sd_relax_edges expects an ase.Atoms or a (Atoms, ...) tuple")

    if n_bot < 0 or n_top < 0:
        raise ValueError("n_bot and n_top must be non-negative")

    n = len(atoms)
    if n == 0:
        raise ValueError("Atoms has no atoms")

    # Rank by z
    z = atoms.get_scaled_positions()[:, 2] if use_scaled else atoms.positions[:, 2]
    order = np.argsort(z)

    # Sets that should be relaxed
    bot = set(order[:n_bot]) if n_bot else set()
    top = set(order[-n_top:]) if n_top else set()
    relax = np.array(sorted(bot | top), dtype=int)

    # Start with all frozen (False False False)
    sd = np.zeros((n, 3), dtype=bool)
    if relax.size:
        sd[relax, :] = True  # edges relax (T T T)

    # Ensure C-contiguous boolean array and attach
    sd = np.ascontiguousarray(sd, dtype=bool)
    atoms.set_array("selective_dynamics", sd)

    # Also set an ASE FixAtoms constraint matching the frozen core
    frozen = np.setdiff1d(np.arange(n), relax, assume_unique=False)
    atoms.set_constraint(FixAtoms(indices=frozen))

    # Verify
    sd2 = atoms.get_array("selective_dynamics")
    if sd2.shape != (n, 3) or sd2.dtype != bool:
        raise RuntimeError(f"Failed to set selective_dynamics properly; got shape={sd2.shape}, dtype={sd2.dtype}")

    return atoms

def write_pos_cart(atoms,outpos):
    """
    write poscar to a specific file in cartesian coordinates
    """
    write(outpos, atoms, format='vasp', sort=True, vasp6=True, direct=False)

def write_pos_direct(atoms,outpos):
    """
    write poscar to a specific file in direct coordinates
    """
    write(outpos, atoms, format='vasp', sort=True, vasp6=True, direct=True)

def build_sc(atoms,nx,ny,nz):
    """
    write poscar to a specific file in direct coordinates
    """
    supercell = atoms.repeat((nx,ny,nz))
    return supercell
# Backwards-compatible alias with corrected spelling
