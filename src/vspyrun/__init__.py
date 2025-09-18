"""vspyrun — helpers for VASP IO and setup
"""
__version__ = "0.1.0"

# Re-export public symbols from submodules
from .file import *
from .incar import *
from .kpoints import *
from .poscar import *
from .potcat import *

# Curated public API (star-import safe)
__all__ = ['__version__', '_strip_blank_edges', 'add_ads', 'add_vacuum', 'append_unique_line', 'apply_sd', 'bool_flag', 'build_inver_sym', 'build_mir_sym', 'center_atoms', 'generate_potcar', 'get_file_name', 'incar_ecorr', 'incar_elec', 'incar_geo', 'incar_para', 'incar_soc', 'incar_solv', 'incar_sys_bs', 'incar_sys_rel', 'incar_sys_scf', 'join_incar_blocks', 'keep_at_el', 'kpts_bs_split', 'kpts_scf', 'rem_at_z', 'rep_lat_para', 'shift_min_z', 'write_file', 'write_incar', 'write_pos_cart', 'write_pos_direct', 'read_pos']
