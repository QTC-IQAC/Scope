import os
import sys
import pickle
import shutil
import readline
import json
from ast          import literal_eval
from typing       import Any, Callable, List, Optional, Type
from platformdirs import user_config_dir

############
## Hidden ##
############
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        return self
    def __exit__(self, *_):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

###########
## Other ##
###########
def center_geom(coords: list, origin_atom_idx: int):
    ref = np.array(coords[origin_atom_idx])
    new_coords = np.array(coords)-ref
    return new_coords

def save_list_as_text(inplist: list, pathfile: str=os.getcwd()+"outfile.txt"):
    with open(pathfile, "w") as fil:
        for l in inplist:
            print(l, file=fil)

def read_user_input(message: str, rtext_options: Optional[List[str]] = None, rtype_options: Optional[List[Type]] = None, validator: Optional[Callable[[Any], bool]] = None, default: Any = None, case_insensitive: bool = True, limit_attempts: bool = True, attempts: int = 3, debug: int = 0) -> Any:
    """
    Read user input with optional validation rules.

    Parameters:
        message (str):                 Prompt shown to the user.
        rtext_options (list | None):   Allowed string values.
        rtype_options (list | None):   Allowed Python types.
        validator (callable | None):   Extra validation callback.
        default:                       Value returned on empty input.
        case_insensitive (bool):       Whether to compare text case-insensitively.
        limit_attempts (bool):         Whether to stop after a fixed number of attempts.
        attempts (int):                Maximum number of attempts.
        debug (int):                   Verbosity level.

    Returns:
        Any: Validated user input.
    """
    rtext_options = rtext_options or []
    rtype_options = rtype_options or []
    att = 0

    while True:
        raw_input_value = input(message).strip()
        if debug > 0: print(f"READ_USER_INPUT: Raw input: {raw_input_value!r}")

        # Return default if user presses Enter
        if raw_input_value == "" and default is not None:
            if debug > 0: print(f"READ_USER_INPUT: Returning default value: {default!r}")
            return default

        opt = raw_input_value
        valid = True

        # --- Type evaluation ---
        if rtype_options:
            try:
                val = literal_eval(opt)
                opt_type = type(val)
                opt = val
            except Exception:
                opt_type = str  # treat as string if eval fails
            if opt_type not in rtype_options:
                valid = False
            if debug > 0: print(f"READ_USER_INPUT: Type check: {opt_type}, valid={valid}")

        # --- Text options check ---
        if rtext_options:
            comp = str(opt)
            if case_insensitive:
                comp = comp.lower()
                allowed = [x.lower() for x in rtext_options]
            else:
                allowed = rtext_options
            if comp not in allowed:
                valid = False
            if debug > 0: print(f"READ_USER_INPUT: Text check: {comp}, allowed={allowed}, valid={valid}")

        # --- Custom validator ---
        if validator is not None:
            try:
                if not validator(opt):
                    valid = False
            except Exception:
                valid = False
            if debug > 0: print(f"READ_USER_INPUT: Validator check, valid={valid}")

        # --- Success ---
        if valid:
            return opt

        # --- Failed attempt ---
        att += 1
        print(f"Invalid input. Please try again. Options: {rtext_options}")
        if limit_attempts and att >= attempts:
            raise UserInputError(f"Maximum attempts exceeded.", attempts_used=att, max_attempts=attempts)

#################
## JSON Config ##
#################
def get_config_path(debug: int=0):
    """
    Return the path to the SCOPE user config file.

    Parameters:
        debug (int):                   Verbosity level.

    Returns:
        str: Path to the config file.
    """
    config_dir = user_config_dir("scope")
    os.makedirs(config_dir, exist_ok=True)
    config_path = os.path.join(config_dir, "scope_config.json")
    if not os.path.isfile(config_path): 
        if debug > 0: print(f"GET_CONFIG_PATH: creating config file at {config_path}")
        save_json(dict(), config_path) 
    return config_path

def save_to_config(data : dict=dict()):
    from scope.read_write import save_json

    # Load existing data if config exists
    config_path = get_config_path()
    config_dict = load_config()

    # Updates data
    config_dict.update(data)

    # Saves
    save_json(config_dict, config_path)
    return config_dict

def load_config(debug: int=0):
    config_path = get_config_path()
    if debug > 0: print(f"LOAD_CONFIG: trying to load config file from {config_path=}")
    return load_json(config_path)

##########
## JSON ##
##########
def save_json(dict, pathfile):
    dir_name = os.path.dirname(pathfile)
    os.makedirs(dir_name, exist_ok=True)
    with open(pathfile, "w") as f:
        json.dump(dict, f, indent=4) 

def load_json(pathfile, debug: int=0):
    with open(pathfile, "r") as f:
        dict = json.load(f) 
    return dict

def load_environment(name: str):
    config_path         = get_config_path()
    env_path            = load_json(config_path)[f"env_{name}_filepath"]
    env                 = load_binary(env_path)
    return env

##########
## Text ##
##########
def save_text(variable, pathfile):
    dir_name = os.path.dirname(pathfile)
    os.makedirs(dir_name, exist_ok=True)
    with open(pathfile, "w") as f:
        f.write(variable)

##############
## Binaries ##
##############
def load_binary(pathfile):
    import pickle
    with open(pathfile, "rb") as pickle_file:
        binary = pickle.load(pickle_file)
    return binary

def save_binary(variable, pathfile):
    import tempfile
    # Write to a temporary file first
    dir_name = os.path.dirname(pathfile)
    os.makedirs(dir_name, exist_ok=True)
    try:
        with tempfile.NamedTemporaryFile(dir=dir_name, delete=False) as tmp_file:
            pickle.dump(variable, tmp_file)
            temp_name = tmp_file.name
        # If pickle.dump succeeds, replace the original file
        shutil.move(temp_name, pathfile)
    except Exception as exc:
        print("Error Saving Binary for pathfile:", pathfile)
        print(exc)
        # Clean up temp file if exists
        if 'temp_name' in locals() and os.path.exists(temp_name):
            os.remove(temp_name)

#####################
## Plain XYZ Files ##
#####################
def read_xyz(xyz_file):
    assert(xyz_file[-4:] == ".xyz")
    labels = []
    coords = []
    try:    xyz = open(xyz_file, "r").readlines()
    except: print("Could not read xyz file: {0}".format(xyz_file))
    for line in xyz[2:]:
        line_data = line.split()
        if len(line_data) == 4:
            label, x, y, z = line.split()
            coords.append([float(x), float(y), float(z)])
            labels.append(label)
        else: print("I can't read the xyz. It has =/ than 4 columns")
    return labels, coords

def write_xyz(path, labels, coords, charge: int=0, spin: int=1, append: bool=False):
    assert len(labels) == len(coords)
    if append:     mode = "a"
    else:          mode = "w"
    with open(path, mode) as fil:
        print(len(labels), file=fil)
        print(charge, spin, file=fil)
        for idx, l in enumerate(labels):
            print("%s  %.6f  %.6f  %.6f" % (l, coords[idx][0], coords[idx][1], coords[idx][2]),file=fil)

#########################
## Custom ExtXYZ Files ##
#########################
def write_xyz_with_forces_and_energy(path: str, labels: list, coords: list, forces: list, energy: str, charge: int=0, spin: int=1, other=None):
    other = other if isinstance(other,str) else ''
    assert len(labels) == len(coords)
    assert len(labels) == len(forces)
    if os.path.isfile(path): mode = 'a'
    else:                    mode = 'w'
    frmt_atoms = 3*'{:>+16.10f}'
    total_fmt  = '{:6}' + frmt_atoms + '{:8}' + frmt_atoms 
    with open(path, mode) as fil:
        print(len(labels), file=fil)
        print(charge, spin, energy, other, file=fil)
        for idx, l in enumerate(labels):
            print(total_fmt.format(l, *coord[idx], 0, *forces[idx]), file=fil)

#####
def write_data_MACE_extxyz(path: str, labels: list, coords: list, forces: list, energy: str, charge: int=0, spin: int=1, other=None):
    other = other if isinstance(other,str) else ''
    assert len(labels) == len(coords)
    assert len(labels) == len(forces)
    if os.path.isfile(path): mode = 'a'
    else:                    mode = 'w'
    frmt_atoms = 3*'{:>+16.10f}'
    total_fmt  = '{:6}' + frmt_atoms + '{:8}' + frmt_atoms 
    with open(path, mode) as fil:
        print(len(labels), file=fil)
        print(f'Properties=species:S:1:pos:R:3:molID:I:1:forces:R:3 Nmols=1 Comp={fname.split("_")[0]}_{other} charge={charge} energy={energy} pbc="F F F"', file=fil)
        for idx, l in enumerate(labels):
            print(total_fmt.format(l, *coord[idx], 0, *forces[idx]), file=fil)

###########
## Print ##
###########
def print_xyz(labels, coord):
    for idx, l in enumerate(labels):
        print("%s  %.6f  %.6f  %.6f" % (l, coord[idx][0], coord[idx][1], coord[idx][2]))

#########
## ASE ##
#########
def write_ase_cell(symbols, positions, cell=None, filename='structure.cif'):
    from ase import Atoms
    from ase.io import write
    """
    Save a crystal structure using ASE.

    Parameters:
    - symbols: list of str, chemical symbols (e.g., ['C', 'H', 'H', 'H'])
    - positions: list of lists or Nx3 array of atomic positions in Å
    - cell: 3x3 list or array representing the unit cell vectors
    - filename: output file name (e.g., 'structure.cif', 'structure.xyz', 'structure.traj')
    - pbc: bool or 3-tuple, periodic boundary conditions (default: True)
    """
    if cell is not None: atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    else:                atoms = Atoms(symbols=symbols, positions=positions, pbc=False)
    write(filename, atoms)

def read_ase_atoms(filename):
    atoms = read(filename)
    print(f" Structure read from '{filename}':")
    print(f" Number of atoms: {len(atoms)}")
    print(f" Chemical symbols: {atoms.get_chemical_symbols()}")
    print(f" Positions (Å):\n{atoms.get_positions()}")
    print(f" Periodic boundary conditions: {atoms.get_pbc()}")
    return atoms

##########################
## For .view() function ##
##########################
def set_scene(fig, positions, padding=1.0, width: int=500, height: int=500):
    xmin, xmax = positions[:,0].min() - padding, positions[:,0].max() + padding
    ymin, ymax = positions[:,1].min() - padding, positions[:,1].max() + padding
    zmin, zmax = positions[:,2].min() - padding, positions[:,2].max() + padding

    fig.update_layout(scene=dict(
        xaxis  = dict(title='X (Å)', range=[xmin, xmax]),
        yaxis  = dict(title='Y (Å)', range=[ymin, ymax]),
        zaxis  = dict(title='Z (Å)', range=[zmin, zmax]),
    ))

    fig.update_layout(width=width,height=height)

def prepare_specie_figure(specie, bond_thr):
    import numpy as np
    import plotly.graph_objects as go
    from scipy.spatial.distance import cdist
    from scope.elementdata import ElementData  
    elemdatabase = ElementData()

    fig             = go.Figure()

    # Gather Data
    positions       = specie.positions
    symbols         = specie.symbols

    # Calculate bonds
    distances       = cdist(positions, positions)
    bond_indices    = np.where((distances < bond_thr) & (distances > 0))
    unique_bonds    = set()

    for i, j in zip(*bond_indices):
        if i < j:
            unique_bonds.add((i, j))

    # Plot atoms as markers
    fig.add_trace(go.Scatter3d(
        x           = positions[:, 0],
        y           = positions[:, 1],
        z           = positions[:, 2],
        mode        ='markers',
        marker      = dict(
            size        = 10,
            color       = [atom_colors.get(sym, 'gray') for sym in symbols],
            line        = dict(color='black', width=1),
        ),
        hoverinfo   = 'text',
        text        = symbols,
        showlegend  = False
    ))

    # Label atom indices
    fig.add_trace(go.Scatter3d(
        x           = positions[:, 0],
        y           = positions[:, 1],
        z           = positions[:, 2],
        mode        = 'text',
        text        = [str(i) for i in range(len(positions))],
        textfont    = dict(color='black', size=12),
        hoverinfo   = 'none',
        showlegend  = False
    ))

    # Plot bonds as lines and calculate midpoints
    midpoints   = []
    bond_pairs  = []

    for i, j in unique_bonds:
        # Add bond trace
        fig.add_trace(go.Scatter3d(
            x           = [positions[i, 0], positions[j, 0]],
            y           = [positions[i, 1], positions[j, 1]],
            z           = [positions[i, 2], positions[j, 2]],
            mode        = 'lines',
            line        = dict(color='gray', width=5),
            hoverinfo   = 'none',
            showlegend  = False
        ))

        # Calculate midpoints
        midpoint = (positions[i] + positions[j]) / 2
        midpoints.append(midpoint)
        bond_pairs.append((i, j))

    midpoints = np.array(midpoints)

    # Return figure, midpoints, and bond pairs (for bond-related data plotting)
    return fig, midpoints, bond_pairs

###################
## Autocompleter ##
###################
def complete_path(text, state):
    import glob, readline
    buffer = readline.get_line_buffer()

    # If buffer already has content and text is empty, do NOT override it
    if not text and buffer: return None

    text = os.path.expanduser(text)
    matches = glob.glob(text + '*')

    # Add trailing slash for directories
    matches = [m + '/' if os.path.isdir(m) else m for m in matches]
    if state < len(matches):
        return matches[state]
    return None

######
def input_with_default(prompt: str, default: str | None = None) -> str:
    default = "" if default is None else str(default)
    if default:
        prompt = f"{prompt}[{default}] "
    value = input(prompt)
    return value.strip() or default

################
## Exceptions ##
################

class UserInputError(Exception):
    """Raised when user fails to provide valid input after allowed attempts."""
    def __init__(self, message, attempts_used=None, max_attempts=None):
        super().__init__(message)
        self.attempts_used = attempts_used
        self.max_attempts  = max_attempts
