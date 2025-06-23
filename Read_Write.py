import pickle
import numpy as np
import os
import shutil
from ast import literal_eval

###########
## Other ##
###########
def center_geom(coord: list, origin_atom_idx: int):
    ref = np.array(coord[origin_atom_idx])
    new_coord = np.array(coord)-ref
    return new_coord

def save_list_as_text(inplist: list, pathfile: str=os.getcwd()+"outfile.txt"):
    with open(pathfile, "w") as fil:
        for l in inplist:
            print(l, file=fil)

def read_user_input(message: str, rtext: bool=False, rtext_options: list=[], rtype: bool=False, rtype_options: list=[], limit_attempts: bool=True, attempts: int=3, debug: int=0):
    att = 0
    correct = False
    if limit_attempts:
        while att < attempts and not correct:
            opt = input(message)
            if debug > 0: print(f"Read opt={opt}, of type={type(opt)}, rtype={rtype}, rtext={rtext}")
            if rtext and not rtype:
                if opt in rtext_options:                                correct = True
                else:                                                   correct = False; att += 1
            elif rtype and not rtext:
                try:        opt = literal_eval(opt); isstr = False
                except:     isstr = True
                if debug > 0: print(f"isstr={isstr}, opt={opt}, type={type(opt)}")

                if type(opt) in rtype_options:                          correct = True
                else:                                                   correct = False; att += 1

            elif rtype and rtext:
                if type(opt) in rtype_options and opt in rtext_options: correct = True
                else:                                                   correct = False; att += 1
            else: correct = True
            if debug > 0: print(f"correct={correct}, #attempt={att}")

            if not correct: print(f"Please, try again. Options are: {rtext_options}")
        if correct: return opt
        else:       return None

##############
## Binaries ##
##############
def load_binary(pathfile):
    import pickle
    with open(pathfile, "rb") as pickle_file:
        binary = pickle.load(pickle_file)
    return binary


def save_binary(variable, pathfile, backup: bool=False):
    pathfile = pathfile.replace("lustre","home")
    if not backup:
        try: 
            file = open(pathfile,'wb')
            pickle.dump(variable,file)
            file.close()
        except Exception as exc:
            print("Error Saving Binary for pathfile:", pathfile)
            print(exc)
    if backup:
        print("Backup not implemented, file not saved")

#####################
## Plain XYZ Files ##
#####################
def read_xyz(xyz_file):
    assert(xyz_file[-4:] == ".xyz")
    labels = []
    coord = []
    try:    xyz = open(xyz_file, "r")
    except: print("Could not read xyz file: {0}".format(xyz_file))
    n_atoms = xyz.readline()
    title = xyz.readline()
    for line in xyz:
        line_data = line.split()
        if len(line_data) == 4:
            label, x, y, z = line.split()
            coord.append([float(x), float(y), float(z)])
            labels.append(label)
        else:
            print("I can't read the xyz. It has =/ than 4 columns")
    xyz.close()
    return labels, coord

def writexyz(fdir, fname, labels, coord, charge: int=0, spin: int=1):
    if fdir[-1] != "/":
        fdir = fdir + "/"
    natoms = len(labels)
    fullname = fdir + fname
    with open(fullname, "w") as fil:
        print(natoms, file=fil)
        print(charge, spin, file=fil)
        for idx, l in enumerate(labels):
            print("%s  %.6f  %.6f  %.6f" % (l, coord[idx][0], coord[idx][1], coord[idx][2]),file=fil)

#########################
## Custom ExtXYZ Files ##
#########################
def write_xyz_forces_energy(fdir, fname, labels, coord, forces, energy, charge: int=0, spin: int=1, other=None):
    if fdir[-1] != "/":
        fdir = fdir + "/"
    natoms = len(labels)
    fullname = fdir + fname
    if os.path.isfile(fullname): mode = 'a'
    else:                        mode = 'w'
    with open(fullname, mode) as fil:
        print(natoms, file=fil)
        if other is None: print(charge, spin, energy, file=fil)
        else:             print(charge, spin, energy, other, file=fil)
        for idx, l in enumerate(labels):
            print("%s  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f" % (l, coord[idx][0], coord[idx][1], coord[idx][2], forces[idx][0], forces[idx][1], forces[idx][2]),file=fil)

#####
def write_data_MACE_extxyz(fdir, fname, labels, coord, forces, energy, charge: int=0, spin: int=1, other=None):
    other = other if isinstance(other,str) else ''
    if fdir[-1] != "/": fdir = fdir + "/"
    natoms = len(labels)
    fullname = fdir + fname
    if os.path.isfile(fullname): mode = 'a'
    else:                        mode = 'w'
    frmt_atoms = 3*'{:>+16.10f}'
    total_fmt  = '{:6}' + frmt_atoms + '{:8}' + frmt_atoms 
    with open(fullname, mode) as fil:
        print(natoms, file=fil)
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

#####################################
## I don't know where this is used ##
#####################################
def prepare_specie_figure(specie, bond_thr):
    import plotly.graph_objects as go
    from scipy.spatial.distance import cdist
    from Scope.Elementdata import ElementData  
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
