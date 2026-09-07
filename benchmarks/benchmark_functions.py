######################################################################################################################
## Set of Functions employed in the Benchmark notebooks: 1-Import, 2-Filtering_and_Execution and 3-Molecule_Overlap ##
## These are not SCOPE functions, but may use some functionality implemented therein                                ##
######################################################################################################################
import os
import pickle
import random
import time

import networkx as nx
import numpy as np
from cell2mol import xyz2mol
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign, rdMolTransforms

from scope.connectivity import labels2formula
from scope.operations.graphs import build_graph
from scope.overlap import overlap_molecules

#############################################
### Select Cell, XYZ and MOL object paths ###
#############################################
def select_cell_paths(number_entries, dataset_folder):
    ## Function to randomly select cell2mol Cell objects located in the metal-dedicated folders.
    metal_folders        = ['1-Iron', '2-Manganese', '3-Ruthenium', '4-Rhenium', '5-Chromium', '6-Cobalt', '7-Nickel', '8-Copper']
    available_cell_paths = []
    for folder_name in metal_folders:
        folder_path = os.path.join(dataset_folder, folder_name)
        if not os.path.isdir(folder_path):
            continue
        for file_name in os.listdir(folder_path):
            if file_name.endswith('.gmol'):
                cell_path = os.path.join(folder_name, file_name)
                available_cell_paths.append(cell_path)
    if number_entries > len(available_cell_paths):
        raise ValueError(f'Requested {number_entries} entries, but only {len(available_cell_paths)} are available')
    return random.sample(available_cell_paths, number_entries)

def select_xyz_paths(number_xyz, dataset_folder):
    ## Function to randomly select XYZ files located in a dedicated folder.
    available_xyz = []
    for file_name in os.listdir(dataset_folder):
        if file_name.endswith('.xyz'):
            xyz_path = file_name
            available_xyz.append(xyz_path)
    if number_xyz > len(available_xyz):
        raise ValueError(f'Requested {number_xyz} entries, but only {len(available_xyz)} are available')
    return random.sample(available_xyz, number_xyz)

def select_mol_paths(number_mols, summary, summary_path, dataset_folder):
    ## Function to randomly select available GEOM-QM9 pickle files using the dataset summary.
    summary_entries = list(summary.values())
    random.shuffle(summary_entries)
    selected_paths = []
    for entry in summary_entries:
        if not isinstance(entry, dict) or 'pickle_path' not in entry: continue
        stored_path = entry['pickle_path']
        if not isinstance(stored_path, str): continue
        if os.path.isabs(stored_path):
            candidate_paths = [stored_path, os.path.join(dataset_folder, 'qm9', os.path.basename(stored_path)), os.path.join(dataset_folder, os.path.basename(stored_path))]
        else:
            candidate_paths = [os.path.join(os.path.dirname(summary_path), stored_path), os.path.join(dataset_folder, stored_path), os.path.join(os.path.dirname(dataset_folder), stored_path)]
        mol_path = next((path for path in candidate_paths if os.path.isfile(path)), None)
        if mol_path is None: continue
        selected_path = os.path.relpath(mol_path, dataset_folder) if os.path.commonpath([os.path.abspath(mol_path), os.path.abspath(dataset_folder)]) == os.path.abspath(dataset_folder) else mol_path
        if selected_path in selected_paths: continue
        selected_paths.append(selected_path)
        if len(selected_paths) == number_mols: break
    if number_mols > len(selected_paths):
        raise ValueError(f'Requested {number_mols} entries, but only {len(selected_paths)} GEOM-QM9 pickle files are currently available')
    return selected_paths

#############################################
### Select Cell, XYZ and MOL object paths ###
#############################################
def load_geom_rdkit_molecule(mol_path, conformer_index=0):
    ## Load one RDKit molecule from a GEOM molecular-species pickle file.
    with open(mol_path, 'rb') as mol_file:
        mol_data = pickle.load(mol_file)
    if not isinstance(mol_data, dict) or 'conformers' not in mol_data:
        raise ValueError('GEOM pickle does not contain a conformer collection')
    conformers = mol_data['conformers']
    if not conformers:
        raise ValueError('GEOM pickle contains no conformers')
    if conformer_index >= len(conformers):
        raise IndexError(f'Requested conformer {conformer_index}, but the GEOM pickle contains {len(conformers)} conformers')
    if not isinstance(conformers[conformer_index], dict) or 'rd_mol' not in conformers[conformer_index]:
        raise ValueError(f'GEOM conformer {conformer_index} does not contain an RDKit molecule')
    return conformers[conformer_index]['rd_mol']

#########################
### Molecular Overlap ###
#########################
def get_distance_metrics(reference_coord, overlapped_coord):
    ## Calculate RMSD and the largest displacement between corresponding atoms.
    atom_distances         = np.linalg.norm(np.asarray(reference_coord) - np.asarray(overlapped_coord), axis=1)
    rmsd                   = np.sqrt(np.mean(atom_distances**2))
    maximum_distance_index = int(np.argmax(atom_distances))
    maximum_distance       = float(atom_distances[maximum_distance_index])
    return float(rmsd), maximum_distance, maximum_distance_index

def scope_overlap(molecule1, molecule2):
    ## Overlap two SCOPE Molecules while excluding connectivity preparation from the elapsed time.
    if not hasattr(molecule1, 'adjmat'): molecule1.get_adjmatrix()
    if not hasattr(molecule2, 'adjmat'): molecule2.get_adjmatrix()
    bond_orders1 = molecule1.get_bond_order_matrix()
    bond_orders2 = molecule2.get_bond_order_matrix()

    tini = time.perf_counter()
    isgood, _, coord1, _, coord2, _ = overlap_molecules(molecule1.labels, molecule1.coord, molecule2.labels, molecule2.coord, adjmat1=molecule1.adjmat, adjmat2=molecule2.adjmat, bond_orders1=bond_orders1, bond_orders2=bond_orders2)
    if not isgood: raise ValueError('SCOPE could not overlap the molecules')

    rmsd, maximum_distance, maximum_distance_index = get_distance_metrics(coord1, coord2)
    elapsed_time = time.perf_counter() - tini
    return elapsed_time, rmsd, maximum_distance, maximum_distance_index, np.asarray(coord2)

def rdkit_overlap(rdkit_mol1, rdkit_mol2):
    ## Overlap two pre-built RDKit Molecules and calculate mapped-atom distance metrics.
    tini = time.perf_counter()
    rmsd, transform, atom_map = rdMolAlign.GetBestAlignmentTransform(rdkit_mol2, rdkit_mol1)
    rdMolTransforms.TransformConformer(rdkit_mol2.GetConformer(), transform)

    reference_coord         = np.asarray(rdkit_mol1.GetConformer().GetPositions())
    overlapped_coord        = np.asarray(rdkit_mol2.GetConformer().GetPositions())
    mapped_reference_coord  = np.asarray([reference_coord[reference_index] for mobile_index, reference_index in atom_map])
    mapped_overlapped_coord = np.asarray([overlapped_coord[mobile_index] for mobile_index, reference_index in atom_map])
    calculated_rmsd, maximum_distance, maximum_distance_index = get_distance_metrics(mapped_reference_coord, mapped_overlapped_coord)
    if not np.isclose(rmsd, calculated_rmsd, atol=1e-6, rtol=0): raise ValueError('RDKit RMSD does not coincide with the mapped atom distances')

    maximum_distance_atom_pair = tuple(int(index) for index in atom_map[maximum_distance_index])
    elapsed_time = time.perf_counter() - tini
    return elapsed_time, float(rmsd), maximum_distance, maximum_distance_atom_pair, overlapped_coord

########################
### Validate Objects ###
########################
def validate_cell_molecule(cell2mol_molecule, scope_molecule):
    """
    Compare a cell2mol molecule with its imported SCOPE representation.

    Parameters:
        cell2mol_molecule (object):     Molecule stored in the original cell2mol Cell.
        scope_molecule (object):        Molecule imported into SCOPE.

    Returns:
        bool:                           True when the representations are consistent.
        str:                            Empty when consistent; otherwise, a report of all differences.
    """
    differences = []

    # Compare the molecular composition and atom order.
    try:
        cell2mol_labels = list(cell2mol_molecule.labels)
        scope_labels    = list(scope_molecule.labels)
        if labels2formula(cell2mol_labels) != scope_molecule.formula:
            differences.append(f'Different formula: cell2mol={labels2formula(cell2mol_labels)}, SCOPE={scope_molecule.formula}')
        if cell2mol_labels != scope_labels:
            if len(cell2mol_labels) != len(scope_labels):
                differences.append(f'Different number of atom labels: cell2mol={len(cell2mol_labels)}, SCOPE={len(scope_labels)}')
            different_atoms = [(index, cell2mol_label, scope_label) for index, (cell2mol_label, scope_label) in enumerate(zip(cell2mol_labels, scope_labels)) if cell2mol_label != scope_label]
            if different_atoms:
                differences.append(f'Different atom labels or order (index, cell2mol, SCOPE): {different_atoms}')
    except Exception as error:
        differences.append(f'Atom or formula comparison failed ({type(error).__name__}): {error}')

    # Compare Cartesian coordinates and the two adjacency representations.
    try:
        cell2mol_coordinates = np.asarray(cell2mol_molecule.coord)
        scope_coordinates    = np.asarray(scope_molecule.coord)
        if cell2mol_coordinates.shape != scope_coordinates.shape:
            differences.append(f'Different coordinate shapes: cell2mol={cell2mol_coordinates.shape}, SCOPE={scope_coordinates.shape}')
        elif not np.allclose(cell2mol_coordinates, scope_coordinates, rtol=0.0, atol=1e-6):
            atom_displacements = np.linalg.norm(cell2mol_coordinates - scope_coordinates, axis=1)
            largest_difference = int(np.argmax(atom_displacements))
            differences.append(f'Different Cartesian coordinates: maximum atom displacement={atom_displacements[largest_difference]:.6f} Angstrom at index {largest_difference}')
    except Exception as error:
        differences.append(f'Coordinate comparison failed ({type(error).__name__}): {error}')

    try:
        cell2mol_adjacency = np.asarray(cell2mol_molecule.conmat)
        scope_adjacency    = np.asarray(scope_molecule.adjmat)
        if cell2mol_adjacency.shape != scope_adjacency.shape:
            differences.append(f'Different covalent adjacency shapes: cell2mol={cell2mol_adjacency.shape}, SCOPE={scope_adjacency.shape}')
        elif not np.array_equal(cell2mol_adjacency, scope_adjacency):
            missing_scope_bonds    = [f'{index1}:{cell2mol_labels[index1]}-{index2}:{cell2mol_labels[index2]}' for index1 in range(len(cell2mol_labels)) for index2 in range(index1 + 1, len(cell2mol_labels)) if cell2mol_adjacency[index1, index2] > 0 and scope_adjacency[index1, index2] == 0]
            additional_scope_bonds = [f'{index1}:{cell2mol_labels[index1]}-{index2}:{cell2mol_labels[index2]}' for index1 in range(len(cell2mol_labels)) for index2 in range(index1 + 1, len(cell2mol_labels)) if cell2mol_adjacency[index1, index2] == 0 and scope_adjacency[index1, index2] > 0]
            if missing_scope_bonds:    differences.append(f'Covalent bonds missing from SCOPE: {missing_scope_bonds}')
            if additional_scope_bonds: differences.append(f'Additional SCOPE covalent bonds: {additional_scope_bonds}')

        cell2mol_metal_adjacency = np.asarray(cell2mol_molecule.mconmat)
        scope_metal_adjacency    = np.asarray(scope_molecule.madjmat)
        if cell2mol_metal_adjacency.shape != scope_metal_adjacency.shape:
            differences.append(f'Different metal adjacency shapes: cell2mol={cell2mol_metal_adjacency.shape}, SCOPE={scope_metal_adjacency.shape}')
        elif not np.array_equal(cell2mol_metal_adjacency, scope_metal_adjacency):
            missing_scope_bonds    = [f'{index1}:{cell2mol_labels[index1]}-{index2}:{cell2mol_labels[index2]}' for index1 in range(len(cell2mol_labels)) for index2 in range(index1 + 1, len(cell2mol_labels)) if cell2mol_metal_adjacency[index1, index2] > 0 and scope_metal_adjacency[index1, index2] == 0]
            additional_scope_bonds = [f'{index1}:{cell2mol_labels[index1]}-{index2}:{cell2mol_labels[index2]}' for index1 in range(len(cell2mol_labels)) for index2 in range(index1 + 1, len(cell2mol_labels)) if cell2mol_metal_adjacency[index1, index2] == 0 and scope_metal_adjacency[index1, index2] > 0]
            if missing_scope_bonds:    differences.append(f'Metal bonds missing from SCOPE: {missing_scope_bonds}')
            if additional_scope_bonds: differences.append(f'Additional SCOPE metal bonds: {additional_scope_bonds}')
    except Exception as error:
        differences.append(f'Adjacency comparison failed ({type(error).__name__}): {error}')

    # Validate the transition-metal classification and required substructures.
    try:
        cell2mol_iscomplex = cell2mol_molecule.type == 'Complex'
        scope_iscomplex    = scope_molecule.iscomplex
        scope_ligands      = getattr(scope_molecule, 'ligands', [])
        scope_metals       = getattr(scope_molecule, 'metals', [])
        if cell2mol_iscomplex != scope_iscomplex:
            differences.append(f'Different TMC classification: cell2mol={cell2mol_iscomplex}, SCOPE={scope_iscomplex}')
        if scope_iscomplex and not scope_ligands:
            differences.append('SCOPE TMC does not contain ligands')
        if scope_iscomplex and not scope_metals:
            differences.append('SCOPE TMC does not contain metals')
    except Exception as error:
        differences.append(f'TMC classification failed ({type(error).__name__}): {error}')
        cell2mol_iscomplex = False
        scope_ligands      = []
        scope_metals       = []

    # Ligand and metal details only apply to cell2mol transition-metal complexes.
    if cell2mol_iscomplex:
        try:
            if len(cell2mol_molecule.ligandlist) != len(scope_ligands):
                differences.append(f'Different number of ligands: cell2mol={len(cell2mol_molecule.ligandlist)}, SCOPE={len(scope_ligands)}')
            else:
                cell2mol_ligands = sorted((labels2formula(list(ligand.labels)), ligand.totcharge) for ligand in cell2mol_molecule.ligandlist)
                scope_ligands    = sorted((ligand.formula, ligand.charge) for ligand in scope_ligands)
                if cell2mol_ligands != scope_ligands:
                    differences.append(f'Different ligand formulas or charges: cell2mol={cell2mol_ligands}, SCOPE={scope_ligands}')
        except Exception as error:
            differences.append(f'Ligand comparison failed ({type(error).__name__}): {error}')

        try:
            if len(cell2mol_molecule.metalist) != len(scope_metals):
                differences.append(f'Different number of metals: cell2mol={len(cell2mol_molecule.metalist)}, SCOPE={len(scope_metals)}')
            else:
                cell2mol_metals = sorted((metal.label, metal.totcharge) for metal in cell2mol_molecule.metalist)
                scope_metals    = sorted((metal.formula, metal.charge) for metal in scope_metals)
                if cell2mol_metals != scope_metals:
                    differences.append(f'Different metal formulas or charges: cell2mol={cell2mol_metals}, SCOPE={scope_metals}')

                cell2mol_coordination = sorted((metal.label, metal.totcharge, len(metal.coord_sphere)) for metal in cell2mol_molecule.metalist)
                scope_coordination    = sorted((metal.formula, metal.charge, len(metal.get_coord_sphere())) for metal in getattr(scope_molecule, 'metals', []))
                if cell2mol_coordination != scope_coordination:
                    differences.append(f'Different metal coordination numbers: cell2mol={cell2mol_coordination}, SCOPE={scope_coordination}')
        except Exception as error:
            differences.append(f'Metal comparison failed ({type(error).__name__}): {error}')

    isvalid        = len(differences) == 0
    report_message = '' if isvalid else '\n'.join(f'- {difference}' for difference in differences)
    return isvalid, report_message

## Compare the chemical information in an RDKit molecule and its SCOPE representation.
def compare_rdkit_scope_molecules(rdkit_molecule, scope_molecule):
    """
    Compare an RDKit molecule with its SCOPE representation.

    Parameters:
        rdkit_molecule (object):       RDKit molecule used as the reference.
        scope_molecule (object):       SCOPE molecule to validate.

    Returns:
        bool:                          True when the representations are consistent.
        str:                           Empty when consistent; otherwise, a report of all differences.
    """
    differences = []

    # SCOPE represents every atom explicitly, so materialize any implicit RDKit hydrogens before comparing the representations.
    try:
        has_3d             = rdkit_molecule.GetNumConformers() > 0 and rdkit_molecule.GetConformer().Is3D()
        implicit_hydrogens = sum(atom.GetNumImplicitHs() for atom in rdkit_molecule.GetAtoms())
        if implicit_hydrogens > 0: rdkit_molecule = Chem.AddHs(Chem.Mol(rdkit_molecule), addCoords=has_3d)
    except Exception as error:
        report_message = f'- RDKit implicit-hydrogen expansion failed ({type(error).__name__}): {error}'
        return False, report_message

    # Fundamental atom information is required by all subsequent RDKit checks.
    try:
        rdkit_atoms       = list(rdkit_molecule.GetAtoms())
        rdkit_natoms      = rdkit_molecule.GetNumAtoms()
        rdkit_labels      = [atom.GetSymbol() for atom in rdkit_atoms]
        rdkit_atomic_nums = [atom.GetAtomicNum() for atom in rdkit_atoms]
    except Exception as error:
        report_message = f'- RDKit atom inspection failed ({type(error).__name__}): {error}'
        return False, report_message

    try:
        scope_natoms      = scope_molecule.natoms
        scope_labels      = list(scope_molecule.labels)
        scope_atomic_nums = list(scope_molecule.get_atomic_numbers())
    except Exception as error:
        report_message = f'- SCOPE atom inspection failed ({type(error).__name__}): {error}'
        return False, report_message

    # Sanitization verifies valence, aromaticity, and related RDKit chemical constraints without modifying the original object.
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(Chem.Mol(rdkit_molecule))
    except Exception as error:
        differences.append(f'RDKit sanitization failed ({type(error).__name__}): {error}')

    # Atom count, identity, and order must be preserved.
    if rdkit_natoms != scope_natoms:
        differences.append(f'Different number of explicit atoms: RDKit={rdkit_natoms}, SCOPE={scope_natoms}')
    if rdkit_atomic_nums != scope_atomic_nums or rdkit_labels != scope_labels:
        different_atoms = [(index, rdkit_label, scope_label) for index, (rdkit_label, scope_label) in enumerate(zip(rdkit_labels, scope_labels)) if rdkit_label != scope_label]
        if different_atoms:
            differences.append(f'Different atom identities or order (index, RDKit, SCOPE): {different_atoms}')

    # Verify that the comparison reference no longer contains implicit hydrogens.
    try:
        with rdBase.BlockLogs():
            implicit_hydrogens = [(atom.GetIdx(), atom.GetSymbol(), atom.GetNumImplicitHs()) for atom in rdkit_atoms if atom.GetNumImplicitHs() > 0]
        if implicit_hydrogens:
            differences.append(f'RDKit contains implicit hydrogens (index, element, count): {implicit_hydrogens}')
    except Exception as error:
        differences.append(f'RDKit implicit-hydrogen inspection failed ({type(error).__name__}): {error}')

    # Coordinates are compared only when the RDKit object contains a conformer.
    try:
        if rdkit_molecule.GetNumConformers() > 0:
            conformer          = rdkit_molecule.GetConformer()
            rdkit_coordinates  = np.array([list(conformer.GetAtomPosition(index)) for index in range(rdkit_natoms)])
            scope_coordinates  = np.asarray(scope_molecule.coord)
            if rdkit_coordinates.shape != scope_coordinates.shape:
                differences.append(f'Different coordinate shapes: RDKit={rdkit_coordinates.shape}, SCOPE={scope_coordinates.shape}')
            elif not np.allclose(rdkit_coordinates, scope_coordinates, rtol=0.0, atol=1e-6):
                atom_displacements = np.linalg.norm(rdkit_coordinates - scope_coordinates, axis=1)
                largest_difference = int(np.argmax(atom_displacements))
                differences.append(f'Different Cartesian coordinates: maximum atom displacement={atom_displacements[largest_difference]:.6f} Angstrom at index {largest_difference}:{rdkit_labels[largest_difference]}')
    except Exception as error:
        differences.append(f'RDKit coordinate inspection failed ({type(error).__name__}): {error}')

    # Compare binary covalent connectivity and identify the individual mismatching bonds.
    rdkit_adjacency = None
    try:
        rdkit_adjacency = Chem.GetAdjacencyMatrix(rdkit_molecule).astype(int)
    except Exception as error:
        differences.append(f'RDKit adjacency extraction failed ({type(error).__name__}): {error}')
    scope_adjacency = getattr(scope_molecule, 'adjmat', None)
    if scope_adjacency is None:
        differences.append('SCOPE adjacency matrix is unavailable')
    else:
        try:
            scope_adjacency = (np.asarray(scope_adjacency) > 0).astype(int)
            if not np.array_equal(scope_adjacency, scope_adjacency.T): differences.append('SCOPE adjacency matrix is not symmetric')
            if np.any(np.diag(scope_adjacency) != 0):                 differences.append('SCOPE adjacency matrix contains self-adjacencies')
            if rdkit_adjacency is not None:
                if rdkit_adjacency.shape != scope_adjacency.shape:
                    differences.append(f'Different adjacency shapes: RDKit={rdkit_adjacency.shape}, SCOPE={scope_adjacency.shape}')
                else:
                    missing_scope_bonds    = [f'{index1}:{rdkit_labels[index1]}-{index2}:{rdkit_labels[index2]}' for index1 in range(rdkit_natoms) for index2 in range(index1 + 1, rdkit_natoms) if rdkit_adjacency[index1, index2] > 0 and scope_adjacency[index1, index2] == 0]
                    additional_scope_bonds = [f'{index1}:{rdkit_labels[index1]}-{index2}:{rdkit_labels[index2]}' for index1 in range(rdkit_natoms) for index2 in range(index1 + 1, rdkit_natoms) if rdkit_adjacency[index1, index2] == 0 and scope_adjacency[index1, index2] > 0]
                    if missing_scope_bonds:    differences.append(f'Bonds missing from SCOPE: {missing_scope_bonds}')
                    if additional_scope_bonds: differences.append(f'Additional SCOPE bonds: {additional_scope_bonds}')
        except Exception as error:
            differences.append(f'SCOPE adjacency inspection failed ({type(error).__name__}): {error}')

    # Compare formal bond orders only when bond objects have been created in SCOPE.
    if getattr(scope_molecule, 'has_bonds', False):
        try:
            rdkit_bond_orders = np.zeros((rdkit_natoms, rdkit_natoms))
            for bond in rdkit_molecule.GetBonds():
                index1 = bond.GetBeginAtomIdx()
                index2 = bond.GetEndAtomIdx()
                rdkit_bond_orders[index1, index2] = bond.GetBondTypeAsDouble()
                rdkit_bond_orders[index2, index1] = bond.GetBondTypeAsDouble()
            scope_bond_orders = scope_molecule.get_bond_order_matrix()
            if scope_bond_orders is None:
                differences.append('SCOPE bond-order matrix is unavailable')
            elif rdkit_bond_orders.shape != np.shape(scope_bond_orders):
                differences.append(f'Different bond-order shapes: RDKit={rdkit_bond_orders.shape}, SCOPE={np.shape(scope_bond_orders)}')
            else:
                different_bond_orders = [(f'{index1}:{rdkit_labels[index1]}-{index2}:{rdkit_labels[index2]}', rdkit_bond_orders[index1, index2], scope_bond_orders[index1, index2]) for index1 in range(rdkit_natoms) for index2 in range(index1 + 1, rdkit_natoms) if not np.isclose(rdkit_bond_orders[index1, index2], scope_bond_orders[index1, index2])]
                if different_bond_orders:
                    differences.append(f'Different bond orders (bond, RDKit, SCOPE): {different_bond_orders}')
        except Exception as error:
            differences.append(f'Bond-order comparison failed ({type(error).__name__}): {error}')

    # Molecular charge and the number of unpaired electrons must coincide.
    try:
        rdkit_charge   = Chem.GetFormalCharge(rdkit_molecule)
        rdkit_radicals = sum(atom.GetNumRadicalElectrons() for atom in rdkit_atoms)
    except Exception as error:
        differences.append(f'RDKit charge or radical inspection failed ({type(error).__name__}): {error}')
    else:
        try:
            scope_charge = scope_molecule.charge
            scope_spin   = scope_molecule.spin
            if scope_charge is None:
                differences.append('SCOPE molecular charge is unavailable')
            elif rdkit_charge != scope_charge:
                differences.append(f'Different molecular charge: RDKit={rdkit_charge}, SCOPE={scope_charge}')
            if scope_spin is None:
                differences.append('SCOPE molecular spin is unavailable')
            elif rdkit_radicals != scope_spin:
                differences.append(f'Different number of unpaired electrons: RDKit={rdkit_radicals}, SCOPE={scope_spin}')
        except Exception as error:
            differences.append(f'SCOPE charge or spin inspection failed ({type(error).__name__}): {error}')

    isconsistent   = len(differences) == 0
    report_message = '' if isconsistent else '\n'.join(f'- {difference}' for difference in differences)
    return isconsistent, report_message

def validate_xyz_molecule(scope_molecule):
    """
    Validate a SCOPE Molecule created from an XYZ file.

    Parameters:
        scope_molecule (object):        SCOPE molecule to validate.

    Returns:
        bool:                           True when the molecule passes all checks.
        str:                            Empty when valid; otherwise, a report of all errors.
    """
    differences = []
    try:
        scope_adjacency = getattr(scope_molecule, 'adjmat', None)
        expected_shape  = (scope_molecule.natoms, scope_molecule.natoms)
        scope_iscomplex = scope_molecule.iscomplex
    except Exception as error:
        report_message = f'- SCOPE molecule inspection failed ({type(error).__name__}): {error}'
        return False, report_message
    adjacency_ready = False

    # Validate the SCOPE adjacency before using it to generate an RDKit molecule.
    if scope_adjacency is None:
        differences.append('SCOPE adjacency matrix is unavailable')
    else:
        try:
            scope_adjacency = np.asarray(scope_adjacency)
            if scope_adjacency.shape != expected_shape:
                differences.append(f'Incorrect adjacency shape: expected={expected_shape}, found={scope_adjacency.shape}')
            else:
                adjacency_ready = True
                if not np.array_equal(scope_adjacency, scope_adjacency.T): differences.append('Adjacency matrix is not symmetric')
                if np.any(np.diag(scope_adjacency) != 0):                 differences.append('Adjacency matrix contains self-adjacencies')
                molecular_graph = build_graph(scope_adjacency)
                if not nx.is_connected(molecular_graph):
                    differences.append('Molecular graph is not connected')
        except Exception as error:
            differences.append(f'Adjacency validation failed ({type(error).__name__}): {error}')

    if scope_iscomplex:
        differences.append('Molecule is classified as a transition-metal complex')

    # xyz2mol assigns a chemically valid RDKit representation to the SCOPE adjacency when possible.
    if adjacency_ready and not scope_iscomplex:
        try:
            with rdBase.BlockLogs():
                rdkit_molecules = xyz2mol.xyz2mol(scope_molecule.get_atomic_numbers(), scope_molecule.coord, scope_adjacency, 1.3, charge=0)
            if not rdkit_molecules:
                raise ValueError('xyz2mol did not return an RDKit molecule')
            scope_molecule.rdkit_obj = rdkit_molecules[0]
            isconsistent, comparison_report = compare_rdkit_scope_molecules(scope_molecule.rdkit_obj, scope_molecule)
            if not isconsistent:
                comparison_report = comparison_report.replace('\n', '; ')
                differences.append(f'RDKit-SCOPE comparison failed: {comparison_report}')
        except Exception as error:
            differences.append(f'RDKit generation failed ({type(error).__name__}): {error}')

    isvalid        = len(differences) == 0
    report_message = '' if isvalid else '\n'.join(f'- {difference}' for difference in differences)
    return isvalid, report_message
