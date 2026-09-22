## This file gathers the classes of quantum chemistry concepts like vibrational normal modes (VNM)

import numpy as np
from scope import constants
from scope.parse_general import search_string
from scope.elementdata import ElementData
elemdatabase = ElementData()

######################
##### QC Objects #####
######################
class VNM(object):
    """
    Represent one vibrational normal mode and its optional displacement vectors.

    The authoritative mode array has shape `(natoms, 3)`. Its coordinate
    convention is recorded by `is_mass_weighted`; conversion methods can either
    update the stored representation or return an independent converted array.

    Attributes:
        object_type (str):              Object category (`"vnm"`).
        index (int):                    Mode index.
        freq_cm (float):                Frequency in cm^-1.
        freq (float):                   Frequency in atomic units.
        red_mass (float):               Reduced mass.
        has_mode (bool):                Whether eigenvectors are available.
        is_mass_weighted (bool):        Whether the stored mode is mass weighted.

    Methods:
        set_mode():                     Store eigenvector data for the mode.
        mass_weight_mode():             Convert the stored mode to mass-weighted coordinates.
        unweight_mode():                Convert the stored mode to Cartesian coordinates.
        get_atomic_participation():     Calculate each atom's contribution to the mode.
        write_dyn():                    Export a displacement trajectory.
        overlap():                      Compute overlap with another mode.
    """
    def __init__(self, index: int, freq: float, red_mass: float=1.0, force_cnt: float=0.0, IR_int: float=0.0, sym: str='A'):
        self.object_type      = "vnm"
        self.index            = index
        self.freq_cm          = freq                     ## In cm-1
        self.freq             = freq*constants.cm2har    ## In atomic units
        self.red_mass         = red_mass                 ## In AMU     as in Gaussian
        self.force_cnt        = force_cnt                ## In mDyne/A as in Gaussian
        self.IR_int           = IR_int                   ## In KM/Mole as in Gaussian
        self.sym              = sym
        self.has_mode         = False
        self.is_mass_weighted = False

    def set_mode(self, atomidxs: list, atnums: list, xs: list, ys: list, zs: list, is_mass_weighted: bool):
        """Store displacement vectors and their mass-weighting convention.

        Parameters:
            atomidxs (list):            Atom identifiers in mode order.
            atnums (list):              Atomic numbers in mode order.
            xs, ys, zs (list):          Cartesian components for every atom.
            is_mass_weighted (bool):    Whether the supplied vectors are mass weighted.

        Returns:
            np.ndarray: Stored mode with shape `(natoms, 3)`.
        """
        if not isinstance(is_mass_weighted, (bool, np.bool_)): raise TypeError("VNM.SET_MODE: is_mass_weighted must be a boolean")
        if not len(atomidxs) == len(atnums) == len(xs) == len(ys) == len(zs): raise ValueError("VNM.SET_MODE: Atom and mode-component arrays must have the same length")
        self.atomidxs         = list(atomidxs)
        self.atnums           = list(atnums)
        self.labels           = [elemdatabase.elementsym[atnum] for atnum in atnums]
        self.masses           = [elemdatabase.elementweight[l] for l in self.labels]
        self.mode             = np.asarray(np.column_stack([xs, ys, zs]), dtype=float)
        self.has_mode         = True
        self.is_mass_weighted = bool(is_mass_weighted)
        return self.mode

    @property
    def mode_format2(self) -> np.ndarray:
        """Return a flattened view of the authoritative mode array."""
        return self.mode.reshape(-1)

    def mass_weight_mode(self, permanent: bool=True) -> np.ndarray:
        """Return the mass-weighted mode, optionally storing the conversion.

        Parameters:
            permanent (bool):           Whether to replace `mode` and update its convention.

        Returns:
            np.ndarray: Mass-weighted mode. A copy is returned when `permanent=False`.
        """
        if not self.has_mode: raise ValueError("VNM.MASS_WEIGHT_MODE: No mode is stored")
        if not isinstance(permanent, (bool, np.bool_)): raise TypeError("VNM.MASS_WEIGHT_MODE: permanent must be a boolean")
        if self.is_mass_weighted: return self.mode if permanent else self.mode.copy()
        weighted_mode = self.mode * np.sqrt(np.asarray(self.masses, dtype=float))[:, np.newaxis]
        if not permanent: return weighted_mode
        self.mode             = weighted_mode
        self.is_mass_weighted = True
        return self.mode

    def unweight_mode(self, permanent: bool=True) -> np.ndarray:
        """Return the Cartesian mode, optionally storing the conversion.

        Parameters:
            permanent (bool):           Whether to replace `mode` and update its convention.

        Returns:
            np.ndarray: Cartesian mode. A copy is returned when `permanent=False`.
        """
        if not self.has_mode: raise ValueError("VNM.UNWEIGHT_MODE: No mode is stored")
        if not isinstance(permanent, (bool, np.bool_)): raise TypeError("VNM.UNWEIGHT_MODE: permanent must be a boolean")
        if not self.is_mass_weighted: return self.mode if permanent else self.mode.copy()
        unweighted_mode = self.mode / np.sqrt(np.asarray(self.masses, dtype=float))[:, np.newaxis]
        if not permanent: return unweighted_mode
        self.mode             = unweighted_mode
        self.is_mass_weighted = False
        return self.mode

    def get_atomic_participation(self) -> np.ndarray:
        """Return each atom's normalized contribution to the mode's kinetic-energy norm.

        Returns:
            np.ndarray: Atomic contributions in mode order, summing to one.
        """
        if not self.has_mode: raise ValueError("VNM.GET_ATOMIC_PARTICIPATION: No mode is stored")
        weighted_mode = self.mass_weight_mode(permanent=False)
        contributions = np.sum(weighted_mode**2, axis=1)
        norm = np.sum(contributions)
        if not np.isfinite(norm) or norm <= 0.0: raise ValueError("VNM.GET_ATOMIC_PARTICIPATION: Mode has an invalid norm")
        return contributions / norm
    
    def write_dyn(self, initial_coord: list, amplitude: int=10, outfolder: str='./', labels: None=list, name: str=None):
        ## Writes a file with a trajectory representing the displacement of the VNM
        from scope.read_write import write_xyz
        if name is None: filename: str="dyn_vnm_"+str(self.index)+".xyz"
        else:            filename: str=name
        if outfolder[-1] != '/': outfolder += '/'
        initial_coord = np.array(initial_coord)
        for f in range(-amplitude,amplitude,1):
            labels = []
            coords = []
            for idx in range(len(initial_coord)):
                vector = self.mode[idx]
                coord  = initial_coord[idx] + vector*f*0.1
                if labels is None: label = elemdatabase.elementsym[vnm.atnums[idx]]
                else:              label = self.labels[idx]
                labels.append(label)
                coords.append(coord)
            write_xyz(outfolder+filename, labels, coords, append=True) 

    def overlap(self, other: object) -> float:
        """Return the absolute normalized overlap in mass-weighted coordinates."""
        if not isinstance(other, type(self)):       return None
        if not self.has_mode or not other.has_mode: return None
        from scope.operations.vecs_and_mats import normalize
        vnm_a = normalize(self.mass_weight_mode(permanent=False).reshape(-1))
        vnm_b = normalize(other.mass_weight_mode(permanent=False).reshape(-1))
        ov    = float(np.abs(np.dot(vnm_a, vnm_b)))
        return ov

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):       return False
        if not self.has_mode or not other.has_mode: return None
        ov = self.overlap(other)
        if ov > 0.99: return True
        else:         return False

    def __ne__(self, other: object) -> bool:
        if not isinstance(other, type(self)):       return False
        if not self.has_mode or not other.has_mode: return None
        ov = self.overlap(other)
        if ov < 0.01: return True
        else:         return False

    def __repr__(self) -> None:
        to_print  = f'----------------------------------\n'
        to_print += f'------   SCOPE VNM Object    -----\n'
        to_print += f'----------------------------------\n'
        to_print += f' Index                  = {self.index}\n'
        to_print += f' Freq (cm-1)            = {self.freq_cm}\n'
        to_print += f' IR Intensity (KM/Mole) = {self.IR_int}\n'
        to_print += f' Reduced Mass (AMU)     = {self.red_mass}\n'
        if hasattr(self,"has_mode"): to_print += f' Has Mode               = {self.has_mode}\n'
        else:                        to_print += f' Has Mode               = False\n'
        if self.has_mode:            to_print += f' Is Mass Weighted       = {getattr(self, "is_mass_weighted", False)}\n'
        return to_print

####################################################
##### Import older versions of VNM Class ###########
####################################################
def import_vnm(old_vnm):
    new_vnm = VNM(old_vnm.index, old_vnm.freq_cm, old_vnm.red_mass, old_vnm.force_cnt, old_vnm.IR_int, old_vnm.sym)
    if hasattr(old_vnm,"haseigenvec"):
        if old_vnm.haseigenvec:
            new_vnm.set_mode(old_vnm.atomidxs, old_vnm.atnums, old_vnm.xs, old_vnm.ys, old_vnm.zs, is_mass_weighted=getattr(old_vnm, "is_mass_weighted", False))
    return new_vnm

##############
## Plotting ##
##############
def plot_overlap_vnms_diagonal(vnmsA: object, vnmsB: object):
    import matplotlib.pyplot as plt
    min_len = min(len(vnmsA), len(vnmsB))
    diagonal_overlaps = [vnmsA[i].overlap(vnmsB[i]) for i in range(min_len)]
    plt.figure(figsize=(8, 6))
    plt.scatter(range(min_len), diagonal_overlaps, c='blue', marker='o')
    plt.xlabel('Mode Index')
    plt.ylabel('Diagonal Overlap')
    plt.title('Diagonal Overlap: vnmsA vs vnmsB')
    plt.grid(True)
    plt.show()

def plot_overlap_vnms(vnmsA: object, vnmsB: object):
    import matplotlib.pyplot as plt
    overlap_matrix = np.array([[a.overlap(b) for b in vnmsB] for a in vnmsA])
    plt.figure(figsize=(10, 8))
    plt.imshow(overlap_matrix, aspect='auto', interpolation='nearest', cmap='viridis')
    plt.colorbar(label='Overlap')
    plt.xlabel('vnmsB index')
    plt.ylabel('vnmsA index')
    plt.title('Overlap Matrix: vnmsA vs vnmsB')
    plt.show()

##############################
## Electronic Excited State ##
##############################
class ExcitedState(object):
    """
    Represent an electronic excited state.

    Attributes:
        object_type (str):              Object category (`"excited_state"`).
        index (int):                    State index.
        energy (float):                 Excitation energy.
        wavelength (float):             Transition wavelength.
        fosc (float):                   Oscillator strength.
        s2 (float):                     Spin expectation value.

    Methods:
        shift_wavelength():             Shift wavelength and update energy.
        shift_energy():                 Shift energy and update wavelength.
        restore():                      Restore original shifted values.
    """
    def __init__(self, index: int, energy: float, wavelength: float, fosc: float, s2: float, debug: int=0) -> None:
        self.object_type       = "excited_state"
        self.index             = index
        self.energy            = energy
        self.wavelength        = wavelength
        self.fosc              = fosc
        self.s2                = s2

    def shift_wavelength(self, shift: float, debug: int=0):
        self.original_wl      = self.wavelength
        self.wavelength       = self.wavelength + shift
        self.original_energy  = self.energy
        self.energy           = constants.hc/self.wavelength
        if debug > 0: print(f"EXC_STATE.SHIFT_WAVELENGTH: wavelength shifted to {self.wavelength} from {self.original_wl}, and energy adapted to {self.energy}")
        return self.wavelength

    def shift_energy(self, shift: float, debug: int=0):
        self.original_energy  = self.energy
        self.energy           = self.energy + shift
        self.original_wl      = self.wavelength
        self.wavelength       = constants.hc/self.energy
        if debug > 0: print(f"EXC_STATE.SHIFT_ENERGY: energy shifted to {self.energy} from {self.original_energy}, and wavelength adapted to {self.wavelength}")
        return self.energy

    def restore(self):
        if hasattr(self,"original_wl"):     self.wavelength  = self.original_wl
        if hasattr(self,"original_energy"): self.energy      = self.original_energy

    def __repr__(self):
        to_print  = f'-------------------------------------------\n'
        to_print +=  '------ SCOPE Elec. Exc. State Object ------\n' 
        to_print  = f'-------------------------------------------\n'
        to_print += f' Index                  = {self.index}\n'
        to_print += f' Energy (eV)            = {self.energy} eV\n'
        to_print += f' Wavelength (nm)        = {self.wavelength} nm\n'
        to_print += f' Oscillator Strength    = {self.fosc}\n'
        return to_print
