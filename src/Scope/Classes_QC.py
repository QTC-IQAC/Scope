import numpy as np
from Scope import Constants
from Scope.Parse_General import search_string
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

######################
##### QC Objects #####
######################
class VNM(object):
    def __init__(self, index: int, freq: float, red_mass: float=1.0, force_cnt: float=0.0, IR_int: float=0.0, sym: str='A'):
        self.index        = index 
        self.freq_cm      = freq                     ## In cm-1
        self.freq         = freq*Constants.cm2har    ## In atomic units 
        self.red_mass     = red_mass                 ## In AMU     as in Gaussian
        self.force_cnt    = force_cnt                ## In mDyne/A as in Gaussian
        self.IR_int       = IR_int                   ## In KM/Mole as in Gaussian
        self.sym          = sym
        self.has_mode     = False

    def set_mode(self, atomidxs: list, atnums: list, xs: list, ys: list, zs: list):
        self.atomidxs     = atomidxs
        self.atnums       = atnums
        self.labels       = [elemdatabase.elementsym[atnum] for atnum in atnums]
        self.masses       = [elemdatabase.elementweight[l] for l in self.labels]
        self.mode         = np.column_stack([xs, ys, zs])
        self.mode_format2 = np.column_stack([xs, ys, zs]).reshape(-1)
        self.has_mode     = True

    def mass_weight_mode(self):
        if not self.has_mode: return None
        mw = np.repeat(np.sqrt(self.masses), 3)
        self.mode_mw      = self.mode * mw[np.newaxis, :]         ## Mass_Weighted Version of the Mode
        return self.mode_mw
    
    def write_dyn(self, initial_coord: list, amplitude: int=10, outfolder: str='./', labels: None=list, name: str=None):
        ## Writes a file with a trajectory representing the displacement of the VNM
        from .Read_Write import write_xyz
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
        if not isinstance(other, type(self)):       return None
        if not self.has_mode or not other.has_mode: return None
        from .Operations.Vecs_and_Mats import normalize
        vnm_a = normalize(self.mode_format2)
        vnm_b = normalize(other.mode_format2)
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
        to_print  = f'-----------------------------\n'
        to_print +=  '   Vibrational Normal Mode   \n'
        to_print += f'-----------------------------\n'
        to_print += f' Index                  = {self.index}\n'
        to_print += f' Freq (cm-1)            = {self.freq_cm}\n'
        to_print += f' IR Intensity (KM/Mole) = {self.IR_int}\n'
        to_print += f' Reduced Mass (AMU)     = {self.red_mass}\n'
        if hasattr(self,"has_mode"): to_print += f' Has Mode               = {self.has_mode}\n'
        else:                        to_print += f' Has Mode               = False\n'
        return to_print

######
class orbital_set(object):
    def __init__(self, name: str) -> None:
        self.name = name
        self.orbitals = []
    
    def get_eigenvalues(self, eigen_search_string: str, lines: list):
        self.eigen_search_string = eigen_search_string
        self.eigen_search_line_numbers, self.eigen_search_line_found = search_string(eigen_search_string, lines, typ='all')
        self.eigen_search_lines = []
        for l in self.eigen_search_line_numbers:
            self.eigen_search_lines.append(lines[l])        
        self.eigenvalues = []
        
    def add_orbital(self, orb):
        self.orbitals.append(orb)

class orbital(object):
    def __init__(self, index: int, occupation: int, energy: float) -> None:
        self.name = index
        self.occupation = occupation
        self.energy = energy

####################################################
##### Import older versions of VNM Class ###########
####################################################
def import_vnm(old_vnm):
    new_vnm = VNM(old_vnm.index, old_vnm.freq_cm, old_vnm.red_mass, old_vnm.force_cnt, old_vnm.IR_int, old_vnm.sym)
    if hasattr(old_vnm,"haseigenvec"):
        if old_vnm.haseigenvec:
            new_vnm.set_mode(old_vnm.atomidxs, old_vnm.atnums, old_vnm.xs, old_vnm.ys, old_vnm.zs)
    return new_vnm

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

def plot_ir_spectrum(vnms, xmin=None, xmax=None, broadening=10.0, points=2000, kind="gaussian"):
    import matplotlib.pyplot as plt
    """
    Simulate and plot IR spectrum from a list of VNM objects.

    Parameters
    ----------
    vnms : list of VNM
        List of vibrational normal mode objects.
    xmin, xmax : float, optional
        Frequency range in cm^-1. If None, automatically set from vnms.
    broadening : float
        Half-width-at-half-maximum (HWHM) for Lorentzian,
        or standard deviation (sigma) for Gaussian (in cm^-1).
    points : int
        Number of points in spectrum.
    kind : str
        "gaussian" or "lorentzian"
    """

    # Extract frequencies and intensities
    freqs       = np.array([v.freq_cm for v in vnms])
    intensities = np.array([v.IR_int for  v in vnms])

    # Define range
    if xmin is None:  xmin = max(0, freqs.min() - 100)
    if xmax is None:  xmax = freqs.max() + 100

    x = np.linspace(xmin, xmax, points)
    spectrum = np.zeros_like(x)

    # Build spectrum
    for f, I in zip(freqs, intensities):
        if kind.lower() == "gaussian":
            spectrum += I * np.exp(-0.5 * ((x - f)/broadening)**2)
        elif kind.lower() == "lorentzian":
            spectrum += I * (broadening**2 / ((x - f)**2 + broadening**2))
        else:
            raise ValueError("kind must be 'gaussian' or 'lorentzian'")

    # Normalize Gaussian so height matches IR_int approximately
    if kind.lower() == "gaussian":
        spectrum *= 1/(broadening*np.sqrt(2*np.pi))
        spectrum *= broadening

    # Plot
    plt.figure(figsize=(8,4))
    plt.plot(x, spectrum, 'k-', lw=1.5)
    plt.fill_between(x, 0, spectrum, color="grey", alpha=0.4)

    # Optional: add sticks
    for f, I in zip(freqs, intensities):
        plt.vlines(f, 0, I, color="r", lw=1, linestyle="--")

    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Intensity (a.u.)")
    plt.title("Simulated IR Spectrum")
    plt.gca().invert_xaxis()  # convention in spectroscopy
    plt.tight_layout()
    plt.show()

    return x, spectrum
