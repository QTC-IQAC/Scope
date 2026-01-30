import numpy as np
from scope                            import constants
from scope.classes_data               import *
from scope.operations.dicts_and_lists import range2list

###########
def get_Svib(freqs: list, temp: float, freq_units: str='au', outunits: str='au', typ: str='default', FR_cutoff: int=100, FR_alpha: int=4, nmol: int=1, debug: int=0):
    ## Temperature must be provided in K

    ## If typ = 'free-rotor' or 'fr': 
        # it uses the Quasi Rigid Rotor Harmonic Oscillator approach of Grimme...
        # Described in: Grimme, S. Chem. Eur. J., 2012, 18, 9955–9964.
        # Equations of this manuscript are references along the function
        # With Chai and Head-Gordon's damping function, using FR_cutoff and FR_alpha 
        # FR_cutoff must be provided in cm-1 units

    ## else, it uses the Harmonic Oscillator expressions

    ## Function works with freqs in s-1 first, and then au
    if   freq_units.lower() == 'cm':  freqs_mod = np.array(freqs) * constants.cm2s_1
    elif freq_units.lower() == 'ev':  freqs_mod = np.array(freqs) * constants.eV2s_1
    elif freq_units.lower() == 'au':  freqs_mod = np.array(freqs) * constants.har2s_1
    elif freq_units.lower() == 's_1': freqs_mod = np.array(freqs)
    else: raise ValueError(f"GET_Svib: can't understand input units of frequencies {freq_units}")

    ## Prepare the parameter
    FR_cutoff *= constants.cm2s_1
        
    total=0.0
    for idx, f in enumerate(freqs_mod):                                                      # Freqs in au
        #if f < 0.0: f = np.abs(f)                                                        # Converts Negative frequencies to positive
        #if np.abs(f) < 1.0*constants.cm2s_1: continue                                    # Ignores frequencies below 1 cm-1

        ## Free Rotor Term
        if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
            bav=1.0000E-44                                                               # Kg·m2, parameter in manuscript
            mu=constants.planck_Js/(8*(np.pi)**2*f)                                      # [J·s]/[s_1] = J·s2 = Kg·m2, equation 4
            mu_prime=(mu*bav)/(mu+bav)                                                   # Kg·m2, equation 5 
            a=(8*np.pi**3)*mu_prime*constants.boltz_J*temp/(constants.planck_Js**2)      # Dimensionless, equation 6 numerator. Internally: [J·s2]·[J·K_1]·[K] / [J·s]^2
            b=np.sqrt(a)                                                                 # Dimensionless, equation 6, square root
            c=np.log(b)                                                                  # Dimensionless, equation 6, ln()
            Svib_FR=constants.boltz_au*(c+0.5)                                           # Hartree/K/molecule, equation 6
            weight_FR = 1-(1/(1+(FR_cutoff/f)**FR_alpha))                                  # Dimensionless
            if debug > 1: 
                print(f"\tGET_Svib: FR Term: {idx=} {f=}: {mu=}, {mu_prime=}, {a=}, {Svib_FR=}, {weight_FR=}")
            elif debug == 1: 
                print(f"\tGET_Svib: FR Term: f={freqs[idx]:.4f}: {Svib_FR:.4e}, {np.round(weight_FR,4)}")
        else: 
            Svib_FR   = 0.0
            weight_FR = 0.0
    
        ## Below, we work in atomic units. We convert the frequency to hartree
        f /= constants.har2s_1

        ## Harmonic Oscillator Term
        exponential_pos=np.exp(f/(constants.boltz_au*temp))                              # Dimensionless
        exponential_neg=np.exp(-f/(constants.boltz_au*temp))                             # Dimensionless
        fstterm=f/(temp*(exponential_pos-1))                                             # Hartree/molecule*K
        scnterm=-constants.boltz_au*np.log(1-exponential_neg)                            # Hartree/molecule*K
        Svib_HO = fstterm+scnterm
        weight_HO = 1 - weight_FR
        if debug > 1: 
            print(f"\tGET_Svib: HO Term: {idx=} {f=}: {exponential_pos=}, {exponential_neg=}, {fstterm=}, {scnterm=}, {Svib_HO=}, {weight_HO=}")
        elif debug == 1:
            print(f"\tGET_Svib: HO Term: f={freqs[idx]:.4f}: {Svib_HO:.4e}, {np.round(weight_HO,4)}")

        ## Combine both terms through weights and divide by number of molecules
        total += (Svib_FR * weight_FR + Svib_HO * weight_HO)/nmol
        if debug > 0: print(f"GET_Svib: {total:.4e} au")
        
    ## Arranges units 
    if   outunits.lower() == 'kj':  total = total*constants.har2kJmol
    elif outunits.lower() == 'eV':  total = total*constants.har2eV
    elif outunits.lower() == 'cm':  total = total*constants.har2cm
    elif outunits.lower() == 'au':  pass
    else: raise ValueError(f"GET_Svib: can't understand desired output units: {outunits}")

    ## Creates data-class object
    new_data = Data("Svib", float(total), outunits.lower(), "scope.thermal_corrections.get_Svib()")
    new_data.add_property("temperature", temp, overwrite=True)

    return new_data

###########
def get_Hvib(freqs: list, temp: float, freq_units: str='au', outunits: str='au', nmol: int=1, debug: int=0):
    # temperature in K
    # function works with freqs in au, so we adapt if needed

    if   freq_units.lower() == 'au':  freqs = np.array(freqs)
    elif freq_units.lower() == 'cm':  freqs = np.array(freqs) * constants.cm2har
    elif freq_units.lower() == 'ev':  freqs = np.array(freqs) * constants.eV2har
    elif freq_units.lower() == 's_1': freqs = np.array(freqs) * constants/har2s_1
    else: raise ValueError("GET_Hvib: can't understand input units of frequencies")
    
    total=0.0
    if debug > 0: print(f"GET_Hvib: Computing Hvib with {len(freqs)} frequencies, and first: {freqs[0]} au")
    for idx, f in enumerate(freqs):
        if f > 0.0:
            if temp > 0: exponential = np.exp(-f/(constants.boltz_au*temp))       # Dimensionless
            else:        exponential = float(0.0)                                 # Dimensionless
            fstterm = f/2.                                           # hartree/molecule
            scnterm = (f*exponential)/(1-exponential)                # hartree/molecule
            total += (fstterm+scnterm)/nmol
            if debug > 0: print(f"GET_Hvib: {idx} {total} {f} {fstterm} {scnterm}")
        else:
            if debug > 0: print(f"GET_Hvib: {idx} {total} {f}")
    if debug > 0: print(f"GET_Hvib: Left loop with {total=} au")

    ## Arranges units 
    if outunits.lower() == 'kj':  total = total*constants.har2kJmol    # kJ/mol

    ## Creates data-class object
    new_data = Data("Hvib", float(total), outunits, "scope.thermal_corrections.get_Hvib()")
    new_data.add_property("temperature", temp, overwrite=True)
    return new_data

def get_Selec(spin_multiplicity, outunits: str='au', nmol: int=1):
    if outunits.lower()     == 'kj': value = float(8.314*np.log(spin_multiplicity)/1000/nmol)
    elif outunits.lower()   == 'au': value = float(8.314*np.log(spin_multiplicity)/constants.har2kJmol/1000/nmol)
    return Data("Selec", value, outunits,  'scope.thermal_corrections.get_Selec()') 

def get_Gibbs(Helec: float, Hvib: float, Selec: float, Svib: float, temp: float):
    return Helec + Hvib - temp*(Svib + Selec)

def find_t12(templist, dGlist: list):
    if type(templist) == range: templist = range2list(templist) 
    if dGlist[0] < 0.0: return None 
    else:
        for idx, g in enumerate(dGlist):
            if g < 0.0: return float(templist[idx])
