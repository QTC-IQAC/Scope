import numpy as np
from scope                            import constants
from scope.classes_data               import *
from scope.operations.dicts_and_lists import range2list

###########
def get_Svib(freqs: list, temp: float, freq_units: str='au', outunits: str='au', typ: str='default', FR_cutoff: int=100, FR_alpha: int=4, nmol: int=1, debug: int=0):
    
    ## Temperature must be provided in K
    ## Function works with freqs in au, except the free-rotor part, which works in s-1 units

    bav=1.0000E-44 # Kg*m**2
    
    if typ.lower() == 'free-rotor' or typ.lower() == 'fr':              
        if   freq_units.lower() == 'cm':  FR_cutoff=FR_cutoff*constants.cm2s_1
        elif freq_units.lower() == 'ev':  FR_cutoff=FR_cutoff*constants.eV2s_1
        elif freq_units.lower() == 'au' : FR_cutoff=FR_cutoff*constants.har2s_1
        elif freq_units.lower() == 's_1': pass
        else: 
            print("Svin: can't understand input units of frequencies")
            sys.exit()
        
    # Converts Frequencies to s-1
    if   freq_units.lower() == 'au':  freqs = np.array(freqs)
    elif freq_units.lower() == 'cm':  freqs = np.array(freqs) * constants.cm2har
    elif freq_units.lower() == 'ev':  freqs = np.array(freqs) * constants.eV2har
    elif freq_units.lower() == 's_1': freqs = np.array(freqs) * constants/har2s_1
    else: raise ValueError("GET_Svib: can't understand input units of frequencies")
        
    total=0.0
    for idx, f in enumerate(freqs):   # Freqs in s-1
        if f > 0.0: 
            ## Free Rotor Term
            if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
                mu=constants.planck_Js/(8*(np.pi)**2*f)                                      # Kg·m2
                mu_prime=(mu*bav)/(mu+bav)                                                   # Kg·m2
                a=(8*(np.pi)**3)*mu_prime*constants.boltz_J*temp/(constants.planck_Js**2)    # Dimensionless
                b=np.sqrt(a)                                                                 # Dimensionless
                c=np.log(b)                                                                  # Dimensionless
                Svib_FR=constants.boltz_au*(c+0.5)                                           # Hartree/molecule*K
                #print(f, FR_cutoff, FR_alpha)
                weight=(1/(1+(FR_cutoff/f)**FR_alpha))                                       # Dimensionless
                Svib_FR=(1-weight)*Svib_FR/nmol                                              # Hartree/molecule*K
            else: 
                Svib_FR=0.0
                weight=0.0
    
            ## Harmonic Oscillator Term
            f = f/constants.har2s_1  # Now frequency in har
            exponential_pos=np.exp(f/(constants.boltz_au*temp))                              # Dimensionless
            exponential_neg=np.exp(-f/(constants.boltz_au*temp))                             # Dimensionless
            fstterm=f/(temp*(exponential_pos-1))                                             # Hartree/molecule*K
            scnterm=-constants.boltz_au*np.log(1-exponential_neg)                            # Hartree/molecule*K
            
            if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
                weight=1/(1+(FR_cutoff/f)**FR_alpha)
            else:
                weight=1.0
            Svib_HO = (fstterm+scnterm)*weight/nmol
    
            ## Sums both Contributions
            total += (Svib_FR + Svib_HO)

    if debug > 0: print(f"GET_Svib: Left main loop with {total=} au")
        
    ## Arranges units 
    if outunits.lower() == 'kj':  total = total*constants.har2kJmol
    ## Creates data-class object
    new_data = Data("Svib", float(total), outunits, "scope.thermal_corrections.get_Svib()")
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
