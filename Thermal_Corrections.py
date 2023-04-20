import numpy as np
from Test_V3 import Constants
from Test_V3.Classes_Data import *
from Test_V3.Other import range2list

###########
def get_Svib(freqs: list, temp: float, freq_units: str='au', outunits: str='au', typ: str='default', FR_cutoff: int=100, FR_alpha: int=4, nmol: int=1):
    
    ## Temperature must be provided in K
    ## Function works with freqs in au, except the free-rotor part, which works in s-1 units

    bav=1.0000E-44 # Kg*m**2
    
    if typ.lower() == 'free-rotor' or typ.lower() == 'fr':              
        if freq_units.lower() == 'cm': FR_cutoff=FR_cutoff*Constants.cm2s_1
        elif freq_units.lower() == 'ev': FR_cutoff=FR_cutoff*Constants.eV2s_1
        elif freq_units.lower() == 'au' : FR_cutoff=FR_cutoff*Constants.har2s_1
        elif freq_units.lower() == 's_1': pass
        else: 
            print("Svin: can't understand input units of frequencies")
            sys.exit()
        
    #if temp == 10: print("temp, idx, f, Svib_FR, Svib_HO, total, weight_HO")
    
    # Converts Frequencies to s-1
    freqs_adapted = []
    for f in freqs:
        if freq_units.lower() == 'cm': freqs_adapted.append(f*Constants.cm2s_1)
        elif freq_units.lower() == 'ev': freqs_adapted.append(f*Constants.eV2s_1)
        elif freq_units.lower() == 'au' : freqs_adapted.append(f*Constants.har2s_1)
        elif freq_units.lower() == 's_1': freqs_adapted.append(f)
        else: 
            print("Svib: can't understand input units of frequencies")
            sys.exit()
        
    total=0.0
    for idx, f in enumerate(freqs_adapted):   # Freqs in s-1
        if f > 0.0: 
            ## Free Rotor Term
            if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
                mu=Constants.planck_Js/(8*(np.pi)**2*f)                                      # Kg·m2
                mu_prime=(mu*bav)/(mu+bav)                                                   # Kg·m2
                a=(8*(np.pi)**3)*mu_prime*Constants.boltz_J*temp/(Constants.planck_Js**2)    # Dimensionless
                b=np.sqrt(a)                                                                 # Dimensionless
                c=np.log(b)                                                                  # Dimensionless
                Svib_FR=Constants.boltz_au*(c+0.5)                                           # Hartree/molecule*K
                #print(f, FR_cutoff, FR_alpha)
                weight=(1/(1+(FR_cutoff/f)**FR_alpha))                                       # Dimensionless
                Svib_FR=(1-weight)*Svib_FR/nmol                                              # Hartree/molecule*K
            else: 
                Svib_FR=0.0
                weight=0.0
    
            ## Harmonic Oscillator Term
            f = f/Constants.har2s_1  # Now frequency in har
            exponential_pos=np.exp(f/(Constants.boltz_au*temp))      # Dimensionless
            exponential_neg=np.exp(-f/(Constants.boltz_au*temp))     # Dimensionless
            fstterm=f/(temp*(exponential_pos-1))                     # Hartree/molecule*K
            scnterm=-Constants.boltz_au*np.log(1-exponential_neg)    # Hartree/molecule*K
            
            if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
                weight=1/(1+(FR_cutoff/f)**FR_alpha)
            else:
                weight=1.0
            Svib_HO = (fstterm+scnterm)*weight/nmol
    
            ## Sums both Contributions
            total += (Svib_FR + Svib_HO)
            #if temp == 10: print(f"{temp} {idx} {freqs[idx]} {Svib_FR:.3e} {Svib_HO:.3e} {total:.3e} {weight}")
        
    ## Arranges units 
    if outunits.lower() == 'kj':  total = total*Constants.har2kJmol
    ## Creates data-class object
    new_data = data("Svib", float(total), outunits, "get_Svib")
    new_data.add_property("temp", temp, overwrite=True)

    return new_data

###########
def get_Hvib(freqs: list, temp: float, freq_units: str='au', outunits: str='au', nmol: int=1):
    # temperature in K
    # function works with freqs in au
    freqs_adapted = []
    if freq_units.lower() != 'au':
        for f in freqs:
            if freq_units.lower() == 'cm': freqs_adapted.append(f*Constants.cm2har)
            elif freq_units.lower() == 'ev': freqs_adapted.append(f*Constants.eV2har)
            elif freq_units.lower() == 's_1': freqs_adapted.append(f*Constants/har2s_1)
            elif freq_units.lower() == 'au' : freqs_adapted.append(f)
            else: print("Svib: can't understand input units of frequencies")
    
    total=0.0
    for f in freqs_adapted:        
        if f > 0.0:
            exponential = np.exp(-f/(Constants.boltz_au*temp))       # Dimensionless
            fstterm = f/2.                                           # hartree/molecule
            scnterm = (f*exponential)/(1-exponential)                # hartree/molecule
            total += (fstterm+scnterm)/nmol

    ## Arranges units 
    if outunits.lower() == 'kj':  total = total*Constants.har2kJmol    # kJ/mol

    ## Creates data-class object
    new_data = data("Hvib", float(total), outunits, "get_Hvib")
    new_data.add_property("temp", temp, overwrite=True)

    return new_data

def get_Selec(spin, outunits: str='au'):
    if   spin == "HS": multiplicity = int(5)
    elif spin == "IS": multiplicity = int(3)
    elif spin == "LS": multiplicity = int(1)
    if outunits.lower()     == 'kj': value = float(8.314*np.log(multiplicity))
    elif outunits.lower()   == 'au': value = float(8.314*np.log(multiplicity)/Constants.har2kJmol)
    return data("Selec", value, outunits,  'get_Selec') 

def get_Gibbs(Helec: float, Hvib: float, Selec: float, Svib: float, temp: float):
    return Helec + Hvib - temp*(Svib + Selec)

#
#def get_dG(templist, dHelec, freqs_HS, freqs_LS, freq_units='cm', outunits='kj', typ='default'):
#
#    assert len(freqs_HS) == len(freqs_LS)
#
#    ## HS ##
#    Hvib_HS = get_Hvib(freqs_HS, templist, freq_units=freq_units, outunits=outunits).value
#    Svib_HS = get_Svib(freqs_HS, templist, freq_units=freq_units, outunits=outunits, typ=typ).value 
#    ## LS ##
#    Hvib_LS = get_Hvib(freqs_LS, templist, freq_units=freq_units, outunits=outunits).value
#    Svib_LS = get_Svib(freqs_LS, templist, freq_units=freq_units, outunits=outunits, typ=typ).value 
#
#    Selec_HS = get_Selec(5, outunits).value
#    Selec_LS = get_Selec(0, outunits).value
#    dSelec = Selec_HS - Selec_LS
#
#    dSvib = []
#    for idx, z in enumerate(zip(Svib_HS,Svib_LS)):
#        dSvib.append(z[0]-z[1])
#
#    dHvib = []
#    for idx, z in enumerate(zip(Hvib_HS,Hvib_LS)):
#        dHvib.append(z[0]-z[1])
#
#    dG = []
#    for idx, t in enumerate(templist):
#        dG.append(dHelec + dHvib[idx] - t*(dSvib[idx] + dSelec)) 
#
#    new_data = data("dG", dG, outunits, "get_dG", notes="trange= "+str(trange)) 
#    return new_data

def find_t12(templist, dGlist: list):
    if type(templist) == range: templist = range2list(templist) 
    #assert len(templist) == len(dGlist)
    if dGlist[0].value < 0.0: return None 
    else:
        for idx, g in enumerate(dGlist):
            if g.value < 0.0: return templist[idx]
