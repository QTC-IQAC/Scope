def get_vibS(freqs: list, temp: float, freq_units: str='au', outunits: str='au', typ: str='default', FR_cutoff: int=100, FR_alpha: int=4, nmol: int=1):
    # temperature most be provided in K   
    #bav=float((1.0000E-44*(1/4.359748E-18)))                # From Kg*m**2 to Hartree*S**2
    bav=1/(1.0000E-44/Constants.planck_Js)                   # From Kg*m**2 to s-1

    if typ.lower() == 'free-rotor' or typ.lower() == 'fr':              
        if freq_units.lower() == 'cm': FR_cutoff=FR_cutoff*Constants.cm2har
        elif freq_units.lower() == 'ev': FR_cutoff=FR_cutoff*Constants.eV2har
        
    if temp == 10: print("temp, idx, f, Svib_FR, Svib_HO, total, weight_HO")

    freqs_adapted = []
    if freq_units.lower() != 'au':
        for f in freqs:
            if freq_units.lower() == 'cm': freqs_adapted.append(f*Constants.cm2har)
            elif freq_units.lower() == 'ev': freqs_adapted.append(f*Constants.eV2har)
            else: print("vibS: can't understand input units of frequencies")
    else: freqs_adapted = freqs.copy()

    total=0.0
    for idx, f in enumerate(freqs_adapted):
        if f > 0: 
            ## Free Rotor Term
            if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
                mu=Constants.planck_hs/(8*(np.pi)**2*f)                                      # s-1
                mu_prime=(mu*bav)/(mu+bav)                                                   # s-1
                a=(8*(np.pi)**3)*mu_prime*Constants.boltz*temp/(Constants.planck_hs**2)      # Dimensionless
                b=np.sqrt(a)                                                                 # Dimensionless
                c=np.log(b)                                                                  # Dimensionless
                Svib_FR=Constants.boltz*(c+0.5)                                              # Hartree/molecule*K
                #print(f, FR_cutoff, FR_alpha)
                weight=(1/(1+(FR_cutoff/f)**FR_alpha))                                       # Dimensionless
                Svib_FR=(1-weight)*Svib_FR                                                   # Hartree/molecule*K
            else: 
                Svib_FR=0.0
                weight=0.0

            ## Harmonic Oscillator Term
            exponential_pos=np.exp(f/(Constants.boltz*temp))
            exponential_neg=np.exp(-f/(Constants.boltz*temp))
            fstterm=f/(temp*(exponential_pos-1))
            scnterm=-Constants.boltz*np.log(1-exponential_neg)

            if typ.lower() == 'free-rotor' or typ.lower() == 'fr':
                weight=1/(1+(FR_cutoff/f)**FR_alpha)
            else:
                weight=1.0
            Svib_HO = (fstterm+scnterm)*weight

            ## Sums both Contributions
            total += (Svib_FR + Svib_HO) #*Constants.Avogadro/nmol
            #if temp == 10: print(f"{temp} {idx} {freqs[idx]} {Svib_FR:.3e} {Svib_HO:.3e} {total:.3e} {weight}")
    
    ## Arranges units
    if outunits.lower() == 'kj': total = total*Constants.har2kJmol
    return total

##########################
def get_vibH(freqs: list, temp: float, freq_units: str='au', outunits: str='au'):
    # temperature in K   
    freqs_adapted = []
    if freq_units.lower() != 'au':
        for f in freqs:
            if freq_units.lower() == 'cm': freqs_adapted.append(f*Constants.cm2har)
            elif freq_units.lower() == 'ev': freqs_adapted.append(f*Constants.eV2har)
            else: print("vibS: can't understand input units of frequencies")
    else: freqs_adapted = freqs.copy()
    
    total=0.0
    for f in freqs_adapted:        
        exponential = np.exp(-f/(Constants.boltz*temp))
        fstterm = f/2.
        scnterm = (f*exponential)/(1-exponential)
        total += (fstterm+scnterm)
    if outunits.lower() == 'kj': total = total*Constants.har2kJmol
    return total
