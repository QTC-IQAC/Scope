import numpy as np

###########
## Notes ##
###########
# Joule is kg·m2·s-2
# 'cm' in variable names below, actually means cm-1

######################
# Physical Constants #
######################
hbar                =1                       ## Atomic Units
boltz_au            =0.00000316682968        ## hartree/K
boltz_J             =1.380649e-23            ## J/K
speed_light         =299792458               ## meter/second
planck_hs           =1.5198300452227972e-16  ## hartree·second
planck_Js           =6.626076e-34            ## J·second == Kg·m2/second
R_J                 =8.3144                  ## J/(K·mol)                   Universal Gas Constant
avogadro            =6.02214076e23           ## Avogadro's number
elem_charge         =1.602176634e-19         ## Coulomb  
electron_mass       =9.1093837139e-31        ## Kg
epsilon_0           =8.8541878188e-12        ## F/m = C^2⋅s^2/(kg−1⋅m−3)    Vacuum permittivity

######################
# Conversion Factors #
######################
# Always to be applied as multiplications
bohr2angs       = 0.529177
angs2bohr       = 1/bohr2angs

# Energy
har2kJmol       = 2625.5002
har2cm          = 219474.6
har2eV          = 27.211399
eV2har          = 0.03674930495
cm2har          = 0.000004556335281
ry2har          = 0.5
amu2au          = 1822.888486  # 1 amu in atomic units
kcal2kJmol      = 4.184
kJmol2kcal      = 1/kcal2kJmol

### Derived Factors. Always to be applied as multiplications
cm2s_1          = 100 * speed_light
eV2s_1          = eV2har * har2cm * cm2s_1
har2s_1         = har2cm * cm2s_1

### Pressure
atm2kpa         = 101.325
