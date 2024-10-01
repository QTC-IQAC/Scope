import numpy as np

boltz_au=0.00000316682968         ## hartree/K
boltz_J=1.380649e-23              ## J/K
speed_light=299792458             ## meter/second
planck_hs=1.5198300452227972e-16  ## hartree/second
planck_Js=6.626076e-34            ## J/second == Kg·m2/second
bohr2angs = 0.529177
angs2bohr = 1/bohr2angs
#################

### Conversion Factor. Always to be applied as multiplications

har2kJmol=2625.5002
har2cm=219474.6
har2eV=27.211399
eV2har=0.03674930495
cm2har=0.000004556335281
ry2har=0.5


### Derived Factors. Always to be applied as multiplications

cm2s_1=100*speed_light
eV2s_1=eV2har*har2cm*cm2s_1
har2s_1=har2cm*cm2s_1

