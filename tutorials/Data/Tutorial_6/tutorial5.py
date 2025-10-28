import Scope
from Scope import *

tutorial_folder = os.path.abspath('.')

## In this example, we create three systems, each with 4 .xyz sources.
from Scope.Classes_System import system
water_sys   = system("water")
benz_sys    = system("benzene")
form_sys    = system("formaldehyde")

## We set the paths. See Tutorial_4 for more options and information
list_of_systems = []
for tsys in [water_sys, benz_sys, form_sys]:
    tsys.sources_path = f"{tutorial_folder}Sources/{tsys.name}/"
    tsys.calcs_path   = f"{tutorial_folder}Calcs/{tsys.name}/"
    tsys.sys_path     = f"{tutorial_folder}Systems/{tsys.name}/"
    tsys.sys_file     = f"{tutorial_folder}Systems/{tsys.name}/{tsys.name}.npy"
    list_of_systems.append(tsys)

env = load_binary("scope_env_tutorials.npy")
print(env)

for tsys in list_of_systems: 
    tsys.set_paths_from_environment(env)
    tsys.load_multiple_xyz(tsys.sources_path)
    print(tsys)
    tsys.save()
