import os
from argparse import ArgumentParser
from Scope_New.Read_Write import load_binary
from Scope_New.Spin_Crossover.SCO_Classes import sco_system

#### INPUT ####
debug       = 1

this_type
overwrite   = False

def env_exist(path):
    if not os.path.isfile(path):
        raise ValueError(f'Path {path} does not exist!')
    env = load_binary(path)
    if env.type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-n', '--env',     type=env_exist,      help='Path to the Environment')
    parser.add_argument('-f', '--force',   action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser.parse_args()

def main():

    args = parse_args()
    print(args)

    #########################
    # Loads the Environment #
    #########################
    env = load_binary(args.name)
    print(f"\tEnviroment {env.name} Loaded. ")
    print(f"\tLoading Sources From Path (env.sources_path): ")

    ###############################
    # Defines Overwrite and Debug #
    ###############################
    if args.verbose: debug = 1  
    else:            debug = 0
    if args.force:   overwrite = 1
    else:            overwrite = 0

    for name in sorted(os.listdir(env.sources_path)):
        if os.path.isdir(path):
            new_sys = sco_system(name, env)
            new_sys.load_multiple_cell2mol_folders(env.sources_path, debug=debug)
            worked1 = new_sys.set_reference_cells(overwrite=overwrite, debug=debug)
            worked2 = new_sys.set_reference_molecs(overwrite=overwrite, debug=debug)
            if worked1 and worked2: 
                print(f"\tCreation of system {new_sys.name} worked. Saving sys_file here: {new_sys.sys_file}. Folders will be created if necessary:")
                new_sys.create_folders()
                new_sys.save()
