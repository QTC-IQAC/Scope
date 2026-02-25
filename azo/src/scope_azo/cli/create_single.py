import os
import sys
from argparse import ArgumentParser
from scope.read_write      import load_binary, clear_screen
from scope_azo.azo_classes import System_azo
  
def path_exists(path):
    if not os.path.exists(path):
        raise ValueError(f'Path {path} does not exist!')
    return path

def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Path {path} does not exist!')
    env = load_binary(path)
    if env.type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def create_single_parser(subparsers):
    parser = subparsers.add_parser("create_single",help="Creates an AZO system from a smiles string",description="Creates an AZO system from a smiles string")
    parser.add_argument('-n',      '--env', type=env_exists,  help='Path to the Environment. Script will load Source data in env.sources_path')
    parser.add_argument('-s',     '--name', type=str,         help='Name of the system to be created')
    parser.add_argument('-smi', '--smiles', type=str,         help='SMILES string of the trans_isomer of that system')
    parser.add_argument('-f',    '--force', action='store_true')
    parser.add_argument('-v',  '--verbose', action='store_true')
    parser.set_defaults(func=create_single)

def create_single(args):
    clear_screen()

    #########################
    # Loads the Environment #
    #########################
    env = load_binary(args.env)
    print(f"\tEnviroment {env.name} Loaded. ")

    ###############################
    # Defines Overwrite and Debug #
    ###############################
    if args.verbose: debug = 1  
    else:            debug = 0
    if args.force:   overwrite = 1
    else:            overwrite = 0
    ###############################

    if debug > 0: print(f"\tArguments_parsed: {args}")

    system_name = args.name
    system_path = env.sources_path+args.name
    smiles      = args.smiles

    exists = path_exists(system_path)
    if exists and not overwrite:
        print(f"\tSystem {system_name} already exists. Use -f or --force to overwrite it.")
        sys.exit()

    print(f"\tCreating System {system_name}")
    new_sys = System_azo(system_name, smiles)
    new_sys.set_paths_from_environment(env)
    print(f"\tCreation TRANS isomer")
    new_sys.create_trans(debug=debug)
    print(f"\tCreation CIS isomer")
    new_sys.create_cis(debug=debug)
    if args.create_ts: 
        print(f"\tCreation of potential TS structures")
        new_sys.create_ts(debug=debug)
        print(f"\tCreation of system {new_sys.name} worked. Saving sys_file here: {new_sys.system_file}. Folders will be created if necessary")
        new_sys.create_folders()
        new_sys.save()

if __name__ == "__main__":
    create_single()
