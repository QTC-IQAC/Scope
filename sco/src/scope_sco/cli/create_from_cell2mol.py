import os
import sys
from argparse import ArgumentParser
from scope.read_write import load_binary
from scope_sco.sco_classes import System_sco
  
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

def parse_args():
    parser = ArgumentParser(prog="create_from_cell2mol", description="Creates a SCO system from cell2mol data")
    parser.add_argument('-n', '--env',     type=env_exists,      help='Path to the Environment. Script will load Source data in env.sources_path')
    parser.add_argument('-s', '--source',  type=str,             help='Name of the Source Folder Inside env.sources_path')
    parser.add_argument('-f', '--force',   action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser.parse_args()

def main():

    args = parse_args()

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

    source_path = env.sources_path+args.source 
    system_name = args.source
    exists = path_exists(source_path)

    if not exists:
        print(f"\t{source_path} does not exist: ")
        sys.exit()

    print(f"\tLoading Source From Path {source_path}: ")
    new_sys = System_sco(system_name)
    new_sys.set_paths_from_environment(env)
    new_sys.load_multiple_cell2mol_folders(new_sys.sources_path, debug=debug)
    worked1 = new_sys.set_reference_cells(overwrite=overwrite, debug=debug)
    worked2 = new_sys.set_reference_molecs(overwrite=overwrite, debug=debug)
    if worked1 and worked2: 
        print(f"\tCreation of system {new_sys.name} worked. Saving sys_file here: {new_sys.system_file}. Folders will be created if necessary")
        new_sys.create_folders()
        new_sys.save()

if __name__ == "__main__":
    main()
