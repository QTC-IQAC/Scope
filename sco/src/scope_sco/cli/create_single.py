import os
import sys
from argparse import ArgumentParser
from scope.read_write import load_binary, clear_screen
from scope_sco.sco_classes import System_sco
  
def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Path {path} does not exist!')
    env = load_binary(path)
    if env.object_type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def create_single_parser(subparsers):
    parser = subparsers.add_parser("create_single", help="Creates one SCO System from cell2mol data", description="Creates a SCO System from the crystal subfolders in its source folder")
    parser.add_argument('-n', '--env', type=env_exists, required=True, help='Path to the Environment')
    parser.add_argument('-s', '--system', '--source', required=True, help='System name and folder name inside env.sources_path')
    parser.add_argument('-f', '--force', action='store_true', help='Replaces the entire saved System')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.set_defaults(func=create_single)

def create_single(args):
    clear_screen()

    #########################
    # Loads the Environment #
    #########################
    env = load_binary(args.env)
    print(f"\tEnvironment {env.name} Loaded. ")

    ###############################
    # Defines Overwrite and Debug #
    ###############################
    if args.verbose: debug = 1  
    else:            debug = 0
    if args.force:   overwrite = 1
    else:            overwrite = 0
    ###############################

    if debug > 0: print(f"\tArguments_parsed: {args}")

    system_name = args.system
    if system_name in ('', '.', '..') or os.path.basename(system_name) != system_name:
        sys.exit('-s expects a System name, e.g. -s ABITEM')
    new_sys = System_sco(system_name)
    new_sys.set_paths_from_environment(env, create_folders=False, debug=debug)
    if os.path.isfile(new_sys.system_file) and not overwrite:
        print(f"\tSkipping {new_sys.system_file} (use -f to replace the System)")
        return
    if not os.path.isdir(new_sys.sources_path):
        sys.exit(f'Source folder not found: {new_sys.sources_path}')

    print(f"\tLoading cell2mol sources from {new_sys.sources_path}")
    new_sys.load_multiple_cell2mol_folders(new_sys.sources_path, debug=debug)
    if not new_sys.sources:
        sys.exit(f'No cell2mol sources found in {new_sys.sources_path}')
    worked1 = new_sys.set_reference_cells(debug=debug)
    worked2 = new_sys.set_reference_molecs(debug=debug)
    if not (worked1 and worked2):
        sys.exit(f'Could not build the HS and LS references for {system_name}')
    print(f"\tSaving System {new_sys.name} in {new_sys.system_file}")
    new_sys.create_folders()
    new_sys.save()
    if not new_sys.should_save:
        new_sys.clear_status_markers()

if __name__ == "__main__":
    create_single()
