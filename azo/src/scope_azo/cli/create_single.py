import os
import sys
from argparse import ArgumentParser
from scope.read_write      import load_binary, clear_screen
from scope_azo.azo_classes import System_azo

def sanitize_smiles(smiles_str: str) -> str:
    """
    Removes duplicated slash escaping introduced by the shell.

    Parameters:
        smiles_str (str):              Raw SMILES string.

    Returns:
        str | None: Sanitized SMILES string.
    """
    if smiles_str is None:
        return None
    # Replace double forward slashes with a single one
    cleaned = smiles_str.replace('//', '/')
    # Replace literal double backslashes with a single one
    cleaned = cleaned.replace(r'\\', '\\')
    return cleaned
  
def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Path {path} does not exist!')
    env = load_binary(path)
    if env.object_type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def create_single_parser(subparsers):
    parser = subparsers.add_parser("create_single", help="Creates an AZO System from SMILES", description="Creates an AZO System with cis and trans sources from a stereodefined azo SMILES string")
    parser.add_argument('-n', '--env', type=env_exists, required=True, help='Path to the Environment')
    parser.add_argument('-s', '--system', '--name', required=True, help='Name of the System to create')
    parser.add_argument('-smi', '--smiles', type=sanitize_smiles, required=True, help='Azo SMILES string, enclosed in quotes')
    parser.add_argument('-ts', '--create_ts', action='store_true', help='Also generates transition-state and triplet guesses')
    parser.add_argument('-f', '--force', action='store_true', help='Replaces the entire saved System')
    parser.add_argument('-v',    '--verbose',  action='store_true')
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
        sys.exit('-s expects a System name, e.g. -s azobenzene')
    system_file = os.path.join(env.systems_path, system_name, f'{system_name}.npy')
    if os.path.isfile(system_file) and not overwrite:
        print(f"\tSkipping {system_file} (use -f to replace the System)")
        return

    print(f"\tCreating System {system_name}")
    new_sys = System_azo(system_name, args.smiles, debug=debug)
    print(f"\tCreation TRANS isomer")
    if new_sys.create_trans(debug=debug) is None:
        sys.exit(f'Could not build the trans isomer for {system_name}')
    print(f"\tCreation CIS isomer")
    if new_sys.create_cis(debug=debug) is None:
        sys.exit(f'Could not build the cis isomer for {system_name}')
    if args.create_ts: 
        print(f"\tCreation of potential TS structures")
        new_sys.create_ts(debug=debug)
    new_sys.set_paths_from_environment(env, create_folders=True, debug=debug)
    print(f"\tSaving System {new_sys.name} in {new_sys.system_file}")
    new_sys.save()
    if not new_sys.should_save:
        new_sys.clear_status_markers()

if __name__ == "__main__":
    create_single()
