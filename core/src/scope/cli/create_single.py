import os
import sys
from argparse import ArgumentParser
from scope.read_write import load_binary, clear_screen
from scope.classes_system import System

def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Path {path} does not exist!')
    env = load_binary(path)
    if env.object_type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def create_single_parser(subparsers):
    parser = subparsers.add_parser("create_single", help="Creates a System from one XYZ source", description="Creates a System from one XYZ file in its source folder")
    parser.add_argument("-n", "--env", type=env_exists, required=True, help='Path to the Environment')
    parser.add_argument('-s', '--system', required=True, help='System name, associated with a folder name inside env.sources_path')
    parser.add_argument('-S', '--source', required=True, help="XYZ filename inside the System's source folder")
    parser.add_argument('-f', '--force', action='store_true', help='Replaces the entire saved System')
    parser.add_argument('-v', '--verbose', action='store_true')
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
    if args.verbose:
        debug = 1
    else:
        debug = 0
    if args.force:
        overwrite = 1
    else:
        overwrite = 0
    ###############################

    if debug > 0:
        print(f"\tArguments_parsed: {args}")

    system_name = args.system
    if system_name in ('', '.', '..') or os.path.basename(system_name) != system_name:
        sys.exit('-s expects a System name, e.g. -s WATER')
    if os.path.basename(args.source) != args.source or not args.source.endswith('.xyz'):
        sys.exit('-S expects an XYZ filename inside the System folder, e.g. -S water1.xyz')
    new_sys = System(system_name)
    new_sys.set_paths_from_environment(env, create_folders=False, debug=debug)
    if os.path.isfile(new_sys.system_file) and not overwrite:
        print(f"\tSkipping {new_sys.system_file} (use -f to replace the System)")
        return

    source_path = os.path.join(new_sys.sources_path, args.source)
    if not os.path.isfile(source_path):
        sys.exit(f'XYZ source file not found: {source_path}')

    print(f"\tLoading Source From Path {source_path}: ")
    new_sys.load_single_xyz(source_path, debug=debug)
    print(f"\tSaving System {new_sys.name} in {new_sys.system_file}")
    new_sys.create_folders()
    new_sys.save()
    if not new_sys.should_save:
        new_sys.clear_status_markers()


if __name__ == "__main__":
    create_single()
