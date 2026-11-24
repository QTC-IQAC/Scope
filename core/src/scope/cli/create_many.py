import os
from argparse               import ArgumentParser
from scope.read_write       import load_binary, clear_screen
from scope.classes_system   import System

def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Path {path} does not exist!')
    env = load_binary(path)
    if env.object_type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def create_many_parser(subparsers):
    parser = subparsers.add_parser("create_many", help="Creates Systems from XYZ files", description="Creates one System per source folder containing XYZ files")
    parser.add_argument("-n", "--env", type=env_exists, required=True, help='Path to the Environment. Script will load all sources in env.sources_path')
    parser.add_argument('-f', '--force', action='store_true', help='Rebuilds each existing System from the XYZ files in its source folder')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.set_defaults(func=create_many)


def create_many(args):
    clear_screen()

    #########################
    # Loads the Environment #
    #########################
    env = load_binary(args.env)
    print(f"\tEnviroment {env.name} Loaded. ")
    print(f"\tLoading All Sources From Path {env.sources_path}: ")

    print(f"\t------------------------------------------------")
    print(f"\tThe script expects the following folder structure")
    print(f"\t\t{env.sources_path}/")
    print(f"\t\t\t...System_1/")
    print(f"\t\t\t\t...*.xyz")
    print(f"\t\t\t...System_2/")
    print(f"\t\t\t\t...*.xyz")
    print(f"\t------------------------------------------------")

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
        print(f"Arguments_parsed: {args}")

    for name in sorted(os.listdir(env.sources_path)):
        sys_path = os.path.join(env.sources_path, name)
        ## In sources_path, I expect a folder for each system.
        if os.path.isdir(sys_path):
            new_sys = System(name)
            new_sys.set_paths_from_environment(env, create_folders=False, debug=debug)
            if os.path.isfile(new_sys.system_file) and not overwrite:
                print(f"\tSkipping {new_sys.system_file} (use -f to replace the System)")
                continue
            new_sys.load_multiple_xyz(new_sys.sources_path, debug=debug)
            if not new_sys.sources:
                print(f"\tSkipping {sys_path} (no XYZ files)")
                continue
            print(f"\tSaving System {new_sys.name} in {new_sys.system_file}")
            new_sys.create_folders()
            new_sys.save()
            if not new_sys.should_save:
                new_sys.clear_status_markers()


if __name__ == "__main__":
    create_many()
