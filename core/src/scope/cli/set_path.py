import os
from argparse import ArgumentParser
from scope.classes_environment import Environment
from scope.read_write import load_binary, clear_screen

def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Environment Path {path} does not exist!')
    env = load_binary(path)
    if env.object_type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def set_path_parser(subparsers):
    parser = subparsers.add_parser("set_path",help="Sets Current Directory as System Main Path",description="Sets Current Directory as System Main Path")
    parser.add_argument('-n', '--env_path',   type=env_exists,   help='Path to the Environment. Script will load Source data in env.sources_path')
    parser.add_argument('-s', '--sys_name',   help="Name of the System. SCOPE will search in the Environment's Systems folder") 
    parser.add_argument('-q', '--quiet',      help='If true, will not print the progress on screen', action='store_true')
    parser.add_argument('-v', '--verbose',    help='If true, will print extra information on screen', action='store_true')
    parser.set_defaults(func=set_path)

def set_path(args):
    clear_screen()
    
    ###############################
    # Defines Overwrite and Debug #
    ###############################
    if   args.quiet:   debug = 0  
    elif args.verbose: debug = 2
    else:              debug = 1
    ###############################

    ## Convert to absolute paths
    args.env_path = os.path.abspath(args.env_path)
    inp_path = os.path.abspath('.')

    ## Loads Environment
    env           = load_binary(args.env_path)

    ## Checks sys_path
    sys_path      = f"{env.systems_path}{args.sys_name}/{args.sys_name}.npy"
    if not os.path.isfile(sys_path):
        raise ValueError(f'System binary file not found in the expected path {sys_path}!')

    ## Prepare System Paths 
    if debug > 0: print(f"Loading System binary in: {sys_path}")
    sys = load_binary(sys_path)
    if debug > 0: print(f"Setting paths for System: {sys.name}")
    if debug > 0: print(f"Current directory is: {inp_path}")
    sys.set_main_path(inp_path, debug=debug)
    if debug > 0: print(f"Updated System in: {sys.system_file}")
    sys.save(sys.system_file)

if __name__ == "__main__":
    set_path()
