#import os
#from argparse import ArgumentParser
#from Scope.Read_Write import load_binary
#from Scope.Spin_Crossover.SCO_Classes import sco_system
#
#def env_exists(path):
#    if not os.path.isfile(path):
#        raise ValueError(f'Path {path} does not exist!')
#    env = load_binary(path)
#    if env.type == "environment":
#        return path
#    else:
#        raise ValueError(f'Path {path} is not an Environment binary file!')
#
#def parse_args():
#    parser = ArgumentParser(prog="get_t12", description="Computes the transition temperature T12 for a given branch of a SCO system")
#    parser.add_argument('-n', '--env',     type=env_exists,      help='Path to the Environment. Script will load all sources in env.sources_path')
#    parser.add_argument('-b', '--branch',  type=str)
#    parser.add_argument('-f', '--force',   action='store_true')
#    parser.add_argument('-v', '--verbose', action='store_true')
#    return parser.parse_args()
#
#def main():
#
#    args = parse_args()
#
#    #########################
#    # Loads the Environment #
#    #########################
#    env = load_binary(args.env)
#    print(f"\tEnviroment {env.name} Loaded. ")
#
#    ###############################
#    # Defines Overwrite and Debug #
#    ###############################
#    if args.verbose: debug = 1
#    else:            debug = 0
#    if args.force:   overwrite = 1
#    else:            overwrite = 0
#    ###############################
#
#    if debug > 0: print(f"\tArguments_parsed: {args}")
#
#    source_path = env.sources_path+args.source
#    system_name = args.source
#    exists = path_exists(source_path)
#
#    if not exists:
#        print(f"\t{source_path} does not exist: ")
#        sys.exit()
#
#        new_sys.save()
#
#if __name__ == "__main__":
#    main()
#