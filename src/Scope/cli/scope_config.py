from argparse import ArgumentParser
from ..Classes_Environment import environment

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-n', '--name',  type=str, help='Name of the environment.')
    return parser.parse_args()

def main():
    ### Environment is initiated
    args = parse_args()
    env = environment(args.name)

    print("")
    print(f"\tEnvironment {env.name} Created. We will now do the following actions:")
    print(f"\t 1) Set Paths for the SOURCE, CALCULATIONS and SYSTEM folders associated with this environment") 
    if env.management_type != "local": 
        print(f"\t 2) Set Software Modules for Quantum Espresso and Gaussian 16")
        print(f"\t 3) Set Available Queues/Partitions")

    env.set_paths()
    env.set_software()
    env.set_queues()

    env.save()
    #config_path = env.save_config()

    print("")
    print(f"\tEnvironment Created and Saved in {env.filepath}. See details below")
    print(env)

    print("")
    #print(f"\t A JSON Config File was also saved in {config_path}")
    #print("")
    print(f"\t You can overwrite your selections by loading the binary and:")
    print(f"\t 1) To change Software Modules:              env.set_software():")
    print(f"\t 2) To change Available Queues/Partitions:   env.set_queues():")
    print(f"\t 3) To change Paths:                         env.set_paths():")
    print(f"\t 4) To change Storage Path:                  env.set_storage_path():")

    print("")
if __name__ == "__main__":
    main()
