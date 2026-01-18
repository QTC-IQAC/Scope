from argparse import ArgumentParser
from scope.classes_environment import Environment
from scope.read_write import clear_screen

#def parse_args():
#    #parser = ArgumentParser(prog="scope config", description="Configure a new scope Environment")
#    parser = ArgumentParser(description="Configure a new scope Environment")
#    parser.add_argument('-n', '--name',  required=True, help='Suffix to be added to the name of the environment. File will be called: scope_env_NAME.npy}')
#    #parser.add_argument('-n', '--name',  type=str, required=True, help='Suffix to be added to the name of the environment. File will be called: scope_env_NAME.npy}')
#    return parser.parse_args()

def config_parser(subparsers):
    parser = subparsers.add_parser("config",help="Configure SCOPE environment",description="Configure a new scope Environment")
    parser.add_argument("-n", "--name",required=True,help="Suffix to be added to the environment name")
    parser.set_defaults(func=config)

def config(args):
    clear_screen()

    ### Environment is initiated
    env = Environment(args.name)

    print("")
    print(f" Creating Environment scope_env_{env.name}.npy. We will now:")
    print(f" 1) Set the Paths for the SOURCE, CALCULATIONS and SYSTEMS folders associated with this environment") 
    print(f" Then, if the Computer has SLURM or SGE scheduler, we will:")
    print(f" 2) Set Software Modules for Quantum Espresso and Gaussian 16")
    print(f" 3) Set Available Queues/Partitions")
    print("")

    env.set_paths()
    env.set_software()
    env.set_queues()

    env.save()

    print("")
    print(f" --> Environment Created:")
    print(env)
    

    print("")
    print(f" --> Environment Saved in {env.filepath} as binary")
    print("")
    print(f" --> You can overwrite your selections by loading the binary and:")
    print(f"     1) To change Paths:                         env.set_paths():")
    print(f"     2) To change Software Modules:              env.set_software():")
    print(f"     3) To change Available Queues/Partitions:   env.set_queues():")
    #print(f"     4) To change Storage Path:                  env.set_storage_path():")
    print("")
    if env.scheduler != 'local': 
        print(f" --> Proceeding to test the scheduler commands:")
        print("")
        env.test_scheduler()
        print("")
        print(f" --> If any of those checks failed, please contact me at sergi.vela@iqac.csic.es")

if __name__ == "__main__":
    main()
