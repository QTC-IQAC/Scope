import os
import json
from argparse import ArgumentParser
from Scope_New.Classes_Environment import environment

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-n', '--name',  type=str, help='Name of the environment.')
    return parser.parse_args()

def main():
    args = parse_args()
    env = environment(args.name)
    env.set_paths()
    env.save()
    env.save_config()

    print("")
    print(f"Environment Created and Saved in {env.filepath}. See details below")
    print("")
    print(env)

if __name__ == "__main__":
    main()
