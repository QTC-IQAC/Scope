import os
from argparse import ArgumentParser

#### INPUT ####
debug           = 1
path            = "blabal"
overwrite       = False
###############

##### ENVIRONMENT #####
env = load_binary(env_path)



def path_exist(path):
    if not os.path.exists(path):
        raise ValueError(f'Path {path} does not exist!')
    
    return path

def class_type(type):
    from Scope_New.Classes_System import system
    from Scope_New.Spin_Crossover.SCO_Classes import sco_system

    types = dict(
        regular = system,
        sco     = sco_system,
    )

    if not type in types.keys():
        raise ValueError('...')

    return types[type]


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-p', '--path',  type=path_exist, help='missatge de help pel path.')
    parser.add_argument('-t', '--type',  type=class_type)
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    

    return parser.parse_args()

def main():

    args = parse_args()

    for name in sorted(os.listdir(args.path)):
        if os.path.isdir(path):
                new_sys.load_multiple_cell2mol_folders(source_path, debug=debug)
            if worked1 and worked2: 
                new_sys.create_folders()
                new_sys.save()
