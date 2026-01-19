import argparse
from scope_sco.cli.create_many   import create_many_parser
from scope_sco.cli.create_single import create_single_parser
#from scope_sco.cli.get_t12       import get_t12_parser

def main():
    parser = argparse.ArgumentParser(prog="scope_sco")
    subparsers = parser.add_subparsers(dest="command", required=True)

    create_many_parser(subparsers)
    create_single_parser(subparsers)
    #get_t12_parser(subparsers)

    args = parser.parse_args()
    args.func(args)