import argparse
#from scope_azo.cli.create_many   import create_many_parser
from scope_azo.cli.create_single import create_single_parser

def main():
    parser     = argparse.ArgumentParser(prog="scope_azo")
    subparsers = parser.add_subparsers(dest="command", required=True)

    #create_many_parser(subparsers)
    create_single_parser(subparsers)

    args = parser.parse_args()
    args.func(args)
