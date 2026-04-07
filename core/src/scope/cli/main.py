import argparse
from scope.cli.create_many      import create_many_parser
from scope.cli.create_single    import create_single_parser
from scope.cli.config           import config_parser
from scope.cli.run              import run_parser
from scope.cli.set_path         import set_path_parser

def main():
    parser = argparse.ArgumentParser(prog="scope")
    subparsers = parser.add_subparsers(dest="command", required=True)

    create_many_parser(subparsers)
    create_single_parser(subparsers)
    config_parser(subparsers)
    run_parser(subparsers)
    set_path_parser(subparsers)

    args = parser.parse_args()
    args.func(args)
