import argparse
from scope.cli.config import config_parser
from scope.cli.run    import run_parser

def main():
    parser = argparse.ArgumentParser(prog="scope")
    subparsers = parser.add_subparsers(dest="command", required=True)

    config_parser(subparsers)
    run_parser(subparsers)

    args = parser.parse_args()
    args.func(args)