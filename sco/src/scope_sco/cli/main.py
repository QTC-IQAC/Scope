# src/scope/cli/main.py
import sys
import argparse
from scope_sco.cli.create_many   import main as create_many
from scope_sco.cli.create_single import main as create_single
from scope_sco.cli.get_t12       import main as get_t12

def main():
    parser = argparse.ArgumentParser(prog="scope_sco")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # scope_sco create_many
    p_many = subparsers.add_parser("create_many", help="Creates a SCO System from a cell2mol cells")
    p_many.set_defaults(func=create_many)

    # scope_sco create_single
    p_single = subparsers.add_parser("create_single", help="Creates a SCO System for all cell2mol cells in env.sources_path")
    p_single.set_defaults(func=create_single)

    # scope_sco get_t12
    p_get = subparsers.add_parser("get_t12", help="Gets the transition temperature (T1/2) for a system using stored data")
    p_get.set_defaults(func=get_t12)

    # Adapt args so subcommand only sees theirs
    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + remaining

    args.func()
