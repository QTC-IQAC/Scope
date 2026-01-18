# src/scope/cli/main.py
#import sys
#import argparse
#from scope.cli.config import main as run_config
#from scope.cli.run    import main as run_task

#def main():
#    parser = argparse.ArgumentParser(prog="scope")
#    subparsers = parser.add_subparsers(dest="command", required=True)
#
#    # scope config
#    p_config = subparsers.add_parser("config", help="Configure SCOPE environment")
#    p_config.set_defaults(func=run_config)
#
#    # scope run
#    p_run = subparsers.add_parser("run", help="Run Tasks")
#    p_run.set_defaults(func=run_task)
#
#    # Adapt args so subcommand only sees theirs
#    args, remaining = parser.parse_known_args()
#
#    # If user asked for help on a subcommand, forward it from command to subcommand
#    if remaining and remaining[0] in ("-h", "--help"):  sys.argv = [sys.argv[0], "--help"]
#    else:                                               sys.argv = [sys.argv[0]] + remaining
#
#    args.func()

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