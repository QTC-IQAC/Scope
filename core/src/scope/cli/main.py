# src/scope/cli/main.py
import sys
import argparse
from scope.cli.config import main as run_config
from scope.cli.run    import main as run_task

def main():
    parser = argparse.ArgumentParser(prog="scope")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # scope config
    p_config = subparsers.add_parser("config", help="Configure SCOPE environment")
    p_config.set_defaults(func=run_config)

    # scope run
    p_run = subparsers.add_parser("run", help="Run a computation")
    p_run.set_defaults(func=run_task)

    # Adapt args so subcommand only sees theirs
    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + remaining

    args.func()
