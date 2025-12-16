# src/scope/cli/main.py
import argparse

def main():
    parser = argparse.ArgumentParser(prog="scope")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # scope sco ...
    p_sco = subparsers.add_parser("sco", help="Spin crossover tools")
    sco_sub = p_sco.add_subparsers(dest="sco_cmd", required=True)

    p_single = sco_sub.add_parser("single", help="Single SCO computation")
    p_single.set_defaults(func=run_sco_single)

    p_many = sco_sub.add_parser("many", help="Multiple SCO computations")
    p_many.set_defaults(func=run_sco_many)

    args = parser.parse_args()
    args.func(args)
