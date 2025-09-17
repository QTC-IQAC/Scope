import os
from argparse import ArgumentParser
from Scope.Classes_Environment import environment
from Scope.Read_Write import load_binary
from Scope.Utils.Run_Job import run_job

def path_exists(path):
    if not os.path.exists(path):
        raise ValueError(f'Path {path} does not exist!')
    return path

def env_exists(path):
    if not os.path.isfile(path):
        raise ValueError(f'Environment Path {path} does not exist!')
    env = load_binary(path)
    if env.type == "environment":
        return path
    else:
        raise ValueError(f'Path {path} is not an Environment binary file!')

def parse_args():
    parser = ArgumentParser(prog="scope_run_job", description="Runs a job for a given system")
    parser.add_argument('-n', '--env_path',   type=env_exists,   help='Path to the Environment. Script will load Source data in env.sources_path')
    parser.add_argument('-s', '--sys_path',   type=path_exists,  help='Path to the System File')
    parser.add_argument('-j', '--job_path',   type=path_exists,  help='Path to the Job File')
    parser.add_argument('-v', '--verbose',    help='If true, will print debug information', action='store_true')
    parser.add_argument('-e', '--errors',     help,'If true, will automatically handle some common errors', action='store_true')
    return parser.parse_args()

def main():
    ### Environment is initiated
    args = parse_args()
    summary = ""
    report = run_job(args.sys_path, args.job_path, args.env_path, handle_errors=args.errors, debug=args.verbose)
    if type(report) == str: summary += report

    print('-----SUMMARY-----')
    print(summary)
    print('-------END-------')

if __name__ == "__main__":
    main()
