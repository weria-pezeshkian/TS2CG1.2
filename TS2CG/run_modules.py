import argparse
from argparse import RawTextHelpFormatter

import subprocess
from pathlib import Path
from importlib.metadata import version

from TS2CG.PointUpdaterClass.PointUpdaterClass import *


def run_python_module(module_name: str, args: list[str]):
    """
    run the specified python module with given arguments.
    """
    if module_name == 'PUC':
        PointUpdaterClass.from_args(args)
    else:
        print(f"Unknown Python module: {module_name}")


def run_cpp_module(module_name: str, args: list[str]) -> int:
    """
    run the specified cpp module with given arguments.
    """
    current_dir = Path(__file__).parent
    module_binary_path = current_dir / module_name

    # verify if module binary exists
    if not module_binary_path.exists():
        raise filenotfounderror(f"binary {binary_name} not found at {binary_path}")

    # run the module and catch possible errors
    result = subprocess.run([module_binary_path] + args)
    if result.stderr:
        print("errors:", result.stderr, file=subprocess.sys.stderr)

def main():
    """
    main entry point for the TS2CG command-line interface.
    """

    # submodules are handled differently based on their type
    # define the c++ submodules
    cpp_modules = ['SOL', 'PLM', 'PCG']

    # define the python based modules
    python_modules = ['PUC']

    # parse arguments before a calling module
    description = """
    ==========================================================================================================

    -- TS2CG:  converts triangulated surfaces (ts) to coarse-grained membrane models
    -- Version:  version 1.2
    -- Groningen Biomolecular Sciences and Biotechnology Institute and Zernike Institute for Advanced Materials,
    -- University of Groningen, Groningen, Netherlands
    -- For more information contact Weria Pezeshkian: weria.pezeshkian@gmail.com
    -- citation: Pezeshkian, W. et al. Backmapping triangulated surfaces to coarse-grained membrane models. Nat Commun 11, 2296 (2020).

    ==========================================================================================================
    """
    parser = argparse.ArgumentParser(
        description=description,
        prog='TS2CG',
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        'module',
        choices=cpp_modules + python_modules,
        help='choice of which module to run'
    )

    parser.add_argument(
        'args',
        nargs=argparse.REMAINDER,
        help='arguments for the chosen module'
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f'%(prog)s {version("ts2cg")}'
    )

    args = parser.parse_args()

    # call the right subroutine based on the module type
    if args.module in cpp_modules:
        return run_cpp_module(args.module, args.args)
    elif args.module in python_modules:
        return run_python_module(args.module, args.args)

if __name__ == '__main__':
    main()

