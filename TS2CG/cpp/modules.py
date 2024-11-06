from pathlib import Path
from typing import Union, List
from .base import CPPModule

class PCG(CPPModule):
    """Wrapper for PCG (Point Cloud Generator) C++ module"""
    def __init__(self):
        super().__init__("PCG")

    def add_arguments(self, parser):
        """Add a catch-all argument for PCG"""
        parser.add_argument('args', nargs='*', help='Arguments to pass to PCG binary')

    def __call__(self, *args: str):
        """Run PCG with given arguments"""
        return self.run(list(args))

class PLM(CPPModule):
    """Wrapper for PLM (Pointillism) C++ module"""
    def __init__(self):
        super().__init__("PLM")

    def add_arguments(self, parser):
        """Add a catch-all argument for PLM"""
        parser.add_argument('args', nargs='*', help='Arguments to pass to PLM binary')

    def __call__(self, *args: str):
        """Run PLM with given arguments"""
        return self.run(list(args))

class SOL(CPPModule):
    """Wrapper for SOL (Solvate) C++ module"""
    def __init__(self):
        super().__init__("SOL")

    def add_arguments(self, parser):
        """Add a catch-all argument for SOL"""
        parser.add_argument('args', nargs='*', help='Arguments to pass to SOL binary')

    def __call__(self, *args: str):
        """Run SOL with given arguments"""
        return self.run(list(args))
