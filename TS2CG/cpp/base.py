import subprocess
from pathlib import Path
from typing import List, Union, Optional
import logging
import os
import sys

logger = logging.getLogger(__name__)

class CPPModule:
    """Base class for C++ module wrappers"""

    def __init__(self, binary_name: str):
        self.binary_name = binary_name
        self._binary_path = self._find_binary()

    def _find_binary(self) -> Path:
        """Find the binary in the package directory"""
        package_dir = Path(__file__).parent.parent
        binary_path = package_dir / self.binary_name

        if not binary_path.exists():
            raise FileNotFoundError(
                f"Binary '{self.binary_name}' not found at {binary_path}. "
                "Make sure the package is properly installed with C++ components built."
            )

        if not os.access(binary_path, os.X_OK):
            raise PermissionError(
                f"Binary '{binary_path}' exists but is not executable. "
                "Check file permissions."
            )

        return binary_path

    def run(self, args: List[str], cwd: Optional[Union[str, Path]] = None) -> subprocess.CompletedProcess:
        """Run the C++ binary with given arguments"""
        try:
            # For help commands, we want to show the output and exit
            if '-h' in args or '--help' in args:
                result = subprocess.run(
                    [str(self._binary_path)] + args,
                    cwd=cwd,
                    text=True
                )
                sys.exit(0)

            # For normal commands, handle errors
            result = subprocess.run(
                [str(self._binary_path)] + args,
                cwd=cwd,
                capture_output=True,
                text=True
            )

            # Print stdout if there is any
            if result.stdout:
                print(result.stdout)

            # Check for errors
            if result.returncode != 0:
                if result.stderr:
                    logger.error(result.stderr)
                raise subprocess.CalledProcessError(
                    result.returncode,
                    result.args,
                    result.stdout,
                    result.stderr
                )

            return result

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running {self.binary_name}: {e.stderr}")
            raise RuntimeError(f"Failed to run {self.binary_name}") from e
