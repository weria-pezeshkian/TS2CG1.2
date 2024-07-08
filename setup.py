import os
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from pathlib import Path
import multiprocessing


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])

class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        cwd = Path().absolute()
        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.absolute()

        cmake_args = [
            f'-DCMAKE_INSTALL_PREFIX={extdir}',
        ]

        num_cores = multiprocessing.cpu_count()
        build_args = ['--', f'-j{num_cores}']

        os.chdir(build_temp)
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
            self.spawn(['cmake', '--install', '.'])
        os.chdir(cwd)

setup(
    name='TS2CG',
    version='1.2.0',
    packages=find_packages(),
    ext_modules=[CMakeExtension('TS2CG')],
    cmdclass={'build_ext': CMakeBuild},
    package_data={'TS2CG': ['SOL', 'PLM', 'PCG']},
    include_package_data=True,
    install_requires=[],
    entry_points={
        'console_scripts': [
            'TS2CG=TS2CG.wrapper:main',
        ],
    },
)
