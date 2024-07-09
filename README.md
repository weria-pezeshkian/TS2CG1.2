<img width="261" alt="Screenshot 2023-03-03 at 12 31 28" src="https://user-images.githubusercontent.com/47776510/222710874-316a7a7a-5401-4e1c-8082-e786fbb5f206.png">

# TS2CG version 1.2

TS2CG converts triangulated surfaces (TS) to coarse-grained membrane models for molecular
simulation.
It also works as a backmapping algorithm from dynamically triangulated surfaces simulations to CG molecular dynamics simulations or
to take electron microscopy tomography data and build structures for molecular dynamics simulations.

## State

**NOTE**: this version is under development and might be associated with errors and bugs.
Please use the previous version if you are not in direct contact with the developers.
Previous version can be found
https://github.com/marrink-lab/TS2CG1.1

## Installation

### Prerequisites

TS2CG is implemented in C++ and includes three separate scripts. *Pointillism* (PLM) and *Membrane
Builder* (PCG) and *Solvate* (SOL).

The minimum installation requirements are:
- Up-to-date *C++ compilers*.
- *Python* 3.6 or later.
- *CMake* version 3.10 or later.

TS2CG builds with the CMake build system, requiring at least version 3.10. You can check whether
CMake is installed, and what version it is, with

```console
cmake --version.
```

If you need to install CMake, then first check whether your platform’s package management system
provides a suitable version, or visit the [CMake installation](https://cmake.org/resources/) page
for pre-compiled binaries, source code and installation instructions.

### Install _TS2CG_
#### From PyPi
```console
pip3 install TS2CG
```
#### Directly from GitHub
```console
pip3 install git+https://github.com/jan-stevens/TS2CG1.2@python_wrapper
```
#### From source
```console
git clone https://github.com/weria-pezeshkian/TS2CG1.2
cd TS2CG1.2
python3 -m venv venv && source venv/bin/activate # Not required, but often convenient.
pip3 install .
```

## Usage

```console
TS2CG --help
```

### Pointillism
...
```console
TS2CG PLM -h
```

### Membrane Builder
...
```console
TS2CG PCG -h
```

### Solvate
...
```console
TS2CG SOL -h
```

## Quick references
[About TS2CG](https://github.com/weria-pezeshkian/TS2CG1.2/wiki/About-TS2CG) \
[Tutorials](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/558-9-ts2cg) \
[New Updates](https://github.com/weria-pezeshkian/TS2CG1.2/wiki/Updates-of-this-version)

