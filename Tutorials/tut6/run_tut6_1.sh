#!/bin/bash

## Tutorial 6: Fixed Shapes
## 6-1. 1D Fourier shape POPC bilayer
## This script generates a 1D Fourier-shaped POPC membrane and runs the TS2CG outputs using GROMACS 
#===============================================================================================
## Function to display help text
#===============================================================================================

function show_help {
    echo "Usage: script.sh [-h]"
    echo "No options will only execute the TS2CG part of the Tutorial."
    echo "Options:"
    echo "  -h         Display this help message and exit."
    exit
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h) 
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

#===============================================================================================
files="files"  
#===============================================================================================
# use PCG to place lipids 
#===============================================================================================

TS2CG PCG -str input.str -Bondlength 0.2 -LLIB $files/Martini3.LIB -defout system -function analytical_shape || { echo "Error: TS2CG PCG command failed"; exit 1; }

#===============================================================================================
## Edit the topology file
#===============================================================================================

> topol.top
printf "#include \"$files/martini3/martini_v3.0.4.itp\"\n" >> topol.top
printf "#include \"$files/martini3/martini_v3.0_phospholipids.itp\"\n" >> topol.top

if [ -f system.top ]; then
    cat system.top >> topol.top
else
    echo "Warning: system.top does not exist. Skipping append."
fi

#===============================================================================================
# Prompt the user to execute the GROMACS part
#===============================================================================================

read -p "Do you want to execute the GROMACS simulation part? (yes/no): " user_input

if [[ "$user_input" == "yes" ]]; then
    if which gmx &> /dev/null; then
        gmx_exec=gmx
    elif which gmx_d &> /dev/null; then
        gmx_exec=gmx_d
    elif which gmx_mpi &> /dev/null; then
        gmx_exec=gmx_mpi
    else
        echo "GROMACS executable not found. Please install it and try again."
        exit 1
    fi

     if [ -d "output_1" ]; then
     mv output_1 "#output_1"
     fi

     mkdir -p output_1

     set -o errexit
     set -o nounset

    #===============================================================================================
    # Vacuum energy minimization
    #===============================================================================================
    ## soft-core
    $gmx_exec grompp -f $files/mdp/1DFourierShape/em_1.mdp -c system.gro -r system.gro -p topol.top -o output_1/em_1.tpr
    $gmx_exec mdrun -s output_1/em_1.tpr -v -deffnm output_1/em_1

    ## regular
    $gmx_exec grompp -f $files/mdp/1DFourierShape/em_2.mdp -c output_1/em_1.gro -r output_1/em_1.gro -p topol.top -o output_1/em_2.tpr
    $gmx_exec mdrun -s output_1/em_2.tpr -v -deffnm output_1/em_2

    #===============================================================================================
    # Vacuum equilibration
    #===============================================================================================

    $gmx_exec grompp -f $files/mdp/1DFourierShape/eq_v.mdp -c output_1/em_2.gro -p topol.top -r output_1/em_2.gro -o output_1/eq_v.tpr
    $gmx_exec mdrun -v -s output_1/eq_v.tpr -deffnm output_1/eq_v

else
    echo "Skipping the GROMACS part."
fi
