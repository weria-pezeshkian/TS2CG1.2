#!/bin/bash

ts2cg="/coarse/melanie/TS2CG_tutorial/TS2CG/TS2CG1.1/"
files="/coarse/melanie/TS2CG_tutorial/TS2CG/TS2CG1.1/Tutorials/files"
pore_radius=2

#==============================================================================================
# increasing the number of vertices
# placing exclusions
#==============================================================================================

$ts2cg/PLM  -TSfile  $files/Sphere.tsi  -bilayerThickness  0  -Mashno 3 -resizebox


### chosing vertices for exclusion placement
### two pores along each axis
### at minimum and maximum coordinates

##  | all the lines between 'vertex' and 'triangle'|
val=($(awk '/triangle/{f=0} f; /vertex/{f=1}' extended.tsi |
	## | for each of the three coordinat coloumns find the vertex with minimum and maxiumum coordinate |
	awk '!init && NF>0 {nf=NF; init=1; for (i = 2; i <= nf; ++i) max[i] = min[i] = $i; minvertex[i] = maxvertex[i] = $1} \
	init==1 && NF>0 {for (i = 2; i <= nf; ++i) {if (max[i] < $i) {max[i] = $i; maxvertex[i] = $1} else if (min[i] > $i) {min[i] = $i; minvertex[i] = $1}}}\
	END{for (i = 2; i <= nf; ++i)  print(minvertex[i], maxvertex[i])}'))
#echo ${val[*]}

### generate file with exclusions
cp extended.tsi extended_exclusions.tsi
printf "exclusion %19s\n" "${#val[@]}" >> extended_exclusions.tsi
for ((i=0;i<${#val[@]};i++)); do
	printf "%8s %10s %6s\n" $i ${val[i]} $pore_radius >> extended_exclusions.tsi
done

#===============================================================================================
# use PLM on modified .tsi file to read the exclusions
# use PCG to place lipids (here only POPC) except for excluded vertices
#===============================================================================================

$ts2cg/PLM  -TSfile  extended_exclusions.tsi  -bilayerThickness  3.8 -rescalefactor 2 2 2
$ts2cg/PCG  -str  $files/input.str  -Bondlength  0.2  -LLIB  $files/Martini3.LIB  -defout  system

#===============================================================================================
# generate a position restrains will using MDAnlysis
# a position restraints file with all the lipid tails collapsed in the vesicle center
# is necessary to keep the porses open during th equilibration
#===============================================================================================

python3 $files/gen_posres_vesicle.py  system.gro tail_posres.gro

## write a proper topology file

> topol.top

for filename in "$files/itp/martini3"/*.itp; do
	printf "#include \"%s\"\n" "$filename" >> topol.top
done

cat system.top >> topol.top

#===============================================================================================
# Vacuum energy minimisation
#===============================================================================================

## soft-core
gmx grompp -f $files/mdp/vesicle/em_1.mdp -c system.gro -r system.gro -p topol.top -o em_1.tpr
gmx mdrun -v -deffnm em_1

## regular
gmx grompp -f $files/mdp/vesicle/em_2.mdp -c em_1.gro -r system.gro -p topol.top -o em_2.tpr
gmx mdrun -v -deffnm em_2


#===============================================================================================
# solvation and another energy minimization
#==============================================================================================

## using homemade solvation script
$files/SOL -in em_2.gro -tem $files/water.gro -o SOL.gro -Rcutoff 0.32
cat info.txt >> topol.top

gmx grompp -f $files/mdp/vesicle/em_2.mdp -c SOL.gro -r SOL.gro -p topol.top -o em_SOL.tpr -maxwarn 1
gmx mdrun -v -deffnm em_SOL

#==============================================================================================
# 6-step equilibration
#==============================================================================================

## Note: mdp files need to be adapted for other system compositions (different lipids, ions ...)
## providing an index file should be the easiest

cntmax=6

for (( cnt=1; cnt<=6; cnt++));do
        pcnt=$(( $cnt - 1 ))
	if [ $cnt -eq 1 ]; then
		gmx grompp -f $files/mdp/vesicle/eq_$cnt.mdp -c em_SOL.gro -r tail_posres.gro -p topol.top -o eq_$cnt.tpr -maxwarn 2
	else
		gmx grompp -f $files/mdp/vesicle/eq_$cnt.mdp -c eq_$pcnt.gro -r tail_posres.gro -p topol.top -o eq_$cnt.tpr -maxwarn 2
	fi
	gmx mdrun -v -deffnm eq_$cnt
done

rm *#*

#==============================================================================================
# production run
#==============================================================================================

gmx grompp -f $files/mdp/vesicle/md.mdp -c eq_$cntmax.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
