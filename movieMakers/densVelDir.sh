#!/bin/bash

# A script to make density and director field movies for a given simulation. 
# Call this script from the directory containing the simulation output.

#args
L=$1 # system size
frameNo=$2 # max number of frames to render
dirEvery=$3 # show a director every dirEvery steps
defectPath=$4 # path to defect topology file

# script
mkdir dens
mkdir dir
mkdir vel
cd dens
python3 /mnt/c/Work/Code/Github/MPCD_ActivityAnalysis/movieMakers/density2Danimated.py ../coarsegrain.dat $L $L 1 0 $frameNo z 0 equal 1 1
cd ../vel
python3 /mnt/c/Work/Code/Github/MPCD_ActivityAnalysis/movieMakers/flowFieldAnimated.py ../flowfield.dat $L $L 1 0 $frameNo $dirEvery $dirEvery z equal 1 $defectPath
cd ../dir
python3 /mnt/c/Work/Code/Github/MPCD_ActivityAnalysis/movieMakers/orientationField2Danimated.py ../directorfield.dat $L $L 1 0 $frameNo $dirEvery $dirEvery z 0.45 equal 1 $defectPath