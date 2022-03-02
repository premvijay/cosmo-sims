#!/bin/bash
# grid=1024
nmpil=1
ncpus=32
echo $2
eval "$2"

printf "PERIODIC                                     # enables periodic boundary condistions
NTYPES=6                                     # number of particle types
GADGET2_HEADER                               # allows reading of snapshots with GADGET-2/3 header format
SELFGRAVITY                                   # switch to enable self-gravity of particles (typically always on)
PMGRID=$grid                                   # basic mesh size for TreePM calculations
NSOFTCLASSES=3                               # number of different softening classes
DOUBLEPRECISION=1                             # if activated and set to 1, use double precision internally, for 2 use mixed precision, otherwise single precision
DOUBLEPRECISION_FFTW                         # if set, carries out FFTs in double precision, otherwise single precision
POWERSPEC_ON_OUTPUT                          # computes a matter power spectrum when the code writes a snapshot output
NUMBER_OF_MPI_LISTENERS_PER_NODE=$nmpil           # set such that the number of MPI-ranks per node and listener is maller than MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY
MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY=32    # default is 64, but can also be set to 32\n" > $1

if [ "$fofsub" = "yes" ]
then
printf "FOF\n" >> $1
fi

if [ "$bary" = "yes" ]
then
printf "COOLING
STARFORMATION\n" >> $1
fi

if [ "$ngenic" = "yes" ]
then
printf "NGENIC=$grid
NGENIC_2LPT
CREATE_GRID\n" >> $1
if [ "$bary" = "yes" ]
then
printf "GENERATE_GAS_IN_ICS\n" >> $1
fi
fi

