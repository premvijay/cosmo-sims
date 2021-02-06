#!/bin/bash
source ~/.bashrc
conda activate conforg
mypython=`which python`
echo $1
cat $PBS_NODEFILE
parallel -j 30  --sshloginfile $PBS_NODEFILE  $mypython -m hdf5_gadget4_to_2 {} ::: $1 