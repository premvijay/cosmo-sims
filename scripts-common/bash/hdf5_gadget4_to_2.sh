#!/bin/bash
source ~/.bashrc
conda activate conforg
mypython=`which python`
echo $1
parallel -j 30  "$mypython -m hdf5_gadget4_to_2 {}" ::: $1 
# cat $PBS_NODEFILE
# echo hi 
# parallel -j 30  --sshloginfile $PBS_NODEFILE  "conda activate conforg; $mypython -m hdf5_gadget4_to_2 {}" ::: $1 