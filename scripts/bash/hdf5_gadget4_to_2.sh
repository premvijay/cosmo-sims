#!/bin/bash
source ~/.bashrc
conda activate conforg
echo $1
parallel -j 25 python -m hdf5_gadget4_to_2 {} ::: $1 