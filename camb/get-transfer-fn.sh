#!/bin/bash
source ~/.bashrc
dir_camb=/mnt/home/student/cprem/cosmo-sims/camb/
info_file=$dir_camb/${cosmology}.info
trans_file=$dir_camb/${cosmology}_transfer_out.dat
if [ -f "$trans_file" ]; then
    echo "transfer camb file exist already"
else
    (conda activate conforg && python $dir_camb/comp-transfer-fn.py --params "`cat $info_file`" --cambfile $trans_file)
fi