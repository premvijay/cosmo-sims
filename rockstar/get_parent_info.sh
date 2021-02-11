#!/bin/bash
export simnm=${simnm:-L200_N256_Cp18} rund=${rund:-r1} L=${L:-200}
export dir_sim_data=/scratch/cprem/sims/$simnm/$rund/
cd ${dir_sim_data}/halos_rs
parallel -j 26 "find_parents out_{}.list $L > out_wp_{}.list" ::: {0..50}
# find_parents out