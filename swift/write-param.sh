#!/bin/bash
dir_root=$HOME/cosmo-sims/swift/
dir_sim=$dir_root/$simnm/
mkdir -p $dir_sim/$rund
$dir_root/gen-param.sh $dir_sim/$rund/param.yml `cat $dir_sim/param_$rund.info`