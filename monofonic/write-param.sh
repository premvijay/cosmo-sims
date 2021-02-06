#!/bin/bash
dir_root=$HOME/cosmo-sims/monofonic/
dir_sim=$dir_root/$simnm/
$dir_root/gen-param.sh $dir_sim/param_$rund.txt  `cat $dir_sim/param_$rund.info`