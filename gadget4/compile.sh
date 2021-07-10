#!/bin/bash
dir_root=$HOME/cosmo-sims/gadget4/
dir_sim=$dir_root/$simnm/
mkdir -p $dir_sim/$rund
$dir_root/gen-config.sh $dir_sim/$rund/Config.sh `cat $dir_sim/$rund/compile.info`
gadget_srcdir=$HOME/tools/gadget4/
# rm -f $dir_sim/$rund/Gadget4
rm -r -f $dir_sim/$rund/build/
echo "make -j 30 -C $gadget_srcdir DIR=$dir_sim/$rund/ &> $dir_sim/$rund/compile.log"
make -j 30 -C $gadget_srcdir DIR=$dir_sim/$rund/ &> $dir_sim/compile.log

