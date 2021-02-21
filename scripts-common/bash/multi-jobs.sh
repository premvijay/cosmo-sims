#!/bin/bash
boxsize=${1:-200}
Npart=${2:-512}
cosmology=${3:-p18}
export simnm="L${boxsize}_N${Npart}_C${cosmology}" rund=${4:-r1} seed=${6-8899}
z_in=${5:-24}

dir_root=$HOME/cosmo-sims/
cd $dir_root

create-siminfo.sh $boxsize $Npart $cosmology $rund $z_in

gadget4/compile.sh

seeds=(4444 5555 6666 7777)
runds=(r4 r5 r6 r7)

for i in {0..3};
do
export rund=${runds[i]} seed=${seeds[i]};
create-siminfo.sh $boxsize $Npart $cosmology $rund $z_in

jidmono=$(qsub monofonic/comp_ics.pbs -v "simnm=$simnm,rund=$rund,seed=$seed")

jidgad=$(qsub gadget4/runsim.pbs -v "simnm=$simnm,rund=$rund"  -W depend=afterok:${jidmono%.*})

jidvr=$(qsub velociraptor/runstf.pbs -v "simnm=$simnm,rund=$rund,space=6d"  -W depend=afterok:${jidgad%.*})
jidrs=$(qsub rockstar/runrstar.pbs -v "simnm=$simnm,rund=$rund"  -W depend=afterok:${jidgad%.*})

jidtf=$(qsub treefrog/runtree.pbs -v "simnm=$simnm,rund=$rund,space=6d"  -W depend=afterok:${jidvr%.*});

done;