#!/bin/bash
boxsize=${1:-200}
Npart=${2:-512}
cosmology=${3:-p18}
export simnm="L${boxsize}_N${Npart}_C${cosmology}" rund=${4:-r1} seed=${6-8899}
z_in=${5:-24}

dir_root=$HOME/cosmo-sims/
cd $dir_root

create-siminfo.sh $1 $2 $3 $4 $5

jidmono=$(qsub monofonic/comp_ics.pbs -v "simnm=$simnm,rund=$rund,seed=$seed")

gadget4/compile.sh
jidgad=$(qsub gadget4/runsim.pbs -v "simnm=$simnm,rund=$rund"  -W depend=afterok:${jidmono%.*})

jidvr=$(qsub velociraptor/runstf.pbs -v "simnm=$simnm,rund=$rund,space=6d"  -W depend=afterok:${jidgad%.*})
# qsub velociraptor/runstf.pbs -v "simnm=$simnm,rund=$rund,space=3d"
jidrs=$(qsub rockstar/runrstar.pbs -v "simnm=$simnm,rund=$rund"  -W depend=afterok:${jidgad%.*})