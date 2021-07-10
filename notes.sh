create-siminfo.sh 40 512 p18 r1 24
simnm=L40_N512_Cp18 rund=r1 gadget4/compile.sh
qsub monofonic/comp_ics.pbs -v "simnm=L40_N512_Cp18,rund=r1,seed=8899"
qsub gadget4/runsim.pbs -v "simnm=L40_N512_Cp18,rund=r1"
qsub velociraptor/runstf.pbs -v "simnm=L40_N512_Cp18,rund=r1,space=6d"
qsub velociraptor/runstf.pbs -v "simnm=L40_N512_Cp18,rund=r1,space=3d"
qsub rockstar/runrstar.pbs -v "simnm=L40_N512_Cp18,rund=r1"
cd rockstar/
simnm=L40_N512_Cp18 rund=r1 ./gen-param.sh $simnm/param_r1.txt
rockstar -c L40_N512_Cp18/param_r1.txt
simnm=L40_N512_Cp18 rund=r1 L=40 ./get_parent_info.sh

rund=r12 z_in=24 bary=yes softlen1=0.02 tstep=0.02 seed=222222 compgad=1 multi-jobs.sh 200 256 p18

rund=r21 z_in=24 bary=yes softlen1=0.02 tstep=0.02 seed=111111 compgad=1 ngenic=yes multi-jobs.sh 200 256 p18

rund=r1 z_in=24 bary=no softlen1=0.02 tstep=0.1 outlston=1 seed=1111 compgad=1 ngenic=no multi-jobs.sh 600 1024 p18

# An important finding is that Inclusive Halo mass calculation. Changing from 0 to 3 affect 200m, 200c, 500c but not Mvir
# However changing from 2 to 3 affect only Mvir (hmf shifts up for 2)
# On the other hand decreasing ncellfac tilts the hmf to be lower at high masses.
# reference frame origin either centre of mass or particle with minimum potential.