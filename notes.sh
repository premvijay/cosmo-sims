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



An important finding is that Inclusive Halo mass calculation. It mainly affect 200m, 200c, 500c but not Mvir
On the other hand Mvir affected by ncellfac
# reference frame origin either centre of mass or particle with minimum potential.