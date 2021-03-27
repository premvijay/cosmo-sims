#!/bin/bash
# Om0=0.3063375
# OmL=0.6936625
# Omb=0.0484103
# h_par=0.6781
# sig8=0.815
# ns=0.9677
# z_in=24
# boxsize=200.0
# Npart=512
# seed=8888
# cosmology=p18
echo $2
eval "$2"

dir_sim_data=/scratch/cprem/sims/$simnm/$rund/
dir_ics=$dir_sim_data/ics/

dir_camb=/mnt/home/student/cprem/cosmo-sims/camb/
trans_file=$dir_camb/${cosmo}_transfer_out.dat

mkdir -p $dir_ics

H0=$(python -c "print($h_par*100)")

level=$(python -c "import math; print(int(math.log($Npart,2)))")

printf '%s\n' "[output]
format			= gadget2
filename		= $dir_ics/music_${z_in}
gadget_num_files	= 8

[random]
seed[8]		= $seed
cubesize	= 64

[setup]
boxlength		= $boxsize
zstart			= $z_in
levelmin		= $level
levelmin_TF		= $level
levelmax		= $level
padding			= 4
overlap			= 4
align_top		= no
periodic_TF		= no
baryons			= $bary
use_2LPT		= yes
use_LLA			= no
kspace_TF		= yes

[cosmology]
Omega_m         = $Om0
Omega_L         = $OmL
Omega_b         = $Omb
H0              = $H0
sigma_8         = $sig8
nspec           = $ns
transfer	= camb_file
transfer_file	= $trans_file

[poisson]
kspace          = yes
fft_fine        = no
accuracy        = 1e-5
smoother        = gs
laplace_order   = 6
grad_order      = 6" > $1
