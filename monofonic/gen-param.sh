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

printf '%s\n' "[setup]
# number of grid cells per linear dimension for calculations = particles for sc initial load
GridRes         = $Npart

# length of the box in Mpc/h
BoxLength       = $boxsize

# starting redshift
zstart          = $z_in

# order of the LPT to be used (1,2 or 3)
LPTorder        = 3

# also do baryon ICs?
DoBaryons       = $bary

# do mode fixing à la Angulo&Pontzen
DoFixing        = yes

# invert phases (for paired simulations)
DoInversion     = no

# particle load, can be 'sc' (1x), 'bcc' (2x) or 'fcc' (4x) (increases number of particles by factor!) or 'glass'
ParticleLoad    = sc

[cosmology]
## transfer = ... specifies the Einstein-Boltzmann plugin module
# transfer        = eisenstein   # Eisenstein&Hu fitting formula

transfer        = CAMB_file    # CAMB file to be specified as 'transfer_file = ...'
transfer_file   = $trans_file

# transfer        = CLASS          # CLASS module (if enabled in CMake file)

# ztarget         = 2.5            # target redshift for CLASS module, output at ztarget will be back-scaled to zstart

Omega_m         = $Om0
Omega_L         = $OmL
Omega_b         = $Omb
H0              = $H0
sigma_8         = $sig8
nspec           = $ns
# A_s             = 2.148752e-09 # can use instead of sigma_8

ZeroRadiation   = true # For Back-scaling: set to false if your simulation code can deal with Omega_r!=0

[random]
## generator = ... specifies the random field generator plugin module
generator       = NGENIC
seed            = $seed

[execution]
NumThreads      = 32


#########################################################################################

[output]
## format = .... specifies the output plugin module

format          = gadget_hdf5
filename        = $dir_ics/monofonic_${z_in}" > $1

