#!/bin/bash
source ~/.bashrc
cosmo_param=p18
z_in=24
boxsize=200
Npart=512

export simnm="cdm_${boxsize}M_${Npart}_${cosmo_param}" rund=r1

dir_root=$HOME/cosmo-sims/

dir_camb=$dir_root/camb/
dir_mono=$dir_root/monofonic/

mkdir -p $dir_mono/$simnm/

dir_gad=$dir_root/gadget4/
# dir_sim=$dir_gad/$simnm/
mkdir -p $dir_gad/$simnm/

echo "grid=$((2*Npart));nmpil=2;ncpus=32" > $dir_gad/$simnm/compile.info
case $cosmo_param in
  p18)
    echo "0.3063375, 0.6936625, 0.0484103, 0.6781, 0.9677, 0.815" > $dir_camb/${cosmo_param}.info
    echo "Om0=0.3063375;OmL=0.6936625;Omb=0.0484103;h_par=0.6781;boxsize=$boxsize" > $dir_gad/$simnm/param.info
    echo "Om0=0.3063375;OmL=0.6936625;Omb=0.0484103;h_par=0.6781;ns=0.9677;sig8=0.815;boxsize=$boxsize;Npart=$Npart;z_in=$z_in" > $dir_mono/$simnm/param.info
    ;;

  wmap7)
    echo "0.276, 0.724, 0.045, 0.7, 0.961, 0.811" > $dir_camb/${cosmo_param}.info
    echo "Om0=0.276;OmL=0.724;Omb=0.045;h_par=0.7;boxsize=$boxsize" > $dir_gad/$simnm/param.info
    echo "Om0=0.276;OmL=0.724;Omb=0.045;h_par=0.7;ns=0.961;sig8=0.811;boxsize=$boxsize;Npart=$Npart;z_in=$z_in" > $dir_mono/$simnm/param.info
    ;;
esac

$dir_gad/compile.sh