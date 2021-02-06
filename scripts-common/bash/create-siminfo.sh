#!/bin/bash
source ~/.bashrc

echo "Create sim config files for given Box size, number of particles along one axis and cosmology."
echo "Also provide rundir name and initial redshift."

boxsize=${1:-200}
Npart=${2:-512}
cosmology=${3:-p18}

export simnm="L${boxsize}_N${Npart}_C${cosmology}" rund=${4:-r1}

z_in=${5:-24}

echo Default: $boxsize, $Npart, $cosmology, $rund, $z_in, $simn

# cat pigoome;

dir_root=$HOME/cosmo-sims/

dir_camb=$dir_root/camb/
dir_mono=$dir_root/monofonic/
dir_gad=$dir_root/gadget4/
dir_vel=$dir_root/velociraptor/

mkdir -p $dir_mono/$simnm/
mkdir -p $dir_gad/$simnm/

echo "grid=$((2*Npart));nmpil=2;ncpus=32" > $dir_gad/$simnm/compile.info
case $cosmology in
  p18)
    Om0=0.3063375;OmL=0.6936625;Omb=0.0484103;h_par=0.6781;ns=0.9677;sig8=0.815
    ;;

  p13)
    Om0=0.315;OmL=0.685;Omb=0.049;h_par=0.673;ns=0.96;sig8=0.829
    ;;
  
  w7)
    Om0=0.276;OmL=0.724;Omb=0.045;h_par=0.7;ns=0.961;sig8=0.811
    ;;
esac

echo "$Om0, $OmL, $Omb, $h_par, $ns, $sig8" > $dir_camb/${cosmology}.info
echo "Om0=$Om0;OmL=$OmL;Omb=$Omb;h_par=$h_par;boxsize=$boxsize;Npart=$Npart;z_in=$z_in;cosmo=$cosmology;" > $dir_gad/$simnm/param_$rund.info
echo "Om0=$Om0;OmL=$OmL;Omb=$Omb;h_par=$h_par;ns=$ns;sig8=$sig8;boxsize=$boxsize;Npart=$Npart;z_in=$z_in;cosmo=$cosmology;" > $dir_mono/$simnm/param_$rund.info

echo "simnm="L${boxsize}_N${Npart}_C${cosmology}" rund=${4:-r1}"

# $dir_gad/compile.sh