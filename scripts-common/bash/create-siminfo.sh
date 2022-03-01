#!/bin/bash
source ~/.bashrc

echo "Create sim config files for given Box size, number of particles along one axis and cosmology."
echo "Also provide rundir name and initial redshift."

boxsize=${1:-200}
Npart=${2:-512}
cosmology=${3:-p18}
bary=${bary:-no}
ngenic=${ngenic:-no}
outlston=${outlston:-0}

export simnm="L${boxsize}_N${Npart}_C${cosmology}" rund=${rund:-r1} seed=${seed-8899}
if [ "$bary" = yes ]; then export simnm=${simnm}_bar; fi

export softlen1=${softlen1:-0.0065} tstep=${tstep:-0.01}

export z_in=${z_in:-24}

echo Using values: $boxsize, $Npart, $cosmology, $rund, $z_in, $simn

# cat pigoome;

dir_root=$HOME/cosmo-sims/

dir_camb=$dir_root/camb/
dir_mono=$dir_root/monofonic/
dir_music=$dir_root/music/
dir_gad=$dir_root/gadget4/
dir_sw=$dir_root/swift/
dir_vel=$dir_root/velociraptor/

mkdir -p $dir_mono/$simnm/
mkdir -p $dir_music/$simnm/
mkdir -p $dir_gad/$simnm/$rund
mkdir -p $dir_sw/$simnm/$rund

#GADGET configuration info
echo "grid=$((2*Npart));nmpil=1;ncpus=$ncpus;bary=$bary;ngenic=$ngenic;fofsub=$fofsub" > $dir_gad/$simnm/$rund/compile.info

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

  lh1)
    Om0=0.3391 ;OmL=0.6609 ;Omb=0.05172 ;h_par=0.6682 ;ns=0.9931 ;sig8=0.8162
    ;;

  lh2)
    Om0=0.3503 ;OmL=0.6497 ;Omb=0.0505 ;h_par=0.7522 ;ns=0.9703 ;sig8=0.7802
    ;;

  lh3)
    Om0=0.3055 ;OmL=0.6945 ;Omb=0.04682 ;h_par=0.685 ;ns=0.9475 ;sig8=0.7562
    ;;

  lh4)
    Om0=0.3167 ;OmL=0.6833 ;Omb=0.04988 ;h_par=0.7186 ;ns=0.9855 ;sig8=0.7682
    ;;

  lh5)
    Om0=0.2607 ;OmL=0.7393 ;Omb=0.04805 ;h_par=0.601 ;ns=0.9323 ;sig8=0.8042
    ;;

  lh6)
    Om0=0.3279 ;OmL=0.6721 ;Omb=0.04743 ;h_par=0.6178 ;ns=0.9399 ;sig8=0.8282
    ;;

  lh7)
    Om0=0.2831 ;OmL=0.7169 ;Omb=0.04866 ;h_par=0.6346 ;ns=1.0007 ;sig8=0.8402
    ;;

  lh8)
    Om0=0.3615 ;OmL=0.6385 ;Omb=0.04621 ;h_par=0.7354 ;ns=0.9627 ;sig8=0.8642
    ;;

  lh9)
    Om0=0.2943 ;OmL=0.7057 ;Omb=0.05111 ;h_par=0.7018 ;ns=0.9551 ;sig8=0.7922
    ;;

  lh10)
    Om0=0.2719 ;OmL=0.7281 ;Omb=0.04927 ;h_par=0.6514 ;ns=0.9779 ;sig8=0.8522
    ;;
esac

echo "$Om0, $OmL, $Omb, $h_par, $ns, $sig8" > $dir_camb/${cosmology}.info   # To compile gadget
#To generate parameter file for GADGET-4
echo "Om0=$Om0;OmL=$OmL;Omb=$Omb;h_par=$h_par;ns=$ns;sig8=$sig8;boxsize=$boxsize;Npart=$Npart;z_in=$z_in;cosmo=$cosmology;softlen1=$softlen1;timestep=$tstep;outlston=$outlston;bary=$bary;ngenic=$ngenic;seed=$seed" > $dir_gad/$simnm/param_$rund.info
#To generate parameter file for swift
echo "Om0=$Om0;OmL=$OmL;Omb=$Omb;h_par=$h_par;ns=$ns;sig8=$sig8;boxsize=$boxsize;Npart=$Npart;z_in=$z_in;cosmo=$cosmology;softlen1=$softlen1;timestep=$tstep;outlston=$outlston;bary=$bary;ngenic=$ngenic;seed=$seed" > $dir_sw/$simnm/param_$rund.info
#To generate parameter file for monofonIC
echo "Om0=$Om0;OmL=$OmL;Omb=$Omb;h_par=$h_par;ns=$ns;sig8=$sig8;boxsize=$boxsize;Npart=$Npart;z_in=$z_in;cosmo=$cosmology;bary=$bary;seed=$seed" > $dir_mono/$simnm/param_$rund.info
#To generate parameter file for MUSIC
echo "Om0=$Om0;OmL=$OmL;Omb=$Omb;h_par=$h_par;ns=$ns;sig8=$sig8;boxsize=$boxsize;Npart=$Npart;z_in=$z_in;cosmo=$cosmology;bary=$bary;seed=$seed" > $dir_music/$simnm/param_$rund.info
# echo "simnm="L${boxsize}_N${Npart}_C${cosmology}" rund=${4:-r1}"

# $dir_gad/compile.sh