#!/bin/bash
# softlen1=0.0065
# timestep=0.01

# z_in=24
# Om0=0.3063375
# OmL=0.6936625
# Omb=0.0484103
# h_par=0.6781
# boxsize=200.0
echo $2
eval "$2"



dir_root=$HOME/cosmo-sims/swift/
dir_sim=$dir_root/$simnm/
dir_sim_data=/scratch/cprem/sims/$simnm/$rund/
dir_snap=$dir_sim_data/snaps_sw/
dir_ics=$dir_sim_data/ics/

mkdir -p $dir_snap

# conda activate conforg

# a_in=$(printf %.5f `echo "1/($z_in+1)" | bc -l`)
a_in=$(python -c "print (1.0/(1+$z_in))")
# $(echo "scale 4;3*2.5" |bc)
Ocdm=$(python -c "print ($Om0-$Omb)")
hMsun10=$(python -c "print (1.98848e43/$h_par)")
hMpc=$(python -c "print (3.08567758e24/$h_par)")

softlen=$(python -c "print(float($boxsize)/$Npart/30)")
# softlen2=$softlen1

# conda deactivate

printf '%s\n' "#
# Define some meta-data about the simulation
MetaData:
  run_name:   ${simnm}_$rund  # The name of the simulation. This is written into the snapshot headers.

# Define the system of units to use internally. 
InternalUnitSystem:
  UnitMass_in_cgs:     $hMsun10    # 10^10 M_sun
  UnitLength_in_cgs:   $hMpc # 1 Mpc
  UnitVelocity_in_cgs: 1e5   	     # 1 km/s
  UnitCurrent_in_cgs:  1   	     # Amperes
  UnitTemp_in_cgs:     1   	     # Kelvin

# Cosmological parameters
Cosmology:
  h:              $h_par        # Reduced Hubble constant
  a_begin:        $a_in     # Initial scale-factor of the simulation
  a_end:          1.0           # Final scale factor of the simulation
  Omega_cdm:        $Ocdm        # C Matter density parameter
  Omega_lambda:   $OmL         # Dark-energy density parameter
  Omega_b:        $Omb     # Baryon density parameter

# Parameters for the self-gravity scheme
Gravity:
  mesh_side_length:              $Npart       # Number of cells along each axis for the periodic gravity mesh.
  eta:                           0.025     # Constant dimensionless multiplier for time integration.
  MAC:                           geometric  # Choice of mulitpole acceptance criterion: 'adaptive' OR 'geometric'.
  # epsilon_fmm:                   0.001     # Tolerance parameter for the adaptive multipole acceptance criterion.
  theta_cr:                      0.7       # Opening angle for the purely gemoetric criterion.
  comoving_DM_softening:         $softlen # Comoving Plummer-equivalent softening length for DM particles (in internal units).
  max_physical_DM_softening:     $softlen    # Maximal Plummer-equivalent softening length in physical coordinates for DM particles (in internal units).
  comoving_baryon_softening:     $softlen # Comoving Plummer-equivalent softening length for baryon particles (in internal units).
  max_physical_baryon_softening: $softlen    # Maximal Plummer-equivalent softening length in physical coordinates for baryon particles (in internal units).
  #softening_ratio_background:    0.04      # Fraction of the mean inter-particle separation to use as Plummer-equivalent softening for the background DM particles.

# Parameters governing the time integration (Set dt_min and dt_max to the same value for a fixed time-step run.)
TimeIntegration:
  dt_min:            1e-7  # The minimal time-step size of the simulation (in internal units).
  dt_max:            $timestep  # The maximal time-step size of the simulation (in internal units).

# Parameters governing the snapshots
Snapshots:
  basename:   snapshot      # Common part of the name of output files.
  subdir:     snapdir         # (Optional) Sub-directory in which to write the snapshots. Defaults to "" (i.e. the directory where SWIFT is run).
  scale_factor_first: $a_in # (Optional) Scale-factor of the first snapshot if cosmological time-integration.
  delta_time: 1.2        # Time difference between consecutive outputs (in internal units)
  invoke_stf: 0           # (Optional) Call VELOCIraptor every time a snapshot is written irrespective of the VELOCIraptor output strategy.
  invoke_fof: 0           # (Optional) Call FOF every time a snapshot is written
  compression: 0          # (Optional) Set the level of GZIP compression of the HDF5 datasets [0-9]. 0 does no compression. The lossless compression is applied to *all* the fields.
  distributed: 0          # (Optional) When running over MPI, should each rank write a partial snapshot or do we want a single file? 1 implies one file per MPI rank.
  UnitMass_in_cgs:     $hMsun10   # (Optional) Unit system for the outputs (Grams)
  UnitLength_in_cgs:   $hMpc  # (Optional) Unit system for the outputs (Centimeters)
  UnitVelocity_in_cgs: 1e5  # (Optional) Unit system for the outputs (Centimeters per second)
  UnitCurrent_in_cgs:  1  # (Optional) Unit system for the outputs (Amperes)
  UnitTemp_in_cgs:     1  # (Optional) Unit system for the outputs (Kelvin)
  output_list_on:      $outlston  # (Optional) Enable the output list
  output_list:         $dir_root/outputs-few.list # (Optional) File containing the output times (see documentation in Parameter File section)

# Parameters governing the conserved quantities statistics
Statistics:
  delta_time:           1.01        # Time between statistics output

# Parameters related to the initial conditions
InitialConditions:
  file_name:  $dir_ics/monofonic_${z_in}_combined.hdf5 # The file to read
  periodic:                    1    # Are we running with periodic ICs?
  generate_gas_in_ics:         0    # (Optional) Generate gas particles from the DM-only ICs (e.g. from panphasia).
" > $1



if [ "$bary" = "yes" ]
then
printf '%s\n' "%----- Star formation
MaxSfrTimescale     1.5         % Gas consumption timescale (multi-phase model)
FactorSN            0.1         % beta, mass fraction of massive stars (multi-phase model)
FactorEVP           1000        % A_0, evaporation parameter (multi-phase model)
TempSupernova       1e+08       % T_SN, effective supernova temperature,sets feedback energy (multi-phase model)
TempClouds          1000        % temperature of cold clouds (multi-phase model)
CritOverDensity     57.7        % overdensity threshold value for cosological sims
CritPhysDensity     0           % critical physical density for star formation (in cm^(-3))
TreecoolFile        $dir_root/TREECOOL" >> $1
fi

if [ "$ngenic" = "yes" ]
then
printf '%s\n' "%----- N-GenIC
NSample                                           $Npart
GridSize                                          $Npart
Seed                                              $seed
SphereMode                                        0
PowerSpectrumType                                 1
ReNormalizeInputSpectrum                          1
PrimordialIndex                                   $ns
ShapeGamma                                        0.21
Sigma8                                            $sig8
PowerSpectrumFile                                 powerspec
InputSpectrum_UnitLength_in_cm                    3.085678e24" >> $1
fi