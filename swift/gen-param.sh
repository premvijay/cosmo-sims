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

softlen0=$softlen1
softlen2=$softlen1

dir_root=$HOME/cosmo-sims/gadget4/
dir_sim=$dir_root/$simnm/
dir_sim_data=/scratch/cprem/sims/$simnm/$rund/
dir_snap=$dir_sim_data/snaps_sw/
dir_ics=$dir_sim_data/ics/

mkdir -p $dir_snap

# a_in=$(printf %.5f `echo "1/($z_in+1)" | bc -l`)
a_in=$(python -c "print (1.0/(1+$z_in))")

printf '%s\n' "#
# Define some meta-data about the simulation
MetaData:
  run_name:   Name of the sim in less than 256 characters.  # The name of the simulation. This is written into the snapshot headers.

# Define the system of units to use internally. 
InternalUnitSystem:
  UnitMass_in_cgs:     1.98848e43    # 10^10 M_sun
  UnitLength_in_cgs:   3.08567758e24 # 1 Mpc
  UnitVelocity_in_cgs: 1e5   	     # 1 km/s
  UnitCurrent_in_cgs:  1   	     # Amperes
  UnitTemp_in_cgs:     1   	     # Kelvin

# Cosmological parameters
Cosmology:
  h:              0.6777        # Reduced Hubble constant
  a_begin:        0.0078125     # Initial scale-factor of the simulation
  a_end:          1.0           # Final scale factor of the simulation
  Omega_m:        0.307         # Matter density parameter
  Omega_lambda:   0.693         # Dark-energy density parameter
  Omega_b:        0.0482519     # Baryon density parameter

# Parameters for the self-gravity scheme
Gravity:
  mesh_side_length:              128       # Number of cells along each axis for the periodic gravity mesh.
  eta:                           0.025     # Constant dimensionless multiplier for time integration.
  MAC:                           adaptive  # Choice of mulitpole acceptance criterion: 'adaptive' OR 'geometric'.
  epsilon_fmm:                   0.001     # Tolerance parameter for the adaptive multipole acceptance criterion.
  theta_cr:                      0.7       # Opening angle for the purely gemoetric criterion.
  use_tree_below_softening:      0         # (Optional) Can the gravity code use the multipole interactions below the softening scale?
  allow_truncation_in_MAC:       0         # (Optional) Can the Multipole acceptance criterion use the truncated force estimator?
  comoving_DM_softening:         0.0026994 # Comoving Plummer-equivalent softening length for DM particles (in internal units).
  max_physical_DM_softening:     0.0007    # Maximal Plummer-equivalent softening length in physical coordinates for DM particles (in internal units).
  comoving_baryon_softening:     0.0026994 # Comoving Plummer-equivalent softening length for baryon particles (in internal units).
  max_physical_baryon_softening: 0.0007    # Maximal Plummer-equivalent softening length in physical coordinates for baryon particles (in internal units).
  comoving_nu_softening:         0.0026994 # Comoving Plummer-equivalent softening length for neutrino particles (in internal units).
  max_physical_nu_softening:     0.0007    # Maximal Plummer-equivalent softening length in physical coordinates for neutrino particles (in internal units).
  softening_ratio_background:    0.04      # Fraction of the mean inter-particle separation to use as Plummer-equivalent softening for the background DM particles.

# Parameters governing the time integration (Set dt_min and dt_max to the same value for a fixed time-step run.)
TimeIntegration:
  time_begin:        0.    # The starting time of the simulation (in internal units).
  time_end:          1.    # The end time of the simulation (in internal units).
  dt_min:            1e-6  # The minimal time-step size of the simulation (in internal units).
  dt_max:            1e-2  # The maximal time-step size of the simulation (in internal units).

# Parameters governing the snapshots
Snapshots:
  basename:   output      # Common part of the name of output files.
  subdir:     dir         # (Optional) Sub-directory in which to write the snapshots. Defaults to "" (i.e. the directory where SWIFT is run).
  scale_factor_first: 0.1 # (Optional) Scale-factor of the first snapshot if cosmological time-integration.
  delta_time: 0.01        # Time difference between consecutive outputs (in internal units)
  invoke_stf: 0           # (Optional) Call VELOCIraptor every time a snapshot is written irrespective of the VELOCIraptor output strategy.
  invoke_fof: 0           # (Optional) Call FOF every time a snapshot is written
  compression: 0          # (Optional) Set the level of GZIP compression of the HDF5 datasets [0-9]. 0 does no compression. The lossless compression is applied to *all* the fields.
  distributed: 0          # (Optional) When running over MPI, should each rank write a partial snapshot or do we want a single file? 1 implies one file per MPI rank.
  UnitMass_in_cgs:     1  # (Optional) Unit system for the outputs (Grams)
  UnitLength_in_cgs:   1  # (Optional) Unit system for the outputs (Centimeters)
  UnitVelocity_in_cgs: 1  # (Optional) Unit system for the outputs (Centimeters per second)
  UnitCurrent_in_cgs:  1  # (Optional) Unit system for the outputs (Amperes)
  UnitTemp_in_cgs:     1  # (Optional) Unit system for the outputs (Kelvin)
  output_list_on:      0  # (Optional) Enable the output list
  output_list:         snaplist.txt # (Optional) File containing the output times (see documentation in "Parameter File" section)
  select_output_on:    0  # (Optional) Enable the output selection behaviour
  select_output:       selectoutput.yml # (Optional) File containing information to select outputs with (see documentation in the "Output Selection" section)
  run_on_dump:         0 # (Optional) Run the dump_command each time that a snapshot is dumped?
  dump_command:        ./submit_velociraptor.sh # (Optional) Command to run each time that a snapshot is dumped.

# Parameters governing the conserved quantities statistics
Statistics:
  delta_time:           1e-2        # Time between statistics output

# Parameters related to the initial conditions
InitialConditions:
  file_name:  SedovBlast/sedov.hdf5 # The file to read
  periodic:                    1    # Are we running with periodic ICs?
  generate_gas_in_ics:         0    # (Optional) Generate gas particles from the DM-only ICs (e.g. from panphasia).
" > $1

printf '%s\n' "%  Relevant files

InitCondFile	 $dir_ics/monofonic_${z_in}
OutputDir	 $dir_snap

#EnergyFile         energy.txt
#InfoFile           info.txt
#TimingsFile        timings.txt
#CpuFile            cpu.txt

#RestartFile        restart
SnapshotFileBase   snapshot

OutputListFilename	 /mnt/home/student/cprem/cosmo-sims/gadget4/outputs-many-lin.list

% CPU time -limit

%TimeLimitCPU      36000  % = 10 hours
TimeLimitCPU      691200
#ResubmitOn        0
#ResubmitCommand   my-scriptfile  


% Code options


ICFormat                 3
SnapFormat               3
ComovingIntegrationOn    1

#TypeOfTimestepCriterion  0
OutputListOn 		  $outlston
#PeriodicBoundariesOn     1

%  Caracteristics of run

TimeBegin		 $a_in
TimeMax		 1.0

Omega0		 $Om0
OmegaLambda	 $OmL
OmegaBaryon	 $Omb
HubbleParam	 $h_par
BoxSize		 $boxsize
Hubble 100.0

% Output frequency

%% start at redshift 4
TimeBetSnapshot		1.08379838673
TimeOfFirstSnapshot	0.2


CpuTimeBetRestartFile     300.0    % here in seconds
TimeBetStatistics         0.02

NumFilesPerSnapshot		 8
#NumFilesWrittenInParallel	 8 
MaxFilesWithConcurrentIO	 8


% Accuracy of time integration

ErrTolIntAccuracy      0.05

#MaxRMSDisplacementFac  0.15

CourantFac             0.15     

MaxSizeTimestep       $timestep
MinSizeTimestep       0.0




% Tree algorithm, force accuracy, domain update frequency

ErrTolTheta            0.5
ErrTolThetaMax 	       1.0            
TypeOfOpeningCriterion 1
ErrTolForceAcc         0.005


%Domain Allocation
ActivePartFracForNewDomainDecomp 0.01
TopNodeFactor 2.5

#TreeDomainUpdateFrequency    0.05


%  Further parameters of SPH

DesNumNgb              33
MaxNumNgbDeviation     2
ArtBulkViscConst       0.8
InitGasTemp            1000.0        % always ignored if set to 0 
#MinGasTemp             50.0    
MinEgySpec 		50.0

% Memory allocation

#PartAllocFactor       1.6
#TreeAllocFactor       0.8
#BufferSize            250          % in MByte
MaxMemSize             4000

% System of units

UnitLength_in_cm         3.08568025e24      ;  1.0 Mpc/h 
UnitMass_in_g            1.989e43           ;  1.0e10 Msun/h
UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec 
GravityConstantInternal  0
 

%---- Gravitational softening length
SofteningComovingClass0      $softlen0           ;
SofteningMaxPhysClass0       $softlen0 

SofteningComovingClass1      $softlen1           ;
SofteningMaxPhysClass1       $softlen1 

SofteningComovingClass2      $softlen2         ;
SofteningMaxPhysClass2       $softlen2



SofteningClassOfPartType0    0
SofteningClassOfPartType1    1
SofteningClassOfPartType2    2
SofteningClassOfPartType3    2
SofteningClassOfPartType4    2
SofteningClassOfPartType5    2
% Softening lengths

#MinGasHsmlFractional 0.25" > $1

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