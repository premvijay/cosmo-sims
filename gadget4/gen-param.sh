#!/bin/bash
softlen1=0.0065
timestep=0.01
# z_in=24
# Om0=0.3063375
# OmL=0.6936625
# Omb=0.0484103
# h_par=0.6781
# boxsize=200.0
echo $2
eval "$2"

dir_root=$HOME/cosmo-sims/gadget4/
dir_sim=$dir_root/$simnm/
dir_sim_data=/scratch/cprem/sims/$simnm/$rund/
dir_snap=$dir_sim_data/snaps/
dir_ics=$dir_sim_data/ics/

# a_in=$(printf %.5f `echo "1/($z_in+1)" | bc -l`)
a_in=$(python -c "print (1.0/(1+$z_in))")

printf '%s\n' "%  Relevant files

InitCondFile	 $dir_ics/monofonic_${z_in}
OutputDir	 $dir_snap

#EnergyFile         energy.txt
#InfoFile           info.txt
#TimingsFile        timings.txt
#CpuFile            cpu.txt

#RestartFile        restart
SnapshotFileBase   snapshot

OutputListFilename	 dummy.txt

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
OutputListOn 		  0
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


CpuTimeBetRestartFile     1800.0    % here in seconds
TimeBetStatistics         0.05

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
MaxMemSize             9000

% System of units

UnitLength_in_cm         3.08568025e24      ;  1.0 Mpc/h 
UnitMass_in_g            1.989e43           ;  1.0e10 Msun/h
UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec 
GravityConstantInternal  0
 

%---- Gravitational softening length
SofteningComovingClass0      0.5           ; 
SofteningMaxPhysClass0       0.5

SofteningComovingClass1      $softlen1           ;
SofteningMaxPhysClass1       $softlen1 

SofteningComovingClass2      0.005         ;
SofteningMaxPhysClass2       0.005



SofteningClassOfPartType0    0
SofteningClassOfPartType1    1
SofteningClassOfPartType2    2
SofteningClassOfPartType3    2
SofteningClassOfPartType4    2
SofteningClassOfPartType5    2
% Softening lengths

#MinGasHsmlFractional 0.25" > $1

