#!/bin/bash

echo $2
eval "$2"

dir_sim_data=/scratch/cprem/sims/$simnm/$rund/
dir_snap=$dir_sim_data/snaps/

mkdir -p ${dir_sim_data}/halos_rs/
# mkdir -p $1

printf '%s\n' "INBASE = $dir_snap
FORCE_RES = 0.0065
NUM_SNAPS = 51
STARTING_SNAP = 0

NUM_WRITERS = 32
FORK_PROCESSORS_PER_MACHINE = 32

OUTBASE = ${dir_sim_data}/halos_rs/
FILENAME = snapdir_<snap>/snapshot_<snap>.<block>.hdf5
NUM_BLOCKS = 8

FILE_FORMAT = AREPO
AREPO_LENGTH_CONVERSION = 1
AREPO_MASS_CONVERSION = 1e+10

PARALLEL_IO = 1
FORK_READERS_FROM_WRITERS = 1
PERIODIC = 1

TEMPORAL_HALO_FINDING = 1
STRICT_SO_MASSES = 0

MIN_HALO_OUTPUT_SIZE = 30
MIN_HALO_PARTICLES = 20
UNBOUND_THRESHOLD = 0.5

OVERLAP_LENGTH = 3" > $1