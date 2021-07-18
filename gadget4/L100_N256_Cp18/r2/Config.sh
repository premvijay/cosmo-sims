PERIODIC                                     # enables periodic boundary condistions
NTYPES=6                                     # number of particle types
GADGET2_HEADER                               # allows reading of snapshots with GADGET-2/3 header format
SELFGRAVITY                                   # switch to enable self-gravity of particles (typically always on)
PMGRID=512                                   # basic mesh size for TreePM calculations
NSOFTCLASSES=3                               # number of different softening classes
DOUBLEPRECISION=1                             # if activated and set to 1, use double precision internally, for 2 use mixed precision, otherwise single precision
DOUBLEPRECISION_FFTW                         # if set, carries out FFTs in double precision, otherwise single precision
POWERSPEC_ON_OUTPUT                          # computes a matter power spectrum when the code writes a snapshot output
NUMBER_OF_MPI_LISTENERS_PER_NODE=1           # set such that the number of MPI-ranks per node and listener is maller than MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY
MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY=32    # default is 64, but can also be set to 32
