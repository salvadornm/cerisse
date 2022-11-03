#ifndef CNS_INDEX_MACROS_H
#define CNS_INDEX_MACROS_H

#if (AMREX_SPACEDIM == 2)

#define URHO  0
#define UMX   1
#define UMY   2
#define UEDEN 3
#define UEINT 4
#define UTEMP 5
#define UFS   6
#define NVAR  6 + NUM_SPECIES
#define NUM_STATE NVAR // * NUM_FIELDS

#define NCONS 6 + NUM_SPECIES
#define NEQNS 4

#define QRHO   0
#define QU     1
#define QV     2
#define QPRES  3
#define QCS    4
#define QEINT  5
#define QTEMP  6
#define QFS    7
#define NPRIM  7 + NUM_SPECIES
#define NUM_PRIM NPRIM // * NUM_FIELDS

#else

#define URHO  0
#define UMX   1
#define UMY   2
#define UMZ   3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define UFS   7
#define NVAR  7 + NUM_SPECIES
#define NUM_STATE NVAR // * NUM_FIELDS

#define NCONS 7 + NUM_SPECIES // TODO: NEED TO REARRANGE VALUES
#define NEQNS 5

#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QPRES  4
#define QCS    5
#define QEINT  6
#define QTEMP  7
#define QFS    8
#define NPRIM  8 + NUM_SPECIES
#define NUM_PRIM NPRIM // * NUM_FIELDS

#endif

#define CMU    0  // dynamic viscosity
#define CXI    1  // bulk viscosity
#define CLAM   2  // heat conductivity
#define CRHOD  3  // species diffusivity
#define NCOEF  3 + NUM_SPECIES

#define NUM_GROW 5

#endif
