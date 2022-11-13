#ifndef CNS_INDEX_MACROS_H
#define CNS_INDEX_MACROS_H

// Grow cell number
#define NUM_GROW 5

// Number of stochastic fields
#ifndef NUM_FIELD
#define NUM_FIELD 0
#endif

// Aux states
#ifndef NUM_AUX
#define NUM_AUX 0
#endif

// State variables 
// (_MeanField, _Field0, _Field1, ..., Aux_State)
// Index per field
#define URHO   0
#define UMX    1
#define UMY    2
#define UMZ    3 // We have UMY and UMZ no matter what AMREX_SPACEDIM
#define UEDEN  4
#define UEINT  5
#define UTEMP  6
#define UFS    7
#define NVAR   int(7 + NUM_SPECIES)
#define UFA    int(NVAR * (NUM_FIELD + 1))
// Total components of state fab
#define LEN_STATE  int(UFA + NUM_AUX)

// Primitive variables (_MeanField, _Field0, _Field1, ...)
#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QPRES  4
#define QCS    5
#define QEINT  6
#define QTEMP  7
#define QG     8
#define QFS    9
#define NPRIM  int(9 + NUM_SPECIES)
#define LEN_PRIM  int(NPRIM * (NUM_FIELD + 1))

// Transport coefficients (_MeanField, _Field0, _Field1, ...)
#define CMU    0  // dynamic viscosity
#define CXI    1  // bulk viscosity
#define CLAM   2  // heat conductivity
#define CRHOD  3  // species diffusivity
#define NCOEF  int(3 + NUM_SPECIES)
#define LEN_COEF  int(NCOEF * (NUM_FIELD + 1))

// Reaction variables (_MeanField, _Field0, _Field1, ...)
#define NREACT     int(NUM_SPECIES + 2)
#define LEN_REACT  int(NREACT * (NUM_FIELD + 1))

#endif
