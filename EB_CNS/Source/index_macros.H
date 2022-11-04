#ifndef CNS_INDEX_MACROS_H
#define CNS_INDEX_MACROS_H

// State variables 
// (_Field0, _Field1, ..., (if NUM_FIELD > 1, Mean_State), Aux_State)
// Index per field
  #define URHO   0
  #define UMX    1
  #define UMY    2
  #define UMZ    3 // We have UMZ no matter what AMREX_SPACEDIM
  #define UEDEN  4
  #define UEINT  5
  #define UTEMP  6
  #define UFS    7
  #define NVAR   int(7 + NUM_SPECIES)

// Aux states
#ifndef NUM_AUX
  #define NUM_AUX  0
#endif

// For NUM_FIELD fields
#ifndef NUM_FIELD
  #define NUM_FIELD  1
#endif
#if (NUM_FIELD > 1)
  #define UMF  NVAR * NUM_FIELD // Mean field starting entry
#else
  #define UMF  URHO // Mean field starting entry
  
#endif
  #define UFA        int(UMF + NVAR)
  #define LEN_STATE  int(UFA + NUM_AUX)

// Primitive variables (_Field0, _Field1, ...)
  #define QRHO   0
  #define QU     1
  #define QV     2
  #define QW     3
  #define QPRES  4
  #define QCS    5
  #define QEINT  6
  #define QTEMP  7
  #define QFS    8
  #define NPRIM  int(8 + NUM_SPECIES)
  #define LEN_PRIM  int(NPRIM * NUM_FIELD)

// Transport coefficients (_Field0, _Field1, ...)
  #define CMU    0  // dynamic viscosity
  #define CXI    1  // bulk viscosity
  #define CLAM   2  // heat conductivity
  #define CRHOD  3  // species diffusivity
  #define NCOEF  int(3 + NUM_SPECIES)
  #define LEN_COEF  int(NCOEF * NUM_FIELD)

// Reaction variables
// (_Field0, _Field1, ..., (if NUM_FIELD > 1, Mean_State))
  #define NREACT     int(NUM_SPECIES + 2)
#if (NUM_FIELD > 1)
  #define LEN_REACT  int(NREACT * (NUM_FIELD + 1))
#else
  #define LEN_REACT  NREACT
#endif

  #define NUM_GROW 5

#endif
