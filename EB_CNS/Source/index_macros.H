#ifndef CNS_INDEX_MACROS_H
#define CNS_INDEX_MACROS_H

#include <mechanism.H>

// Do EB?
#ifndef CNS_USE_EB
#define CNS_USE_EB 0
#endif

// Grow cell number
#if CNS_USE_EB
#define NUM_GROW 6 // add 3 for redistribution
#else
#define NUM_GROW 3
#endif

// Number of stochastic fields
#ifndef NUM_FIELD
#define NUM_FIELD 0
#endif

// Aux states
#ifndef NUM_AUX
#define NUM_AUX 0
#endif

// State variables (S)
// (_MeanField, _Field0, _Field1, ..., Aux_State)
// Index per field
constexpr int URHO = 0;
constexpr int UMX = 1;
constexpr int UMY = 2;
constexpr int UMZ = 3;
constexpr int UEDEN = 4;
constexpr int UTEMP = 5;
constexpr int UFS = 6;
constexpr int NVAR = 6 + NUM_SPECIES;       // Number of variables per field
constexpr int UFA = NVAR * (NUM_FIELD + 1); // Aux variables
constexpr int LEN_STATE = UFA + NUM_AUX;    // Total components of state fab

// Primitive variables (Q)
constexpr int QRHO = 0;
constexpr int QU = 1;
constexpr int QV = 2;
constexpr int QW = 3;
constexpr int QPRES = 4;
constexpr int QC = 5; // speed of sound
constexpr int QG = 6; // gamma
constexpr int QEINT = 7;
constexpr int QTEMP = 8;
constexpr int QWBAR = 9; // average molecular weight
constexpr int QFS = 10;
constexpr int NPRIM = QFS + NUM_SPECIES;
constexpr int LEN_PRIM = NPRIM * (NUM_FIELD + 1);

// Characteristic variables (W)
constexpr int WRHO = 0;                // material/entropy wave (rho)
constexpr int WY = 1;                  // Y, simply advected
constexpr int WACO = NUM_SPECIES + 1;  // acoustic waves
constexpr int WTHER = NUM_SPECIES + 3; // thermodynamic variable (sos or gamma)
constexpr int WUT = NUM_SPECIES + 4;   // tangential velocities
constexpr int NCHAR = NUM_SPECIES + AMREX_SPACEDIM + 3;

// Transport coefficients (_MeanField, _Field0, _Field1, ...)
constexpr int CMU = 0;   // dynamic viscosity
constexpr int CXI = 1;   // bulk viscosity
constexpr int CLAM = 2;  // heat conductivity
constexpr int CRHOD = 3; // species diffusivity
constexpr int NCOEF = 3 + NUM_SPECIES;
constexpr int LEN_COEF = NCOEF * (NUM_FIELD + 1);

// Reaction variables (_MeanField, _Field0, _Field1, ...)
constexpr int NREACT = NUM_SPECIES + 1;
constexpr int LEN_REACT = NREACT * (NUM_FIELD + 1);

#endif
