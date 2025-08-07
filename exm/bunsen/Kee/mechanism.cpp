#include "mechanism.H"
const int rmap[57] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56};

// Returns 0-based map of reaction order
void GET_RMAP
(int * _rmap)
{
for (int j=0; j<57; ++j)
{
_rmap[j] = rmap[j];
}
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int& nspec, int ki[], int nu[])
{
const int ns[57] =
     {3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,4,3,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,3,4,3,2,2,2,3,3,4,3,2,4,4};
const int kiv[285] =
     {1,9,0,0,0,0,10,1,13,0,0,9,1,8,0,0,11,1,12,0,0,12,1,15,0,1,11,4,9,0,1,12,4,8,0,1,12,2,15,0,1,9,2,8,0,2,9,3,8,0,2,12,4,9,0,2,12,3,15,0,3,10,5,11,0,3,11,7,9,0,3,12,9,5,0,3,6,7,5,0,2,6,4,7,0,2,11,7,9,0,2,11,7,8,0,2,10,6,9,0,2,10,4,11,0,2,10,6,8,0,2,10,7,15,0,2,10,7,9,12,2,10,5,12,0,4,12,15,5,0,4,9,8,5,0,4,9,5,0,0,4,11,5,12,0,5,12,7,15,0,5,7,9,0,0,9,5,7,8,0,5,11,6,9,0,5,10,7,13,0,7,11,6,0,0,7,12,6,9,0,7,10,6,11,0,7,13,6,12,0,8,10,12,0,0,8,12,9,15,0,9,10,11,12,0,8,11,9,12,0,9,10,13,0,0,13,12,15,10,0,9,13,12,0,0,13,11,10,12,0,12,15,11,0,0,9,8,0,0,0,9,8,0,0,0,9,8,0,0,0,9,12,15,0,0,9,11,12,0,0,9,13,8,10,0,13,14,10,0,0,14,12,0,0,0,9,14,8,13,0,14,12,15,13,0};
const int nuv[285] =
     {-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-2,1,1,0,0,-2,1,0,0,0,-2,1,0,0,0,-2,1,0,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-2,1,1,0,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0};
if (i < 1) {
// Return max num species per reaction
nspec = 5;
} else {
if (i > 57) {
nspec = -1;
} else {
nspec = ns[i-1];
for (int j=0; j<nspec; ++j) {
ki[j] = kiv[(i-1)*5 + j] + 1;
nu[j] = nuv[(i-1)*5 + j];
}
}
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real x[], amrex::Real q_f[], amrex::Real q_r[])
{
amrex::Real c[17]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 17; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, T);

// convert to chemkin units
for (int id = 0; id < 57; ++id) {
q_f[id] *= 1.0e-6;
q_r[id] *= 1.0e-6;
}
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };// temperature cache
amrex::Real invT = 1.0 / tc[1];
// compute the Gibbs free energy
amrex::Real g_RT[17];
gibbs(g_RT, tc);

amrex::Real sc_qss[1];
comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 1.008000; // H
awt[1] = 15.999000; // O
awt[2] = 12.011000; // C
awt[3] = 14.007000; // N
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int kd = 4; 
// Zero ncf
for (int id = 0; id < kd * 17; ++ id) {
 ncf[id] = 0; 
}

// CH4
ncf[ 0 * kd + 2 ] = 1; // C
ncf[ 0 * kd + 0 ] = 4; // H

// CH3
ncf[ 1 * kd + 2 ] = 1; // C
ncf[ 1 * kd + 0 ] = 3; // H

// CH2
ncf[ 2 * kd + 2 ] = 1; // C
ncf[ 2 * kd + 0 ] = 2; // H

// CH
ncf[ 3 * kd + 2 ] = 1; // C
ncf[ 3 * kd + 0 ] = 1; // H

// CH2O
ncf[ 4 * kd + 2 ] = 1; // C
ncf[ 4 * kd + 0 ] = 2; // H
ncf[ 4 * kd + 1 ] = 1; // O

// HCO
ncf[ 5 * kd + 2 ] = 1; // C
ncf[ 5 * kd + 0 ] = 1; // H
ncf[ 5 * kd + 1 ] = 1; // O

// CO2
ncf[ 6 * kd + 2 ] = 1; // C
ncf[ 6 * kd + 1 ] = 2; // O

// CO
ncf[ 7 * kd + 2 ] = 1; // C
ncf[ 7 * kd + 1 ] = 1; // O

// H2
ncf[ 8 * kd + 0 ] = 2; // H

// H
ncf[ 9 * kd + 0 ] = 1; // H

// O2
ncf[ 10 * kd + 1 ] = 2; // O

// O
ncf[ 11 * kd + 1 ] = 1; // O

// OH
ncf[ 12 * kd + 0 ] = 1; // H
ncf[ 12 * kd + 1 ] = 1; // O

// HO2
ncf[ 13 * kd + 0 ] = 1; // H
ncf[ 13 * kd + 1 ] = 2; // O

// H2O2
ncf[ 14 * kd + 0 ] = 2; // H
ncf[ 14 * kd + 1 ] = 2; // O

// H2O
ncf[ 15 * kd + 0 ] = 2; // H
ncf[ 15 * kd + 1 ] = 1; // O

// N2
ncf[ 16 * kd + 3 ] = 2; // N

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(4);
ename[0] = "H";
ename[1] = "O";
ename[2] = "C";
ename[3] = "N";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(17);
kname[0] = "CH4";
kname[1] = "CH3";
kname[2] = "CH2";
kname[3] = "CH";
kname[4] = "CH2O";
kname[5] = "HCO";
kname[6] = "CO2";
kname[7] = "CO";
kname[8] = "H2";
kname[9] = "H";
kname[10] = "O2";
kname[11] = "O";
kname[12] = "OH";
kname[13] = "HO2";
kname[14] = "H2O2";
kname[15] = "H2O";
kname[16] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<18; k++) {
for (int l=0; l<18; l++) {
if(Jac[ 18 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<18; k++) {
for (int l=0; l<18; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 18 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the simplified (for preconditioning) system Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<18; k++) {
for (int l=0; l<18; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 18 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

nJdata[0] = nJdata_tmp;
}


// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 18;
int offset_col = nc * 18;
for (int k=0; k<18; k++) {
for (int l=0; l<18; l++) {
if(Jac[18*k + l] != 0.0) {
rowVals[nJdata_tmp] = l + offset_row; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
colPtrs[offset_col + (k + 1)] = nJdata_tmp;
}
}
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 18;
for (int l=0; l<18; l++) {
for (int k=0; k<18; k++) {
if(Jac[18*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 18;
for (int l=0; l<18; l++) {
for (int k=0; k<18; k++) {
if(Jac[18*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 18;
for (int l=0; l<18; l++) {
for (int k=0; k<18; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[18*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 18;
for (int l=0; l<18; l++) {
for (int k=0; k<18; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[18*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU
// BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<18; k++) {
for (int l=0; l<18; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 18*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[18*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 18*k + l;
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
colPtrs[k+1] = nJdata_tmp;
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
amrex::GpuArray<amrex::Real,324> Jac = {0.0};
amrex::GpuArray<amrex::Real,17> conc = {0.0};
for (int n=0; n<17; n++) {
    conc[n] = 1.0/ 17.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<18; l++) {
for (int k=0; k<18; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[18*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int l=0; l<18; l++) {
for (int k=0; k<18; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[18*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
