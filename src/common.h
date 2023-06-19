#ifndef HLSMET_COMMON_H
#define HLSMET_COMMON_H

// Configurations for the testbench program
//#define INTONLY
#define NTEST 1000
#define NPART 128
#define NJET 10
#define DEBUG 0

#include "hls_math.h"
#include "ap_int.h"
#include "ap_fixed.h"

// Define datatypes for pT
//  pT is uint where 1 bit = 1/4 GeV, up to 4096 (16 bits)
//    px, py need to be signed
//    pT^2 needs double precision as pT
// TDR: 16 bits, 1/4 GeV. For us, 12 bits probably OK (1024 GeV) but would need to add overflow checks
#define PT_SIZE 14
// bits used to represent the decimal: 2->1/2^2 GeV precision
#define PT_DEC_BITS 2
#ifdef INTONLY
typedef ap_uint<PT_SIZE> pt_t;
typedef ap_int<PT_SIZE+1> pxy_t;
typedef ap_uint<PT_SIZE*2> pt2_t;
inline pt_t f2hwPt(const float pt) { return int(pt*(1<<PT_DEC_BITS)); }
inline float hw2fPt(const pt_t pt) { return float(pt)/(1<<PT_DEC_BITS); }
#else
typedef ap_ufixed<PT_SIZE, PT_SIZE-PT_DEC_BITS> pt_t;
typedef ap_fixed<PT_SIZE+2, PT_SIZE+2-PT_DEC_BITS> pxy_t;
typedef ap_ufixed<(PT_SIZE+2)*2, (PT_SIZE+2)*2-PT_DEC_BITS> pt2_t;
#endif

// Define datatypes for phi
// phi size = 10bits in TDR. For reference, 2pi/(2^10)=0.0006
#define PHI_SIZE 10

#ifdef INTONLY
typedef ap_int<PHI_SIZE> phi_t;
inline phi_t f2hwPhi(const float phi) { return int(phi*(1<<PHI_SIZE)/(2*M_PI)); }
inline float hw2fPhi(const phi_t phi) { return float(phi)*(2*M_PI)/(1<<PHI_SIZE); }
#else
typedef ap_fixed<PHI_SIZE, 3> phi_t;
#endif

// or use manual float<->int conversion

#endif
