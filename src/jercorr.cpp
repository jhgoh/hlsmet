/*
HLS implementation of MET correction with Jet energy resolution
*/
#include "common.h"
#include "jercorr.h"

#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

// pt, phi are integers
void jercorrmet_hw(pt_t data_met,
                   pt_t in_pt[NJET], phi_t in_phi[NJET], phi_t in_eta[NJET],
                   pt_t& out_met)
{
  #pragma HLS ARRAY_PARTITION variable=in_pt complete
  #pragma HLS ARRAY_PARTITION variable=in_phi complete
  #pragma HLS ARRAY_PARTITION variable=in_eta complete
  #pragma HLS ARRAY_PARTITION variable=dpt complete

  #pragma HLS pipeline ii=54

  pt_t kfactor = 1.2; // FIXME: check the k-factor.

  // Obtain the Jet energy resolution
  // FIXME: have to use proper LUT implementation, maybe.
  //        this first implementation is based on if-else statements
  //        inside of the loop which maybe inefficient...
  pt_t dpt[NJET];

  phi_t var1s[] = {0.016, 0.216, 1, 0.083, 1}; // NOTE: used phi_t intentionally, these coeffs are small enough
  pt_t var2s[] = {6.6, 5.5, 1, 13.2, 1}; 
  phi_t etaedges[] = {1.3, 1.7, 2.5, 3.0};

  LOOP_JER:
  for ( int i=0; i<NJET; ++i ) {
    #pragma HLS UNROLL
    phi_t abseta = hls::abs(in_eta[i]);

    int etabin = 0;
    if ( abseta < etaedges[0] ) etabin = 0;
    else if ( abseta < etaedges[1] ) etabin = 1;
    else if ( abseta < etaedges[2] ) etabin = 2;
    else if ( abseta < etaedges[3] ) etabin = 3;
    else etabin = 4;

    // The MET significance correction formula is obtained with dPt/Pt
    dpt[i] = (var1s[etabin]*in_pt[i] + var2s[etabin])*in_pt[i];
  }

  // Get the sum of dx^2 and dy^2, 
  // treating JER variations of each jet to be independent.
  pt2_t dpt2 = 0;
  SUM_dpt2:
  for(int i=0; i<NJET; ++i) {
    #pragma HLS unroll
    dpt2 += dpt[i]*dpt[i];
  }

  // Apply the MET significance 'correction', being permissive to trigger the event
  out_met = data_met + hls::sqrt(dpt2);

  return;
}

