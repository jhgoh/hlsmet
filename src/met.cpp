/*
HLS implementation of MET calculation from PF objects
*/
#include "common.h"
#include "met.h"

#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

// pt, phi are integers
void met_hw(pt_t in_pt[NPART], phi_t in_phi[NPART], pt2_t& out_pt2, phi_t& out_phi)
{
#pragma HLS ARRAY_PARTITION variable=in_pt complete
#pragma HLS ARRAY_PARTITION variable=in_phi complete

#pragma HLS pipeline ii=54
      
#if(DEBUG==1)
  std::cout << "  HW Begin" << std::endl;
#endif

  pxy_t met_x = 0;
  pxy_t met_y = 0;
  LOOP_PROJECT: for(int i=0; i<NPART; ++i) {
#ifdef INTONLY
    met_x -= in_pt[i] * f2hwPhi(cos(hw2fPhi(in_phi[i])));
    met_y -= in_pt[i] * f2hwPhi(sin(hw2fPhi(in_phi[i])));
#if(DEBUG==1)
std::cout << "DEBUG(HW/INT) " << i << ' ' << hw2fPt(in_pt[i]) << ' ' << hw2fPhi(in_phi[i]) << std::endl;
#endif
#else
#pragma HLS UNROLL
    met_x -= in_pt[i] * hls::cos(in_phi[i]);
    met_y -= in_pt[i] * hls::sin(in_phi[i]);
#if(DEBUG==1)
std::cout << "DEBUG(HW/HLS): " << i << ' ' << in_pt[i] << ' ' << in_phi[i] << ' ' << hls::cos(in_phi[i]) << ' ' << hls::sin(in_phi[i]) << std::endl;
#endif
#endif
  }

  out_pt2 = met_x*met_x + met_y*met_y;
  out_phi = hls::atan2(met_y, met_x);
#if(DEBUG==1)
  std::cout << "DEBUG(HW) MET pT2=" << out_pt2 << " Phi=" << out_phi << " px=" << met_x << " py=" << met_y << std::endl;
#endif

  return;
}

