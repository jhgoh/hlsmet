/*
HLS implementation of MET calculation from PF objects
*/
#include "corr.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#include "hls_math.h"
#endif

// pt, phi are integers
// **** The Problem is started from here. Because the testbench code(corr_test.cpp) gives the values which was shifted.
void corr_hw(pt_t data_pt, phi_t data_phi, 
		pt_t jet_pt[NPART], phi_t jet_phi[NPART], phi_t jet_eta[NPART], 
		 pt2_t& corr_pt2){

    #pragma HLS ARRAY_PARTITION variable=jet_pt complete
    #pragma HLS ARRAY_PARTITION variable=jet_phi complete
    #pragma HLS ARRAY_PARTITION variable=jet_eta complete

    #pragma HLS pipeline ii=15
    

    // calc signed components first
    pxy_t met_x; 
    pxy_t met_y;
    // Get x, y components of pippiMET
    ProjX(data_pt, data_phi, met_x); // **** met_x {2bits shifted from TB}
    ProjY(data_pt, data_phi, met_y); // **** met_y {2bits shifted from TB}

    pxy_t jet_x[NPART];
    pxy_t jet_y[NPART];


    LOOP_PROJECT: for(int i=0; i<NPART;i++){
	#pragma HLS UNROLL
        // Get x, y Jet components
        ProjX(jet_pt[i], jet_phi[i], jet_x[i]); // **** jet_x {2bits shifted from TB}
        ProjY(jet_pt[i], jet_phi[i], jet_y[i]); // **** jet_j {2bits shifted from TB}
		jet_eta[i] = hls::absf(jet_eta[i]);
		
    }

    pt2_t sumJ_x = 0;
    pt2_t sumJ_y = 0;

	phi_t eta13 = 1.3;
	phi_t eta17 = 1.7;
	phi_t eta25 = 2.5;
	phi_t eta30 = 3.0;


    SUM: for(int i=0; i<NPART;i++){
	#pragma HLS UNROLL
		pt_t var1 = 0;
		pt_t var2 = 0;

		if ( jet_eta[i] < eta13){
			var1 = 0.106;
			var2 = 6.6;
		}

		if ( (jet_eta[i] < eta17 && jet_eta[i] > eta13 )){ 
			var1 = 0.216;
			var2 = 5.5;
		}

		if ( jet_eta[i] < eta30 && jet_eta[i] > eta25 ){
			var1 = 0.083;
			var2 = 13.2;
		}

		sumJ_x += (var1*(jet_x[i])+var2) * (var1*(jet_x[i])+var2);
		sumJ_y += (var1*(jet_y[i])+var2) * (var1*(jet_y[i])+var2);
    }


	// correctd MET
	pt_t K_value = 1.5; // **** K_value {2 bits shifted}
	pt2_t corr_x;
	pt2_t corr_y;

    corr_x = met_x + (K_value * pt_t(hls::sqrtf(sumJ_x))); // **** corr_x {2 bits shifted}
	corr_y = met_y + (K_value * pt_t(hls::sqrtf(sumJ_y))); // **** corr_y {2 bits shifted}

	corr_pt2 = corr_x * corr_x + corr_y * corr_y;  // **** corr_pt2 {4 bits shifted}


    return;
}

