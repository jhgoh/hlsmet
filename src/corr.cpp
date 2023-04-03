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
void corr_hw(pt_t data_pt, phi_t data_phi, 
		pt_t jet_pt[NPART], phi_t jet_phi[NPART], phi_t jet_eta[NPART], 
		 pt2_t& corr_pt2){
//    #pragma HLS ARRAY_PARTITION variable=data_pt complete
//    #pragma HLS ARRAY_PARTITION variable=data_phi complete
    #pragma HLS ARRAY_PARTITION variable=jet_pt complete
    #pragma HLS ARRAY_PARTITION variable=jet_phi complete
    #pragma HLS ARRAY_PARTITION variable=jet_eta complete

    #pragma HLS pipeline ii=54
    
    if(DEBUG) std::cout << "  HW Begin" << std::endl;

    pt_t data_pt_int;
    phi_t data_phi_int;
//    #pragma HLS ARRAY_PARTITION variable=data_pt_int complete
 //   #pragma HLS ARRAY_PARTITION variable=data_phi_int complete

    pt_t data_J_pt_int[NPART];
    phi_t data_J_phi_int[NPART];
    phi_t data_J_eta_int[NPART];
    #pragma HLS ARRAY_PARTITION variable=data_J_pt_int complete
    #pragma HLS ARRAY_PARTITION variable=data_J_phi_int complete
    #pragma HLS ARRAY_PARTITION variable=data_J_eta_int complete

    data_pt_int=data_pt;
    data_phi_int=data_phi;

    FILL_INTERNAL: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL

        data_J_pt_int[i]=jet_pt[i];
        data_J_phi_int[i]=jet_phi[i];
        data_J_eta_int[i]=jet_eta[i];
    }

    // calc signed components first
    pxy_t met_x;
    pxy_t met_y;
    pxy_t jet_x[NPART];
    pxy_t jet_y[NPART];
    LOOP_PROJECT: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        // Get x, y components of pippiMET
        ProjX(data_pt_int, data_phi_int, met_x);
        ProjY(data_pt_int, data_phi_int, met_y);
        // Get x, y Jet components
        ProjX(data_J_pt_int[i], data_J_phi_int[i], jet_x[i]);
        ProjY(data_J_pt_int[i], data_J_phi_int[i], jet_y[i]);
    }


    pxy_t sumJ_x = 0;
    pxy_t sumJ_y = 0;

	phi_t eta13 = int(1.3 * (1<<PHI_SIZE));
	phi_t eta17 = int(1.7 * (1<<PHI_SIZE));
	phi_t eta25 = int(2.5 * (1<<PHI_SIZE));
	phi_t eta30 = int(3.0 * (1<<PHI_SIZE));

    SUM: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
		// Minus eta range is not done.
		if ( data_J_eta_int[i] < int(1.3 * (1<<PHI_SIZE)) && data_J_eta_int[i] > int(-1.3 * (1<<PHI_SIZE)) ){
			pxy_t var1 = int(0.106 * (1<<PT_DEC_BITS)); // float to fixed type
			pxy_t var2 = int(6.6 * (1<<PT_DEC_BITS));   // same bit with J_pt ?
        	sumJ_x += (var1*jet_x[i]+var2);
        	sumJ_y += (var1*jet_y[i]+var2);
		}
		if ( data_J_eta_int[i] < int(1.7 * (1<<PHI_SIZE)) && 
				data_J_eta_int[i] > int(1.3 * (1<<PHI_SIZE)) ){
			pxy_t var1 = int(0.216 * (1<<PT_DEC_BITS));
			pxy_t var2 = int(5.5 * (1<<PT_DEC_BITS));
        	sumJ_x += (var1*jet_x[i]+var2);
        	sumJ_y += (var1*jet_y[i]+var2);
		}
		if ( data_J_eta_int[i] < int(3. * (1<<PHI_SIZE)) && 
				data_J_eta_int[i] > int(2.5 * (1<<PHI_SIZE)) ){
			pxy_t var1 = int(0.083 * (1<<PT_DEC_BITS));
			pxy_t var2 = int(13.2 * (1<<PT_DEC_BITS));
        	sumJ_x += (var1*jet_x[i]+var2);
        	sumJ_y += (var1*jet_y[i]+var2);
		}
         if(!DEBUG){
             std::cout << " sigma(Jet) x,y = (" << sumJ_x << ", " << sumJ_y << ") \t";
             std::cout << " \n";
         }
    }


	// correctd MET
	// K factor, HW unit convert?
	pt_t K_value = int(1.5 * (1<<PT_DEC_BITS)); //
	pt2_t corr_x;
	pt2_t corr_y;

	// Now choose one MET because input is N array.
	// Will Fixed the input as one number ?
	corr_x = met_x + (K_value * sumJ_x);
	corr_y = met_y + (K_value * sumJ_y);

	corr_pt2 = corr_x*corr_x + corr_y*corr_y; 
    //PhiFromXY(sum_x,sum_y,res_phi);

	std::cout << " corr_pt2 = " << corr_pt2 << " \n";


    return;
}

