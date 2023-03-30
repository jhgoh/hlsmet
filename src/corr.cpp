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
void corr_hw(pt_t data_pt[NPART], phi_t data_phi[NPART], 
		pt_t jet_pt[NPART], phi_t jet_phi[NPART], phi_t jet_eta[NPART], 
		 pt2_t& corr_pt2){
    #pragma HLS ARRAY_PARTITION variable=data_pt complete
    #pragma HLS ARRAY_PARTITION variable=data_phi complete
    #pragma HLS ARRAY_PARTITION variable=jet_pt complete
    #pragma HLS ARRAY_PARTITION variable=jet_phi complete
    #pragma HLS ARRAY_PARTITION variable=jet_eta complete

    #pragma HLS pipeline ii=54
    
    if(DEBUG) std::cout << "  HW Begin" << std::endl;

    pt_t data_pt_int[NPART];
    phi_t data_phi_int[NPART];
    #pragma HLS ARRAY_PARTITION variable=data_pt_int complete
    #pragma HLS ARRAY_PARTITION variable=data_phi_int complete

    pt_t data_J_pt_int[NPART];
    phi_t data_J_phi_int[NPART];
    phi_t data_J_eta_int[NPART];
    #pragma HLS ARRAY_PARTITION variable=data_J_pt_int complete
    #pragma HLS ARRAY_PARTITION variable=data_J_phi_int complete
    #pragma HLS ARRAY_PARTITION variable=data_J_eta_int complete

    FILL_INTERNAL: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        data_pt_int[i]=data_pt[i];
        data_phi_int[i]=data_phi[i];

        data_J_pt_int[i]=jet_pt[i];
        data_J_phi_int[i]=jet_phi[i];
        data_J_eta_int[i]=jet_eta[i];
    }

    // calc signed components first
    pxy_t met_x[NPART];
    pxy_t met_y[NPART];
    pxy_t jet_x[NPART];
    pxy_t jet_y[NPART];
    LOOP_PROJECT: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        // Get x, y components
        ProjX(data_pt_int[i], data_phi_int[i], met_x[i]);
        ProjY(data_pt_int[i], data_phi_int[i], met_y[i]);
        // Get x, y Jet components
        ProjX(data_J_pt_int[i], data_J_phi_int[i], jet_x[i]);
        ProjY(data_J_pt_int[i], data_J_phi_int[i], jet_y[i]);
    }


    pxy_t sumJ_x = 0;
    pxy_t sumJ_y = 0;
    SUM: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
		if ( data_J_eta_int[i] < int(1.3 * (1<<PHI_SIZE)) || 
				data_J_eta_int[i] > -int(1.3 * (1<<PHI_SIZE)) ){
        	sumJ_x += (0.106*jet_x[i]+6.6);
        	sumJ_y += (0.106*jet_y[i]+6.6);
		}
         if(!DEBUG){
             std::cout << " sum x,y = (" << sumJ_x << ", " << sumJ_y << ") \t";
             std::cout << " eta  = (" << data_eta[i] << ", ) \t";
             std::cout << " phi  = (" << data_phi[i] << ", ) \t";
             std::cout << " \n";
         }
    }



	// correctd MET
	// K factor?
	pt_t K_value = 1.5;
	pt2_t corr_x;
	pt2_t corr_y;

	corr_x = met_x + (K_value * sumJ_x);
	corr_y = met_y + (K_value * sumJ_y);

	corr_pt2 = corr_x*corr_x + corr_y*corr_y; 
    //PhiFromXY(sum_x,sum_y,res_phi);
    res_phi=0;

	std::cout << " corr_pt2 = " << corr_pt2 << " \n";


    return;
}

