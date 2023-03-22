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
void corr_hw(pt_t data_pt[NPART], phi_t data_phi[NPART], phi_t data_eta[NPART], 
		pt_t jet_pt[NPART], phi_t jet_phi[NPART], 
		pt2_t& res_pt2, phi_t& res_phi, pt2_t& corr_pt2){
    #pragma HLS ARRAY_PARTITION variable=data_pt complete
    #pragma HLS ARRAY_PARTITION variable=data_phi complete
    #pragma HLS ARRAY_PARTITION variable=data_eta complete

    #pragma HLS pipeline ii=54
    
    if(DEBUG) std::cout << "  HW Begin" << std::endl;

    pt_t data_pt_int[NPART];
    phi_t data_phi_int[NPART];
    phi_t data_eta_int[NPART];
    #pragma HLS ARRAY_PARTITION variable=data_pt_int complete
    #pragma HLS ARRAY_PARTITION variable=data_phi_int complete
    #pragma HLS ARRAY_PARTITION variable=data_eta_int complete

    pt_t data_J_pt_int[NPART];
    phi_t data_J_phi_int[NPART];
    #pragma HLS ARRAY_PARTITION variable=data_J_pt_int complete
    #pragma HLS ARRAY_PARTITION variable=data_J_phi_int complete

    FILL_INTERNAL: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        data_pt_int[i]=data_pt[i];
        data_phi_int[i]=data_phi[i];
        data_eta_int[i]=data_eta[i];

        data_J_pt_int[i]=jet_pt[i];
        data_J_phi_int[i]=jet_phi[i];
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


    pxy_t sum_x = 0;
    pxy_t sum_y = 0;
    pxy_t sumJ_x = 0;
    pxy_t sumJ_y = 0;
    SUM: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        // Add to x, y sums
        sum_x -= met_x[i];
        sum_y -= met_y[i];
        sumJ_x += jet_x[i];
        sumJ_y += jet_y[i];
         if(!DEBUG){
             std::cout << "     met x,y = (" << -met_x[i] << ", " << -met_y[i] << ") \t";
             std::cout << " sum x,y = (" << sum_x << ", " << sum_y << ") \t";
             std::cout << " eta  = (" << data_eta[i] << ", ) \t";
             std::cout << " phi  = (" << data_phi[i] << ", ) \t";
             std::cout << " \n";
         }
    }
	// mean of jet pt	
	pxy_t mean_x = sumJ_x / NPART;
	pxy_t mean_y = sumJ_y / NPART;

	// SUM of (x - mean)*(x -mean)
	pt2_t dev_x = 0;
	pt2_t dev_y = 0;
	for(int i=0; i<NPART; i++){
		dev_x += ( (jet_x[i] - mean_x) * (jet_x[i] - mean_x) );
		dev_y += ( (jet_y[i] - mean_y) * (jet_y[i] - mean_y) );
	}

	// stdev of jet
	pt2_t sigma2_x = dev_x / (NPART-1);
	pt2_t sigma2_y = dev_y / (NPART-1);

	// sqrt(sigma2)
	pt_t sigma_x = hls::sqrt(sigma2_x); 
	pt_t sigma_y = hls::sqrt(sigma2_y);

	// MET * MET
    res_pt2 = sum_x*sum_x + sum_y*sum_y;

	// sqrt(MET * MET)
	pt_t res_pt = hls::sqrt(res_pt2);

	// correctd MET
	// K factor?
	pt_t K_value = 1.;
	pt2_t corr_x;
	pt2_t corr_y;

	corr_x = sum_x + (K_value * sigma_x);
	corr_y = sum_y + (K_value * sigma_y);

	corr_pt2 = corr_x*corr_x + corr_y*corr_y; 
    //PhiFromXY(sum_x,sum_y,res_phi);
    res_phi=0;

	std::cout << " corr_pt2 = " << corr_pt2 << " \n";


    return;
}

