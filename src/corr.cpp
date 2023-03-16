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
void corr_hw(pt_t data_pt[NPART], phi_t data_phi[NPART], phi_t data_eta[NPART], pt2_t& res_pt2, phi_t& res_phi, pt2_t& corr_pt){
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

    FILL_INTERNAL: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        data_pt_int[i]=data_pt[i];
        data_phi_int[i]=data_phi[i];
        data_eta_int[i]=data_eta[i];
    }

    // calc signed components first
    pxy_t met_x[NPART];
    pxy_t met_y[NPART];
    LOOP_PROJECT: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        // Get x, y components
        ProjX(data_pt_int[i], data_phi_int[i], met_x[i]);
        ProjY(data_pt_int[i], data_phi_int[i], met_y[i]);
    }


    pxy_t sum_x = 0;
    pxy_t sum_y = 0;
    SUM: for(int i=0; i<NPART;i++){
#pragma HLS UNROLL
        // Add to x, y sums
        sum_x -= met_x[i];
        sum_y -= met_y[i];
         if(!DEBUG){
             std::cout << "     met x,y = (" << -met_x[i] << ", " << -met_y[i] << ") \t";
             std::cout << " sum x,y = (" << sum_x << ", " << sum_y << ") \t";
             std::cout << " eta  = (" << data_eta[i] << ", ) \t";
             std::cout << " phi  = (" << data_phi[i] << ", ) \t";
             std::cout << " \n";
         }
    }
	
	pxy_t mean_x = -sum_x / NPART;
	pxy_t mean_y = -sum_y / NPART;

	pt2_t dev_x = 0;
	pt2_t dev_y = 0;
	for(int i=0; i<NPART; i++){
		dev_x += ( (met_x[i] - mean_x) * (met_x[i] - mean_x));
		dev_y += ( (met_y[i] - mean_y) * (met_y[i] - mean_y));
	}

	pt2_t sigma2_x = dev_x / (NPART-1);
	pt2_t sigma2_y = dev_y / (NPART-1);

	pt_t sigma_x = hls::sqrt(sigma2_x); 
	pt_t sigma_y = hls::sqrt(sigma2_y);

    res_pt2 = sum_x*sum_x + sum_y*sum_y;

	pt_t res_pt = hls::sqrt(res_pt2);

	pt_t K_value = 1.;
	pt2_t corr_x;
	pt2_t corr_y;

	corr_x = sum_x + (K_value * sigma_x);
	corr_y = sum_y + (K_value * sigma_y);

	corr_pt = corr_x*corr_x + corr_y*corr_y; 
    //PhiFromXY(sum_x,sum_y,res_phi);
    res_phi=0;

	for(int i=0; i<NPART;i++){
		std::cout << "data_pt_int = " << data_pt_int[i] << "\t";
		std::cout << " corr_pt = " << corr_pt << "\n";
	}
	//if data_eta_int

    return;
}

