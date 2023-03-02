/*
HLS implementation of MET calculation from PF objects
*/
#include "corr.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
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

    res_pt2 = sum_x*sum_x + sum_y*sum_y;
    //PhiFromXY(sum_x,sum_y,res_phi);
    res_phi=0;

	phi_t eta1 = 1.3*(1<<PHI_SIZE);
	float p1 = 0.106*(1<<PT2_SIZE)*(1<<PT_DEC_BITS);
	float p2 = 6.6*(1<<PT2_SIZE)*(1<<PT_DEC_BITS);
	std::cout << "p " << eta1 << "\t" << p1 << "\t" << p2 << "\n";
	pt_t para1 = int(p1);
	pt_t para2 = int(p2);
	std::cout << "para " << eta1 << "\t" << para1 << "\t" << para2 << "\n";
	for(int i=0; i<NPART;i++){
		if (data_eta_int[i]<eta1){
			std::cout << para1*data_pt_int[i] << "\t";
			std::cout << (para1*data_pt_int[i])+para2 << "\t";
			std::cout << data_pt_int[i]*((para1*data_pt_int[i])+para2) << "\n";
			corr_pt[i] = data_pt_int[i]*((para1*data_pt_int[i])+para2);
			std::cout << "data_pt_int = " << data_pt_int[i] << "\t";
			std::cout << " corr_pt = " << corr_pt[i] << "\n";
		}
	}
	//if data_eta_int

    return;
}

