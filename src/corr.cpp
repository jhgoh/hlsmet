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

    #pragma HLS pipeline ii=54
    
    if(DEBUG) std::cout << "  HW Begin" << std::endl;

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
		if(DEBUG){
			std::cout << "\tjetx = " << (jet_x[i] >> 2) << " \t ";
        	std::cout << "\tjety = " << (jet_y[i] >> 2) << " \n";	
		};
    }

    pxy_t sumJ_x = 0;
    pxy_t sumJ_y = 0;

	pt_t eta13 = int(1.3 * (1<<PHI_SIZE)/6);
	pt_t eta17 = int(1.7 * (1<<PHI_SIZE)/6);
	pt_t eta25 = int(2.5 * (1<<PHI_SIZE)/6);
	pt_t eta30 = int(3.0 * (1<<PHI_SIZE)/6); // **** eta30 was overflowed, so I changed all the type of eta values.


    SUM: for(int i=0; i<NPART;i++){
	#pragma HLS UNROLL

		// **** initial value for var1, var2
		pxy_t var1 = 0;
		pxy_t var2 = 0;

		// Minus eta range is not done.
        if(DEBUG) std::cout << "\tr eta "<< jet_eta[i];

		if ( hls::abs(jet_eta[i]) < eta13 ){
			var1 = int(0.106 * (1<<PT_SIZE) * (1<<PT_DEC_BITS));
			var2 = int(6.6   * (1<<PT_SIZE) * (1<<PT_DEC_BITS));
    
			if(DEBUG) std::cout << "\tCase 1\t";
		}

		if ( hls::abs(jet_eta[i]) < eta17 && hls::abs(jet_eta[i]) > eta13 ){ // || 
				/*(jet_eta[i] < int(-1.3 * (1<<PHI_SIZE)/6)) &&*/
				/*jet_eta[i] > int(-1.7 * (1<<PHI_SIZE)/6) */ 
			var1 = int(0.216 * (1<<PT_SIZE) * (1<<PT_DEC_BITS));
			var2 = int(5.5   * (1<<PT_SIZE) * (1<<PT_DEC_BITS));

          
			if(DEBUG) std::cout << "\tCase 2\t";
		}

		if ( hls::abs(jet_eta[i]) < eta30 && hls::abs(jet_eta[i]) > eta25 ){
			var1 = int(0.083 * (1<<PT_SIZE) * (1<<PT_DEC_BITS));
			var2 = int(13.2  * (1<<PT_SIZE) * (1<<PT_DEC_BITS));

      
			if(DEBUG) std::cout << "\tCase 3\t";
		}

		// **** Here's the problem point. var1 & var2 was shifted 16 bits.
		// **** The jet_x & jet_y was shifted 2 bits from TB, so var1 * jet_x(or y) had 18 bits(16bits from var1 & 2 bits from jet_x).
		// **** But var2 had 16bits. So if we try to sum var1*jet_x(which was 18 bits) and var2(which was 16 bits), there comes problem.
		// **** So I solved the issue by decreasing the number of bits for jet_x bits. 
		// **** However, I am not sure if this is the optimal approach in terms of precision.
		sumJ_x += (var1*(jet_x[i] >> 2)+var2); // **** sumJ_x {16bits shifted}
        sumJ_y += (var1*(jet_y[i] >> 2)+var2); // **** sumJ_y {16bits shifted}

        if(DEBUG){
            std::cout << " sigma(Jet) x,y = (" << (sumJ_x >> 16) << ", " << (sumJ_y>>16) << ") \t";
            std::cout << " \n";
        }
    }


	// correctd MET
	// K factor, HW unit convert?
	pt_t K_value = int(1.5 * (1<<PT_DEC_BITS)); // **** K_value {2 bits shifted}
	pt2_t corr_x;
	pt2_t corr_y;

	// Move the bits of sumJ_x and y. For the 0.106 * (1<<PT_SIZE) ... ?

	// **** This is the same problem with sumJ.
	// **** Because K_value was shifted 2 bits and sumJ was shifted 16bits, so K_value * (sumJ >> PT_SIZE) had 4 bits shifted.
	// **** But met value was shifted 2 bits from TB, so I solved them with same way.
    corr_x = met_x + (K_value * (sumJ_x >> (PT_SIZE + PT_DEC_BITS))); // **** corr_x {2 bits shifted}
	corr_y = met_y + (K_value * (sumJ_y >> (PT_SIZE + PT_DEC_BITS))); // **** corr_y {2 bits shifted}

	// **** Because corr value was shifted 2 bits, corr^2 was shifted 4 bits .
	corr_pt2 = corr_x*corr_x + corr_y*corr_y;  // **** corr_pt2 {4 bits shifted}
    //PhiFromXY(sum_x,sum_y,res_phi);

	std::cout << " corr_pt2 (HW unit) = " << (corr_pt2>>4) << " \n";


    return;
}

