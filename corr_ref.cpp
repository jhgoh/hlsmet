#include <math.h>
#include <cmath>
#include <algorithm>
#include "ap_int.h"
#include "ap_fixed.h"
#include "src/corr.h"

void corr_ref(float in_pt, float in_phi, float in_jet_pt[NPART], float in_jet_phi[NPART], float in_jet_eta[NPART], double& out_corr){
    if(DEBUG) std::cout << " REF Begin" << std::endl;
    
	// MET project
	float met_x;
    float met_y;
	met_x = in_pt * cos(in_phi);
	met_y = in_pt * sin(in_phi);
    if(DEBUG) std::cout << "\tmetx = " << met_x << " \t ";
    if(DEBUG) std::cout << "\tmety = " << met_y << " \n\n";

	// Jet project
    float jet_x[NPART];
    float jet_y[NPART];
    for(int i=0;i<NPART;i++){
        jet_x[i] = in_jet_pt[i] * cos(in_jet_phi[i]);
        jet_y[i] = in_jet_pt[i] * sin(in_jet_phi[i]);
        if(DEBUG) std::cout << "\tjetx = " << jet_x[i] << " \t ";
        if(DEBUG) std::cout << "\tjety = " << jet_y[i] << " \n";
    }

	// Sum Jet x, y
	double sumJ_x=0.;
	double sumJ_y=0.;
	
	for(int i=0; i<NPART;i++){
		if(DEBUG) std::cout << "\tr eta "<< in_jet_eta[i];
		if (in_jet_eta[i] < 0) std::cout << "Minus\t" << in_jet_eta[i] << "\n";
		if ( abs(in_jet_eta[i]) < 1.3 ){
			sumJ_x += pow(0.106*jet_x[i]+6.6,2);
			sumJ_y += pow(0.106*jet_y[i]+6.6,2);
			if(DEBUG) std::cout << "\tCase 1\t";
		}
		if ( abs(in_jet_eta[i]) > 1.3 && abs(in_jet_eta[i]) < 1.7 ){
			sumJ_x += pow(0.216*jet_x[i]+5.5,2);
			sumJ_y += pow(0.216*jet_y[i]+5.5,2);
			if(DEBUG) std::cout << "\tCase 2\t";
		}
		if ( abs(in_jet_eta[i]) > 2.5 && abs(in_jet_eta[i]) < 3. ){
			sumJ_x += pow(0.083*jet_x[i]+13.2,2);
			sumJ_y += pow(0.083*jet_y[i]+13.2,2);
			if(DEBUG) std::cout << "\tCase 3\t";
		}
		if(DEBUG) std::cout << " sigma(Jet) x,y = (" << sumJ_x << ", " << sumJ_y << ") \t";
		if(DEBUG) std::cout <<"\n";
	}


	double K_val = 1.5;

	double corr_x = met_x + (K_val * sumJ_x);
	double corr_y = met_y + (K_val * sumJ_y);

	if(DEBUG) std::cout << "corr_x, y = " << corr_x << ", " << corr_y <<"\n";

	//out_corr = sqrt(corr_x*corr_x + corr_y*corr_y);
	double out_corr2 = pow(corr_x, 2) + pow(corr_y,2);
	out_corr = sqrt(out_corr2);
	std::cout << "out_corr2 = " << out_corr2 << "\n";
	std::cout << "out_corr = " << out_corr << "\n";
	//float out_phi
    //if(out_corr<1e-10) out_corr = 1e-10; // guard divide by zero
    //out_phi = met_y>=0 ? acos(met_x/out_corr) : -acos(met_x/out_corr);

    //if(DEBUG){
    //    std::cout << "     x/tot = " << met_x/out_corr << " \t ";
    //    std::cout << "(x/tot)^2 = " << pow(met_x/out_corr,2) << " \t ";
    //    std::cout << "acos(x/tot) = " << acos(met_x/out_corr) << " \t ";
    //    std::cout << "rotated = " << out_phi << " \n ";

    //    std::cout << "corr met = " << out_corr << " \t ";
    //    //std::cout << "met2 = " << met2 << " \t ";
    //    std::cout << "phi = " << out_phi << " \n";
    //}
    return;
}

