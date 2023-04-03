#include <math.h>
#include <cmath>
#include <algorithm>
#include "ap_int.h"
#include "ap_fixed.h"
#include "src/corr.h"

void corr_ref(float in_pt, float in_phi, float in_jet_pt[NPART], float in_jet_phi[NPART], float in_jet_eta[NPART], float out_corr){
    if(DEBUG) std::cout << " REF Begin" << std::endl;
    double met_x=0.;
    double met_y=0.;
    double jet_x[NPART];
    double jet_y[NPART];
	double sumJ_x=0.;
	double sumJ_y=0.;
	met_x -= in_pt * cos(in_phi);
	met_y -= in_pt * sin(in_phi);
    if(DEBUG) std::cout << "     metx = " << met_x << " \t ";
    if(DEBUG) std::cout << "mety = " << met_y << " \n";
    for(int i=0;i<NPART;i++){
        jet_x[i] = in_jet_pt[i] * cos(in_jet_phi[i]);
        jet_y[i] = in_jet_pt[i] * sin(in_jet_phi[i]);
        if(DEBUG) std::cout << "     jetx = " << jet_x[i] << " \t ";
        if(DEBUG) std::cout << "jety = " << jet_y[i] << " \n";
    }
	for(int i=0; i<NPART;i++){
		std::cout << "jet eta = " << in_jet_eta[i] << "\n";
		if ( in_jet_eta[i] < 1.3 && in_jet_eta[i] > -1.3 ){
			sumJ_x += 0.106*jet_x[i]+6.6;
			sumJ_y += 0.106*jet_y[i]+6.6;
		}
		if ( (in_jet_eta[i] > 1.3 && in_jet_eta[i] < 1.7) 
//				|| (in_jet_eta[i] < -1.3 && in_jet_eta[i] > -1.7) 
				){
			sumJ_x += 0.216*jet_x[i]+5.5;
			sumJ_y += 0.216*jet_y[i]+5.5;
		}
		if ( (in_jet_eta[i] > 2.5 && in_jet_eta[i] < 3.) 
//				|| (in_jet_eta[i] < -2.5 && in_jet_eta[i] > -3.) 
				){
			sumJ_x += 0.083*jet_x[i]+13.2;
			sumJ_y += 0.083*jet_y[i]+13.2;
		}
	}


	double K_val = 1.5;

	double corr_x = met_x + (K_val * sumJ_x);
	double corr_y = met_y + (K_val * sumJ_y);

	std::cout << "corr_x, y = " << corr_x << ", " << corr_y <<"\n";

	out_corr = sqrt(corr_x*corr_x + corr_y*corr_y);

	std::cout << "	out_corr = " << out_corr << "\n";
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

