#include <math.h>
#include <cmath>
#include <algorithm>
#include "ap_int.h"
#include "ap_fixed.h"
#include "src/corr.h"

void corr_ref(float in_pt[NPART], float in_phi[NPART], float in_eta[NPART], float in_jet_pt[NPART], float in_jet_phi[NPART], float& out_pt, float& out_phi, float& out_corr){
    if(DEBUG) std::cout << " REF Begin" << std::endl;
    double met_x=0.;
    double met_y=0.;
    double jet_x[NPART];
    double jet_y[NPART];
	double sumJ_x=0.;
	double sumJ_y=0.;
    for(int i=0;i<NPART;i++){
        met_x -= in_pt[i] * cos(in_phi[i]);
        met_y -= in_pt[i] * sin(in_phi[i]);
        jet_x[i] = in_jet_pt[i] * cos(in_jet_phi[i]);
        jet_y[i] = in_jet_pt[i] * sin(in_jet_phi[i]);
        if(DEBUG) std::cout << "     in eta = " << in_eta[i] << " \n ";
        if(DEBUG) std::cout << "     metx = " << met_x << " \t ";
        if(DEBUG) std::cout << "mety = " << met_y << " \n";
    }
	for(int i=0; i<NPART;i++){
		sumJ_x += jet_x[i];
		sumJ_y += jet_y[i];
	}
	double mean_x = sumJ_x / NPART;
	double mean_y = sumJ_y / NPART;

	double dev_x = 0.;
	double dev_y = 0.;
	for(int i=0; i<NPART; i++){
		dev_x += ( (jet_x[i] - mean_x) * (jet_x[i] - mean_x) );
		dev_y += ( (jet_y[i] - mean_y) * (jet_y[i] - mean_y) );
	}
	double sigma_x = sqrt(dev_x/(NPART-1));
	double sigma_y = sqrt(dev_y/(NPART-1));

    double met2 = pow(met_x,2)+pow(met_y,2);
    out_pt = sqrt(met2);
    if(out_pt<1e-10) out_pt = 1e-10; // guard divide by zero
    out_phi = met_y>=0 ? acos(met_x/out_pt) : -acos(met_x/out_pt);

	double K_val = 1.;

	double corr_x = met_x + (K_val * sigma_x);
	double corr_y = met_y + (K_val * sigma_y);

	out_corr = sqrt(corr_x*corr_x + corr_y*corr_y);

    if(DEBUG){
        std::cout << "     x/tot = " << met_x/out_pt << " \t ";
        std::cout << "(x/tot)^2 = " << pow(met_x/out_pt,2) << " \t ";
        std::cout << "acos(x/tot) = " << acos(met_x/out_pt) << " \t ";
        std::cout << "rotated = " << out_phi << " \n ";

        std::cout << "     met = " << out_pt << " \t ";
        std::cout << "met2 = " << met2 << " \t ";
        std::cout << "phi = " << out_phi << " \n";
    }
    return;
}

