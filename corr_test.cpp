/*
MET calculation from PF objects
*/
#include <vector>
#include <cstdio>
#include <utility>
#include <random>

#include "ap_int.h"
#include "ap_fixed.h"
#include "src/corr.h"

int alg_test() {
    // The conventional test algorithm that compares the final outputs 
    // of HLS and floating point calculations

    // calculate met for NPART particles
    pt_t in_pt_hw[NPART];
    pt2_t out_pt2_hw;
	pt2_t out_corr_hw;
    phi_t in_phi_hw[NPART], out_phi_hw;
	phi_t in_eta_hw[NPART];
    float in_pt[NPART], in_phi[NPART], in_eta[NPART], in_jet_pt[NPART], in_jet_phi[NPART];
    float out_pt, out_phi, out_corr;

    pt_t in_J_pt_hw[NPART];
	phi_t in_J_phi_hw[NPART];


    //setup random number generator
    std::default_random_engine generator(1776); // seed
    // random pt uniformly distributed between 10 and 100 GeV for each particle
    std::uniform_real_distribution<float> pt_dist(10.,100.); 
    // random uniform phi
    std::uniform_real_distribution<float> phi_dist(-M_PI,M_PI);
    // random uniform eta
    std::uniform_real_distribution<float> eta_dist(-4.,4.);

	// random jet
	std::uniform_real_distribution<float> jet_pt_dist(10.,100.);
	// random jet phi
	std::uniform_real_distribution<float> jet_phi_dist(-M_PI,M_PI);

    // fill random test data
    std::vector<std::vector<std::pair<float,float> > > vals; 
    std::vector<std::vector<std::pair<float,float> > > Jetvals; 
	//std::vector<std::vector<std::pair<float,float> > > etavals; 
    std::vector<std::vector<float> > etavals; 
    // Dimensions: #events=NTEST x #particles=NPART
    // type is a pair: (pt,phi)
    vals.resize(NTEST);
	Jetvals.resize(NTEST);
    etavals.resize(NTEST);
    for(int i=0; i<NTEST; i++){
        vals[i].resize(NPART);
		Jetvals[i].resize(NPART);
        etavals[i].resize(NPART);
        for(int j=0; j<NPART; j++){
            vals[i][j].first  = pt_dist(generator);
            vals[i][j].second = phi_dist(generator);

			Jetvals[i][j].first  = jet_pt_dist(generator);
			Jetvals[i][j].second = jet_phi_dist(generator);

            etavals[i][j] = eta_dist(generator);
        }
    }
	std::cout<<"Flag" <<  std::endl;

    //write results to file
    FILE *f;
    f=fopen("results.txt","w");
    FILE *f2;
    f2=fopen("results2.txt","w");

    for (int i=0; i<NTEST; ++i) {
        if(DEBUG) std::cout << "\n\n\n\nEvent " << i << std::endl;
        for(int j=0; j<NPART; j++){
            // convert float to hw units
            in_pt_hw[j]  = vals[i][j].first * (1<<PT_DEC_BITS); // 0.25 GeV precision
            in_phi_hw[j] = int(vals[i][j].second * (1<<PHI_SIZE)/(2*M_PI));
            in_eta_hw[j] = int(etavals[i][j] * (1<<PHI_SIZE));
            in_J_pt_hw[j]  = Jetvals[i][j].first * (1<<PT_DEC_BITS); // 0.25 GeV precision
            in_J_phi_hw[j] = int(Jetvals[i][j].second * (1<<PHI_SIZE)/(2*M_PI));

            // keep test vals as float
            in_pt[j]  = vals[i][j].first;
            in_phi[j] = vals[i][j].second;
            in_jet_pt[j]  = Jetvals[i][j].first;
            in_jet_phi[j] = Jetvals[i][j].second;
            in_eta[j] = etavals[i][j];
            if(i%100 == 0){
                std::cout << " \t part pt " << in_pt[j];
                std::cout << "\t phi " << in_phi[j];
                std::cout << "\t eta " << in_eta[j]<<"\n";
				std::cout << "\t jet pt " << in_jet_pt[j];
				std::cout << "\t jet phi " << in_jet_phi[j];
                std::cout << std::endl;
            }
        }
        out_pt2_hw=0.; out_phi_hw=0.; out_corr_hw=0.;
        out_pt=0.; out_phi=0.; out_corr=0.;
        
        // run reference alg
        corr_ref(in_pt, in_phi, in_eta, in_jet_pt, in_jet_phi, out_pt, out_phi, out_corr);

        // run HW alg
        corr_hw(in_pt_hw, in_phi_hw, in_eta_hw, in_J_pt_hw, in_J_phi_hw, out_pt2_hw, out_phi_hw, out_corr_hw);


        if(DEBUG) std::cout << " REF : met(pt = " << out_pt << ", phi = "<< out_phi << "corrMET = " << out_corr<<")\n";
        // for HW alg, convert back to nice units for printing
        int out_phi_hw_int = float(out_phi_hw);
        float out_phi_hw_rad = float(out_phi_hw) * (2*M_PI)/(1<<PHI_SIZE);
        float out_pt_hw = sqrt(float(out_pt2_hw)) / (1<<PT_DEC_BITS); // 0.25GeV to GeV
		float out_co_hw = sqrt(float(out_corr_hw)) / (1<<PT_DEC_BITS);
        std::cout << "  HW : met(pt = " << out_pt_hw << ", phi = "<< out_phi_hw_rad << ", corrMET = "<< out_co_hw<< ")\n";

        //if not debugging the full event details, print a compact output (in nice units)
        if(true && !DEBUG && NTEST<=100)
            std::cout << "Event " << i
                      << " (REF vs HW) met " << out_pt << " vs " << out_pt_hw
                      << ", phi "<< out_phi << " vs "<< out_phi_hw_rad << "\n";
			std::cout << " corrMET (REF vs HW) " << out_corr << " vs " << out_co_hw <<"\n";
        fprintf(f, "%f %f %f %f \n", out_pt, out_phi, out_pt_hw, out_phi_hw_rad);
        fprintf(f2, "%f %f %f %f \n", out_corr, out_phi, out_co_hw, out_phi_hw_rad);
        //fprintf(f, "%f %f %f %d \n", out_pt, out_phi, out_pt_hw, out_phi_hw_int);
    }
    fclose(f);
    fclose(f2);

    return 0;
}


int main() {

    // test the algorithm
    alg_test();

    return 0;
}
