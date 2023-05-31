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
    pt_t in_pt_hw;
    phi_t in_phi_hw, out_phi_hw;

    pt2_t out_pt2_hw;
	pt2_t out_corr_hw;
    
	float in_pt, in_phi, in_jet_eta[NPART], in_jet_pt[NPART], in_jet_phi[NPART];
    float out_pt, out_phi; 
	double out_corr_ref;

    pt_t in_J_pt_hw[NPART];
	phi_t in_J_phi_hw[NPART];
	phi_t in_J_eta_hw[NPART];


    //setup random number generator
    std::default_random_engine generator(1776); // seed
    // random pt uniformly distributed between 10 and 100 GeV for each particle
    std::uniform_real_distribution<float> pt_dist(10.,100.); 
    // random uniform phi
    std::uniform_real_distribution<float> phi_dist(-M_PI,M_PI);

	// random jet
	std::uniform_real_distribution<float> jet_pt_dist(10.,100.);
	// random jet phi
	std::uniform_real_distribution<float> jet_phi_dist(-M_PI,M_PI);
    // random ejt eta
    std::uniform_real_distribution<float> jet_eta_dist(-3.,3.);

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

            etavals[i][j] = jet_eta_dist(generator);
        }
    }
	std::cout<<"Random generate Flag" <<  std::endl;

    //write results to file
    FILE *f;
    f=fopen("results.txt","w");
    FILE *f2;
    f2=fopen("results2.txt","w");

    for (int i=0; i<NTEST; ++i) {
        if(DEBUG) std::cout << "\n\n\n\nEvent " << i << std::endl;

        // convert float to hw units
        in_pt_hw  = pt_t(vals[i][0].first); // 0.25 GeV precision {2bits shifted}
        in_phi_hw = phi_t(vals[i][0].second); // {10bits shifted}
        
        // keep test vals as float
        in_pt  = vals[i][0].first;
        in_phi = vals[i][0].second;

		for(int j=0; j<NPART; j++){
			// convert Jet float to hw units
            in_J_eta_hw[j] = phi_t(etavals[i][j]); // **** {10bits shifted}
            in_J_pt_hw[j]  = pt_t(Jetvals[i][j].first); // 0.25 GeV precision **** {2bits shifted}
            in_J_phi_hw[j] = phi_t(Jetvals[i][j].second); // **** {10bits shifted}

			// Jet test vals
            in_jet_pt[j]  = Jetvals[i][j].first;
            in_jet_phi[j] = Jetvals[i][j].second;
            in_jet_eta[j] = etavals[i][j];
            if(DEBUG && i%100 == 0){
                std::cout << " \t part pt " << in_pt;
                std::cout << "\t phi " << in_phi;
                std::cout << "\t eta " << in_jet_eta[j]<<"\n";
				std::cout << "\t jet pt " << in_jet_pt[j];
				std::cout << "\t jet phi " << in_jet_phi[j];
                std::cout << std::endl;
            }
        }

        out_pt=0.; out_phi=0.; out_pt2_hw=0.; out_phi_hw=0.; // Not use
		//out_corr_hw=0.; out_corr_ref=0.;
        
        // run reference alg
        corr_ref(in_pt, in_phi, in_jet_pt, in_jet_phi, in_jet_eta, out_corr_ref);

        // run HW alg
        corr_hw(in_pt_hw, in_phi_hw, in_J_pt_hw, in_J_phi_hw, in_J_eta_hw, out_corr_hw);


        std::cout << " REF : in METpt = " << in_pt << ", METphi = "<< in_phi << " corrMET = " << out_corr_ref<<"\n";
        // for HW alg, convert back to nice units for printingM
        int out_phi_hw_int = float(out_phi_hw);
        float out_phi_hw_rad = float(out_phi_hw) * (2*M_PI)/(1<<PHI_SIZE);
        float out_pt_hw = sqrt(float(out_pt2_hw)) / (1<<PT_DEC_BITS); // 0.25GeV to GeV  // Not use

        // **** We have to sqrt the restored out_corr_hw value, not the 4 bits shifted out_corr_hw.
		// float out_co_hw = sqrt(float(out_corr_hw)) / (1<<(PT_DEC_BITS)); // Corrected MET // 0.25GeV to GeV
        float out_co_hw = sqrt(float(out_corr_hw)); // Corrected MET // 0.25GeV to GeV 

        std::cout << "  HW  : in METpt = " << float(in_pt_hw) << ", METphi = "<< in_phi_hw << ", corrMET = "<< out_co_hw<< "\n";

        //if not debugging the full event details, print a compact output (in nice units)
        if(DEBUG && i%50 == 0){
            std::cout << "Event " << i << "\n";
			std::cout << " corrMET (REF vs HW) " << out_corr_ref << " vs " << out_co_hw <<"\n";
		}
        fprintf(f, "%f %f %f %f \n", out_corr_ref, out_phi, out_co_hw, out_phi_hw_rad);
        fprintf(f2, "%f %f %f \n", out_corr_ref, out_co_hw, in_pt);
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
