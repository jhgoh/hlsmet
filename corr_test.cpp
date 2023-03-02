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
    float in_pt[NPART], in_phi[NPART], in_eta[NPART];
    float out_pt, out_phi;

    //setup random number generator
    std::default_random_engine generator(1776); // seed
    std::uniform_real_distribution<float> pt_dist(10.,100.); 
    // random pt uniformly distributed between 10 and 100 GeV for each particle
    std::uniform_real_distribution<float> phi_dist(-M_PI,M_PI);
    // random uniform phi
    std::uniform_real_distribution<float> eta_dist(-4.,4.);
    // random uniform eta

    // fill random test data
    std::vector<std::vector<std::pair<float,float> > > vals; 
	//std::vector<std::vector<std::pair<float,float> > > etavals; 
    std::vector<std::vector<float> > etavals; 
    // Dimensions: #events=NTEST x #particles=NPART
    // type is a pair: (pt,phi)
    vals.resize(NTEST);
    etavals.resize(NTEST);
    for(int i=0; i<NTEST; i++){
        vals[i].resize(NPART);
        etavals[i].resize(NPART);
        for(int j=0; j<NPART; j++){
            vals[i][j].first  = pt_dist(generator);
            vals[i][j].second = phi_dist(generator);
            //etavals[i][j].first = eta_dist(generator);
            etavals[i][j] = eta_dist(generator);
        }
    }
	std::cout<<"Flag" <<  std::endl;

    //write results to file
    FILE *f;
    f=fopen("results.txt","w");

    for (int i=0; i<NTEST; ++i) {
        if(DEBUG) std::cout << "\n\n\n\nEvent " << i << std::endl;
        for(int j=0; j<NPART; j++){
            // convert float to hw units
            in_pt_hw[j]  = vals[i][j].first * (1<<PT_DEC_BITS); // 0.25 GeV precision
            in_phi_hw[j] = int(vals[i][j].second * (1<<PHI_SIZE)/(2*M_PI));
            in_eta_hw[j] = int(etavals[i][j] * (1<<PHI_SIZE));
            //in_eta_hw[j] = int(etavals[i][j].first * (1<<PT_DEC_BITS));
            // keep test vals as float
            in_pt[j]  = vals[i][j].first;
            in_phi[j] = vals[i][j].second;
            in_eta[j] = etavals[i][j];
            //in_eta[j] = etavals[i][j].first;
            if(i%100 == 0){
                std::cout << " \t part pt " << in_pt[j];
                std::cout << "\t phi " << in_phi[j];
                std::cout << "\t eta " << in_eta[j];
                std::cout << std::endl;
            }
        }
        out_pt2_hw=0.; out_phi_hw=0.; out_corr_hw=0.;
        out_pt=0.; out_phi=0.;
        
        // run reference alg
        corr_ref(in_pt, in_phi, in_eta, out_pt, out_phi);

        // run HW alg
        corr_hw(in_pt_hw, in_phi_hw, in_eta_hw, out_pt2_hw, out_phi_hw, out_corr_hw);


        if(DEBUG) std::cout << " REF : met(pt = " << out_pt << ", phi = "<< out_phi << ")\n";
        // for HW alg, convert back to nice units for printing
        int out_phi_hw_int = float(out_phi_hw);
        float out_phi_hw_rad = float(out_phi_hw) * (2*M_PI)/(1<<PHI_SIZE);
        float out_pt_hw = sqrt(float(out_pt2_hw)) / (1<<PT_DEC_BITS); // 0.25GeV to GeV
		float out_co_hw = float(out_corr_hw) / (1<<PT_DEC_BITS);
        if(DEBUG) std::cout << "  HW : met(pt = " << out_pt_hw << ", phi = "<< out_phi_hw_rad << ")\n";

        //if not debugging the full event details, print a compact output (in nice units)
        if(true && !DEBUG && NTEST<=100)
            std::cout << "Event " << i
                      << " (REF vs HW) met " << out_pt << " vs " << out_pt_hw
                      << ", phi "<< out_phi << " vs "<< out_phi_hw_rad << "\n";
        fprintf(f, "%f %f %f %f \n", out_pt, out_phi, out_pt_hw, out_phi_hw_rad);
        //fprintf(f, "%f %f %f %d \n", out_pt, out_phi, out_pt_hw, out_phi_hw_int);
    }
    fclose(f);

    return 0;
}


int main() {

    // test the algorithm
    alg_test();

    return 0;
}
