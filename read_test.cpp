/*
MET calculation from PF objects
*/
#include <vector>
#include <cstdio>
#include <utility>
#include <random>
#include <istream>
#include <fstream>

#include "ap_int.h"
#include "ap_fixed.h"
#include "hls_stream.h"
#include "src/corr.h"
#include "DiscretePFInputs.h"

#define FLOATPI 3.141593

typedef ap_uint<64> word_t;
typedef ap_int<16> var_t;
typedef struct { word_t data[NPART];} PFInputWords;

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

	std::string baseDir = "/home/jhong/PFMET/hlsmet";
	std::ifstream infile(baseDir+"/TTbar_1000evt_54part_v2.dump");
	std::string line;
	std::vector<std::vector<word_t> > word_list;
	while(std::getline(infile, line)){
		//if(line[0]=='#') continue; // comment handling
		std::vector<word_t> words;
		std::istringstream iss(line);
		for(std::string s; iss >> s; ){
			word_t w( s.c_str(), 16 /*hex*/ );
			//std::cout << s << " " << w.to_string(16) << std::endl;
			words.push_back(w);
		}
		words.resize(NPART); // add zeros or truncate as desired
		word_list.push_back(words);
	}
	word_list.resize(NTEST); // add zeros or truncate as desired

	for (int i=0; i<NTEST; ++i){
		for (int j=0; j<NPART; j++){
			if(DEBUG) std::cout<<"reading "<< word_list[i][j] << "\n";
		}
	}

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

    //write results to file
    FILE *f;
    f=fopen("results.txt","w");
    FILE *f2;
    f2=fopen("results2.txt","w");

	float in_pt_dump[NPART], in_phi_dump[NPART], in_eta_dump[NPART];
    pt_t inDump_pt_hw;
    phi_t inDump_phi_hw;

    pt_t inDump_J_pt_hw[NPART];
	phi_t inDump_J_phi_hw[NPART];
	phi_t inDump_J_eta_hw[NPART];
	float inDump_pt, inDump_phi, inDump_jet_eta[NPART], inDump_jet_pt[NPART], inDump_jet_phi[NPART];

    for (int i=0; i<NTEST; ++i) {
		if(DEBUG) std::cout << "\n\n\n\nEvent " << i << std::endl;
		// pack words into stream
		hls::stream<PFInputWords> input_stream;
		PFInputWords input_array;
		for(int j=0; j<NPART+1; j++){
			input_array.data[j] = word_list[i][j];
			if (j == 0){
				in_pt_dump[j] = word_list[i][j](63,48);
				if(in_pt_dump[j] > (1<<16)) in_pt_dump[j] = in_pt_dump[j] - (1<<16);
				in_pt_dump[j] = in_pt_dump[j];// / (1<<2);
				inDump_pt_hw = pt_t(in_pt_dump[j]);
				inDump_pt = in_pt_dump[j]; 
				if(DEBUG) std::cout<<j<<" input 0th pt, phi "<<in_pt_dump[j]<<", ";
			
				in_phi_dump[j] = word_list[i][j](47,32);
				if(in_phi_dump[j] > (1<<16)) in_phi_dump[j] = in_phi_dump[j] - (1<<16);
				in_phi_dump[j] = in_phi_dump[j] * (2*FLOATPI)/(1<<10) / (4*(180/FLOATPI));
				inDump_phi_hw = phi_t(in_phi_dump[j]);
				inDump_phi = in_phi_dump[j];
				if(DEBUG) std::cout<<j<<" input 0th pt, phi "<<in_phi_dump[j]<<"\n ";
			}
			else {
				in_pt_dump[j] = word_list[i][j](63,48);
				if(in_pt_dump[j] > (1<<16)) in_pt_dump[j] = in_pt_dump[j] - (1<<16);
				in_pt_dump[j] = in_pt_dump[j];// / (1<<2);
				inDump_J_pt_hw[j-1] = pt_t(in_pt_dump[j]);
				inDump_jet_pt[j-1] = in_pt_dump[j];
				if(DEBUG) std::cout<<j<<" input pt, phi, eta "<<in_pt_dump[j]<<", ";

				in_phi_dump[j] = word_list[i][j](47,32);
				if(in_phi_dump[j] > (1<<16)) in_phi_dump[j] = in_phi_dump[j] - (1<<16);
				in_phi_dump[j] = in_phi_dump[j] * (2*FLOATPI)/(1<<10) / (4*(180/FLOATPI));
				inDump_J_phi_hw[j-1] = phi_t(in_phi_dump[j]);
				inDump_jet_phi[j-1] = in_phi_dump[j];
				if(DEBUG) std::cout<<in_phi_dump[j]<<", ";

				in_eta_dump[j] = word_list[i][j](31,16);
				if(in_eta_dump[j] > (1<<16)) in_eta_dump[j] = in_eta_dump[j] - (1<<16);
				in_eta_dump[j] = in_eta_dump[j] *(2*FLOATPI) / (1<<10) / (4*(180/FLOATPI));
				inDump_J_eta_hw[j-1] = phi_t(in_eta_dump[j]);
				inDump_jet_eta[j-1] = in_eta_dump[j];
				if(DEBUG) std::cout<<in_eta_dump[j]<<", \n";
			}
		}
        // convert float to hw units
        in_pt_hw  = pt_t(vals[i][0].first); // 0.25 GeV precision {2bits shifted}
        in_phi_hw = phi_t(vals[i][0].second); // {10bits shifted}
        
        // keep test vals as float
        in_pt  = vals[i][0].first;
        in_phi = vals[i][0].second;

		if(DEBUG){
                std::cout << " \t part pt " << inDump_pt;
                std::cout << "\t phi " << inDump_phi<<"\n";
		}
		for(int j=0; j<NPART; j++){
			// convert Jet float to hw units
            in_J_eta_hw[j] = phi_t(etavals[i][j]); // **** {10bits shifted}
            in_J_pt_hw[j]  = pt_t(Jetvals[i][j].first); // 0.25 GeV precision **** {2bits shifted}
            in_J_phi_hw[j] = phi_t(Jetvals[i][j].second); // **** {10bits shifted}

			// Jet test vals
            in_jet_pt[j]  = Jetvals[i][j].first;
            in_jet_phi[j] = Jetvals[i][j].second;
            in_jet_eta[j] = etavals[i][j];
            if(DEBUG){
                std::cout << "\t eta (ref, hw) " << inDump_jet_eta[j] << ", " << inDump_J_eta_hw[j];
				std::cout << "\t jet pt (ref, hw) " << inDump_jet_pt[j] << ", " << inDump_J_pt_hw[j];
				std::cout << "\t jet phi (ref), hw " << inDump_jet_phi[j] << ", " << inDump_J_phi_hw[j] << "\n";
                std::cout << std::endl;
            }
        }

        out_pt=0.; out_phi=0.; out_pt2_hw=0.; out_phi_hw=0.; // Not use
		//out_corr_hw=0.; out_corr_ref=0.;
        
        // run reference alg
        //corr_ref(in_pt, in_phi, in_jet_pt, in_jet_phi, in_jet_eta, out_corr_ref);

        // run HW alg
        //corr_hw(in_pt_hw, in_phi_hw, in_J_pt_hw, in_J_phi_hw, in_J_eta_hw, out_corr_hw);

        // run reference alg
        corr_ref(inDump_pt, inDump_phi, inDump_jet_pt, inDump_jet_phi, inDump_jet_eta, out_corr_ref);

        // run HW alg
        corr_hw(inDump_pt_hw, inDump_phi_hw, inDump_J_pt_hw, inDump_J_phi_hw, inDump_J_eta_hw, out_corr_hw);

        //std::cout << " REF : in METpt = " << in_pt << ", METphi = "<< in_phi << " corrMET = " << out_corr_ref<<"\n";
        std::cout << " REF : in METpt = " << inDump_pt << ", METphi = "<< inDump_phi << " corrMET = " << out_corr_ref<<"\n";
        // for HW alg, convert back to nice units for printingM
        int out_phi_hw_int = float(out_phi_hw);
        float out_phi_hw_rad = float(out_phi_hw) * (2*M_PI)/(1<<PHI_SIZE);
        float out_pt_hw = sqrt(float(out_pt2_hw)) / (1<<PT_DEC_BITS); // 0.25GeV to GeV  // Not use

        // **** We have to sqrt the restored out_corr_hw value, not the 4 bits shifted out_corr_hw.
		// float out_co_hw = sqrt(float(out_corr_hw)) / (1<<(PT_DEC_BITS)); // Corrected MET // 0.25GeV to GeV
        float out_co_hw = sqrt(float(out_corr_hw)); // Corrected MET // 0.25GeV to GeV 

        //std::cout << "  HW  : in METpt = " << float(in_pt_hw) << ", METphi = "<< in_phi_hw << ", corrMET = "<< out_co_hw<< "\n";
        std::cout << "  HW  : in METpt = " << float(inDump_pt_hw) << ", METphi = "<< inDump_phi_hw << ", corrMET = "<< out_co_hw<< "\n";

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
