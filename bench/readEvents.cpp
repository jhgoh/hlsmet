#include "readEvents.h"
#include "DiscretePFInputs.h"
#include <random>
#include <sstream>
#include <fstream>

#include "ap_int.h"
#include "ap_fixed.h"
typedef ap_uint<64> word_t;

std::vector<eventF_t> generateEvents(const int nEvent, const int nPart)
{
  //setup random number generator
  std::default_random_engine generator(1776); // seed
  // random pt uniformly distributed between 10 and 100 GeV for each particle
  std::uniform_real_distribution<float> pt_dist(10.,100.); 
  // random uniform phi
  std::uniform_real_distribution<float> phi_dist(-M_PI,M_PI);

  // fill random test data
  std::vector<eventF_t> out;
  out.reserve(nEvent);

  for ( int i=0; i<nEvent; ++i ) {
    std::vector<float> pts, phis;
    for ( int j=0; j<nPart; ++j ) {
      pts.push_back(pt_dist(generator));
      phis.push_back(phi_dist(generator));
    }
    out.push_back(std::make_pair(pts, phis));
  }

  return out;
}

std::vector<eventF_t> readEventsFromHex(const std::string fileName, const int nEvent, const int nPart)
{
  std::cout << "* ReadEvent from the file \"" << fileName << "\"\n"
            << "* Note that the input is a dump file converted to HEX format,\n"
            << "* but we are converting them to float numbers\n"
            << "* Therefore, if our pt_t or anything else in the MET calculation\n"
            << "* is compatible with the original HEX file, there should be no\n"
            << "* loss in the information.\n";
  std::vector<eventF_t> events;

  // Read from the input file and fill up the event;
  std::ifstream fin(fileName);
  if ( !fin.is_open() ) {
    std::cout << "Could not open the file..." << std::endl;
    return events;
  }
  for ( int i=0; i<nEvent; ++i ) {
    eventF_t event;

    std::string line;
    if ( !std::getline(fin, line) ) break;
    std::istringstream ss(line);
    std::string tok;
    for ( int j=0; j<nPart; ++j ) {
      if ( !(ss >> tok) ) break;
      word_t w(tok.c_str(), 16);

      int pt_hw = w(63,48);
      int phi_hw = w(47,32);
      int eta_hw = w(31,16); // not used at this moment.

      // Subtract by 2^16 because negative integers are done in the hex file
      if ( pt_hw  > (1<<15) ) pt_hw -= (1<<16);
      if ( phi_hw > (1<<15) ) phi_hw -= (1<<16);
      if ( eta_hw > (1<<15) ) eta_hw -= (1<<16);

      // Have to divided by factors, defined in the "DiscretePFInputs.h",
      // but also divide by yet another factor used during the conversion step.
      // I think the code for the conversion is this one: 
      //    https://github.com/ParticleChef/hlsmet/blob/64words/convertDump.cpp
      const float pt = float(pt_hw)/l1tpf_int::CaloCluster::PT_SCALE;
      const float phi = float(phi_hw)/l1tpf_int::CaloCluster::ETAPHI_SCALE;
      const float eta = float(eta_hw)/l1tpf_int::CaloCluster::ETAPHI_SCALE;

      event.first.push_back(pt);
      event.second.push_back(phi);
    }
    events.emplace_back(event);
  }

  return events;
}

eventHW_t to_HW(const eventF_t& in_event)
{
  eventHW_t out_event;
  const size_t nPart = in_event.first.size();
  for ( size_t i=0; i<nPart; ++i ) {
    const auto pt = in_event.first[i];
    const auto phi = in_event.second[i];
#ifdef INTONLY
    out_event.first.push_back(f2hwPt(pt));
    out_event.second.push_back(f2hwPhi(phi));
#else
    out_event.first.emplace_back(pt);
    out_event.second.emplace_back(phi);
#endif
  }
  return out_event;
}

