#include "generateEvents.h"
#include <random>

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

