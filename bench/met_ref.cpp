/*
Reference implementation of MET calculation from PF objects
*/
#include "met_ref.h"
#include <cmath>

void met_ref(float in_pt[NPART], float in_phi[NPART], float& out_pt, float& out_phi)
{
  double met_x = 0.;
  double met_y = 0.;
  for(int i=0; i<NPART; ++i) {
    met_x -= in_pt[i] * std::cos(in_phi[i]);
    met_y -= in_pt[i] * std::sin(in_phi[i]);
  }

  out_pt = std::hypot(met_x, met_y);
  out_phi = std::atan2(met_y, met_x);

  return;
}

