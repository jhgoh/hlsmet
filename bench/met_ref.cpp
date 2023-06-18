/*
Reference implementation of MET calculation from PF objects
*/
#include "met_ref.h"
#include <cmath>

void met_ref(float in_pt[NPART], float in_phi[NPART], float& out_pt2, float& out_phi)
{
#if DEBUG==1
  std::cout << "--- REF Begin ---\n";
#endif

  double met_x = 0.;
  double met_y = 0.;
  for(int i=0; i<NPART; ++i) {
    met_x -= in_pt[i] * std::cos(in_phi[i]);
    met_y -= in_pt[i] * std::sin(in_phi[i]);
#if(DEBUG==1)
std::cout << "DEBUG(REF) " << i << ' ' << in_pt[i] << ' ' << in_phi[i] << std::endl;
#endif
  }

  out_pt2 = met_x*met_x + met_y*met_y;
  out_phi = std::atan2(met_y, met_x);
#if DEBUG==1
  std::cout << "REF MET pT2=" << out_pt2 << " Phi=" << out_phi << " px=" << met_x << " py=" << met_y << std::endl;
#endif

  return;
}

