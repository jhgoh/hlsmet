#ifndef CORRMET_H
#define CORRMET_H

void corrmet_hw(pt_t in_met, phi_t in_metphi,
                pt_t in_pt[NJET], phi_t in_phi[NJET], phi_t in_eta[NJET],
                pt_t& out_corrmet, phi_t& out_corrmetphi);

#endif
