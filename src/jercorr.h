#ifndef JERCORR_H
#define JERCORR_H

void jercorrmet_hw(pt_t in_met, phi_t in_metphi,
                   pt_t in_pt[NJET], phi_t in_phi[NJET], phi_t in_eta[NJET],
                   pt_t& out_corrmet, phi_t& out_corrmetphi);

#endif
