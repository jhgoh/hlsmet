#ifndef REF_GETRANDOMEVENT_H
#define REF_GETRANDOMEVENT_H

#include "../src/common.h"
#include <vector>
#include <utility>

// Dataformats to store (pt[NPART], phi[NPART])
typedef std::pair<std::vector<float>, std::vector<float> > eventF_t;
typedef std::pair<std::vector<pt_t>, std::vector<phi_t> > eventHW_t;

std::vector<eventF_t> generateEvents(const int nEvent=NTEST, const int nPart=NPART);
std::vector<eventF_t> readEvents(const std::string fileName, const int nEvent=NTEST, const int nPart=NPART);
std::vector<eventF_t> readEventsFromHex(const std::string fileName, const int nEvent=NTEST, const int nPart=NPART);

eventHW_t to_HW(const eventF_t& in_event);

#endif
