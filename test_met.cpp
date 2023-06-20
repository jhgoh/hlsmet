/*
MET calculation from PF objects
*/
#include "src/common.h"
#include "src/met.h"
#include "src/jercorr.h"
#include "bench/readEvents.h"
#include "bench/met_ref.h"

#include <fstream>

int alg_test();
int main()
{
  alg_test();
  return 0;
}

int alg_test()
{
  using namespace std;
  //auto eventsRef = readEvents("/users/junwonoh/out_TTbar.dump");
  auto eventsRef = readEventsFromHex("/users/jhgoh/work/FPGA/hlsmet/TTbar_1000evt_54part_v2.dump");
  //auto eventsRef = generateEvents(NTEST, NPART);
  cout << "Produced " << eventsRef.size() << " events" << endl;

  std::ofstream fout("met.csv");
  fout << "MET_Ref,MET_HW,CORRMET_HW\n";

  double metMaxRef = 0, metMaxHW = 0;
  double metMinRef = 0, metMinHW = 1e9;
  double metDiff = 0, metDiff2 = 0;
  for ( size_t i=0; i<eventsRef.size(); ++i ) {
    cout << "Event " << (i+1) << "/" << eventsRef.size();
#if(DEBUG==0)
    cout << "\r";
#else
    cout << "\n";
#endif

    const auto eventRef = eventsRef[i];
    const auto eventHW = to_HW(eventRef);

    auto ptsRef = eventRef.first;
    auto phisRef = eventRef.second;
    auto ptsHW = eventHW.first;
    auto phisHW = eventHW.second;

    float metRef, metphiRef;
    pt_t metHW;
    phi_t metphiHW;
    met_ref(&ptsRef[0], &phisRef[0], metRef, metphiRef);
    met_hw(&ptsHW[0], &phisHW[0], metHW, metphiHW);

    pt_t corrmetHW;
    phi_t corrmetphiHW;
    // FIXME: we need jet pts, phis and etas and pass them into the following function.
    jercorrmet_hw(metHW, metphiHW, &ptsHW[0], &phisHW[0], &phisHW[0], corrmetHW, corrmetphiHW);

    const float dmet = metRef-float(metHW);
    metDiff += dmet;
    metDiff2 += dmet*dmet;
    if ( std::abs(dmet) > std::abs(metMaxRef-float(metMaxHW)) ) {
      metMaxRef = metRef;
      metMaxHW = metHW;
    }
    else if ( std::abs(dmet) < std::abs(metMinRef-float(metMinHW)) ) { 
      metMinRef = metRef;
      metMinHW = metHW;
    }
#if(DEBUG==1)
    cout << "Ref=" << metRef << " metHW=" << metHW << " Difference=" << (metHW-metRef) << endl;
#endif
    fout << metRef << "," << metHW << "," << corrmetHW << endl;
  }
  cout << endl;

  cout << "MET with maximum difference:\n"
       << " Ref=" << metMaxRef << " HW=" << metMaxHW << " difference=" << (metMaxHW-metMaxRef) << "\n"
       << "MET with minimum difference:\n"
       << " Ref=" << metMinRef << " HW=" << metMinHW << " difference=" << (metMinHW-metMinRef) << "\n"
       << "Difference between HW and Reference:\n"
       << " <HW-Ref>=" << (metDiff/eventsRef.size())
       << " RMS=" << std::sqrt((metDiff2-metDiff*metDiff/eventsRef.size())/eventsRef.size()) << endl;

  return 0;
}

