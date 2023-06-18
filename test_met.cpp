/*
MET calculation from PF objects
*/
#include "src/common.h"
#include "src/met.h"
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
  cout << "Generated " << eventsRef.size() << " events" << endl;

  std::ofstream fout("met.csv");
  fout << "MET_Ref,MET_HW,HW-REf\n";

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

    float met2Ref, metphiRef;
    pt2_t met2HW;
    phi_t metphiHW;
    met_ref(&ptsRef[0], &phisRef[0], met2Ref, metphiRef);
    met_hw(&ptsHW[0], &phisHW[0], met2HW, metphiHW);

    std::cout << "           >>>>> " << met2Ref << ' ' << met2HW << std::endl;
    break;
    const double metRef = std::sqrt(met2Ref);
    const double metHW = std::sqrt(float(met2HW));
    const double dmet = metRef-metHW;
    metDiff += dmet;
    metDiff2 += dmet*dmet;
    if ( std::abs(dmet) > std::abs(metMaxRef-metMaxHW) ) {
      metMaxRef = metRef;
      metMaxHW = metHW;
    }
    else if ( std::abs(dmet) < std::abs(metMinRef-metMinHW) ) { 
      metMinRef = metRef;
      metMinHW = metHW;
    }
#if(DEBUG==1)
    cout << "Ref=" << metRef << " metHW=" << metHW << " Difference=" << (metHW-metRef) << endl;
#endif
    fout << metRef << "," << metHW << "," << (metHW-metRef) << endl;
  }
  cout << endl;

  cout << "MET with maximum difference:\n"
       << " Ref=" << metMaxRef << " HW=" << metMaxHW << " difference=" << (metMaxHW-metMaxRef) << "\n"
       << "MET with minimum difference:\n"
       << " Ref=" << metMinRef << " HW=" << metMinHW << " difference=" << (metMinHW-metMinRef) << "\n"
       << "Difference between HW and Reference:\n"
       << " <HW-Ref>=" << (metDiff/eventsRef.size())
       << " RMS=" << (metDiff2-metDiff*metDiff/eventsRef.size())/eventsRef.size() << endl;

  return 0;
}

