
#ifndef __D0_Analysis_h__
#define __D0_Analysis_h__

#include <vector>
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "StMessMgr.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TNtuple.h"
#include "../../../Utils_run14/.sl64_gcc447/obj/StRoot/StD0Event/StD0Event.h"


//____________________________________________________________________________________________________
// Class to correct z-vertex dependence of refmult
class D0_Analysis {
  public:
    D0_Analysis();
    virtual ~D0_Analysis(); /// Default destructor

    // Functions
    void setInputDir(const TString inputdir);
    void setAnalysis(Int_t nAnalysis) {eAnalysis = nAnalysis;}
    void setOutputfile(const TString outputfile);
    void setStartEvent_SE(const Long64_t StartEvent_SE);
    void setStopEvent_SE(const Long64_t StopEvent_SE);
    void setStartEvent_ME(const Long64_t StartEvent_ME);
    void setStopEvent_ME(const Long64_t StopEvent_ME);
    void setSEList(const TString iSEList);
    void setMEList(const TString iMEList);
    void setBeamTimeNum(Int_t nBeamTimeNum) {eBeamTimeNum = nBeamTimeNum; cout << "Beam time number = " << eBeamTimeNum << endl;}
    void setEP_method(Int_t nEP_method) {eEP_method = nEP_method; cout << "Event plane method = " << eEP_method << endl;};
    void setReCenteringFile(TString inrecentfile) {recentfile = inrecentfile; cout << "Re-centering file = " << inrecentfile.Data() << endl;};
    void setShiftFile(TString inshiftfile) {shiftfile = inshiftfile; cout << "Shift method file = " << inshiftfile.Data() << endl;};
    void init();
    void loop();
    void finalize();

private:
    // Functions
    void clear() ; /// Clear all arrays

    // Data members
    Int_t eBeamTimeNum,eEP_method;
    Int_t eAnalysis;
    TString recentfile,shiftfile;
    TString pinputdir, SEList, MEList;
    TString poutputfile;
    TChain* input_SE;
    TChain* input_ME;
    Long64_t nentries;
    TFile* Outputfile;
    Long64_t nStartEvent_SE, nStopEvent_SE, nStartEvent_ME, nStopEvent_ME;
    StD0Event *D0_event;
    StD0Track *D0_track;
    Long64_t num_events_SE;
    Long64_t num_events_ME;

    TString HistName;

    ClassDef(D0_Analysis, 0)
};
#endif

