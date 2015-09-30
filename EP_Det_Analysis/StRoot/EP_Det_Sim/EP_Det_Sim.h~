#ifndef EP_Det_Sim_hh
#define EP_Det_Sim_hh

#include "StarClassLibrary/SystemOfUnits.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TTree.h"
#include "StMessMgr.h"
//#include <iostream.h>
#include "TChain.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TRandom.h"
#include "TRandom3.h"

#include "./StRoot/EP_Det_Sim_Event/EP_Det_Sim_Event.h"


class EP_Det_Sim {
public:
    EP_Det_Sim();
    virtual ~EP_Det_Sim();

    Int_t Init();
    Int_t Make();
    Int_t Finish();
    void  SetNEvents(Long64_t nevents) {nEvents = nevents; cout << "Number of analyzed events set to " << nEvents << endl;}
    void  SetBeamtime(Int_t nBeamtime) {eBeamtime = nBeamtime; cout << "Beamtime " << eBeamtime << endl;}
    void  SetOutFileName(TString outfileName) {output_file_name = outfileName; cout << "Output filename = " << output_file_name << endl;}
    void  SetEtaRange(Double_t eta_low_in, Double_t eta_high_in) {eta_low = eta_low_in; eta_high = eta_high_in; cout << "eta_low: " << eta_low
        << "eta_high: " << eta_high << endl;}

private:

    Int_t eBeamtime;
    Long64_t nEvents;
    TString output_file_name;
    TFile* dNdeta_file;
    Double_t eta_low, eta_high;
    ClassDef(EP_Det_Sim,1)

};

#endif
