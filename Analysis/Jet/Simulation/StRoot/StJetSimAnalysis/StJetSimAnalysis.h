#ifndef StJetSimAnalysis_hh
#define StJetSimAnalysis_hh

#include "StarClassLibrary/SystemOfUnits.h"
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
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
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "/global/homes/a/aschmah/STAR/EP_Det_Analysis/StRoot/EP_Det_Sim_Event/EP_Det_Sim_Event.h"

#include "fastjet/config.h"             // will allow a test for FJ3
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

#include "StPhysicalHelixD.hh"

using namespace std;
using namespace fastjet;

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

static const Double_t Pi = TMath::Pi();

class StJetSimAnalysis {
public:
    StJetSimAnalysis();
    virtual ~StJetSimAnalysis();

    void  setSEList(TString iSEList) {SEList = iSEList;};
    void  setListName(TString iListName) {ListName = iListName;}
    void  setIndir(TString iIndir) {eIndir = iIndir;};
    void  setOutdir(TString iOutdir) {eOutdir = iOutdir;};
    void  setStartEvent(Long64_t iStartEvent) {eStartEvent = iStartEvent;}
    void  setNEvents(Long64_t iN_Events) {N_Events = iN_Events;}
    void  setJetR(Double_t iJet_R) {Jet_R = iJet_R;}

    void Init();
    void Make();
    void Finish();

private:

    TH1D* h_jet_pt_sub, *h_Delta_phi, *h_Delta_phi_trigger;
    TH2D* h_jet_pt_vs_Delta_phi, *h_jet_pt_vs_trigger_pt, *h_trigger_pt_vs_delta_phi;
    TH2D* h_track_pt_vs_delta_phi, *h_jet_recoil_pt_vs_Delta_phi;
    Long64_t N_Events, eStartEvent;
    TString SEList, eOutdir, eIndir, ListName;
    TChain* input_chain;
    Double_t Jet_R;
    TFile* Outputfile;
    ClassDef(StJetSimAnalysis,1)

};

#endif
