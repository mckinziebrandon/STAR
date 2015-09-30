
#include <TSystem>

void JetSimAnalysis_Macro(Long64_t StartEvent, Long64_t N_events, Double_t Jet_R, TString infilelist)
{
    // root4star -b -q JetSimAnalysis_Macro.cc\(0,10000,0.3,\"JetSim_list14_14\"\)

    cout << "StJetSimAnalysis_Macro started" << endl;

    gSystem ->Load("/global/homes/a/aschmah/local/fastjet_3_0_6/lib/libfastjet.so");
    gSystem ->Load("/global/homes/a/aschmah/local/fastjet_3_0_6/lib/libfastjettools.so");
    gSystem ->Load("/global/homes/a/aschmah/STAR/EP_Det_Analysis/.sl64_gcc447/lib/EP_Det_Sim_Event.so");
    gSystem ->Load("/global/homes/a/aschmah/STAR/EP_Det_Analysis/.sl64_gcc447/lib/EP_Det_Sim.so");
    gSystem ->Load("StJetSimAnalysis");

    //**************************** Set graphic style ***************************************
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,
                                                 greens, blues, NCont);
    gStyle->SetNumberContours(NCont);
    //**************************************************************************************

    cout << "Start of analysis" << endl;

    TString InListName = "File_lists/";
    InListName += infilelist;
    TString Outdir = "/project/projectdirs/star/pwg/starjetc/aschmah/Jet/Simulation/flow/out/";
    TString Indir  = "/project/projectdirs/star/pwg/starjetc/aschmah/Jet/Simulation/flow/";

    StJetSimAnalysis *StJetSimAnalysis_Ana = new StJetSimAnalysis();

    StJetSimAnalysis_Ana->setSEList(InListName.Data());
    StJetSimAnalysis_Ana->setListName(infilelist.Data());
    StJetSimAnalysis_Ana->setIndir(Indir);
    StJetSimAnalysis_Ana->setOutdir(Outdir);
    StJetSimAnalysis_Ana->setStartEvent(StartEvent);
    StJetSimAnalysis_Ana->setNEvents(N_events);
    StJetSimAnalysis_Ana->setJetR(Jet_R);
    StJetSimAnalysis_Ana->Init();
    StJetSimAnalysis_Ana->Make();
    StJetSimAnalysis_Ana->Finish();

    cout << "End of analysis" << endl;
}