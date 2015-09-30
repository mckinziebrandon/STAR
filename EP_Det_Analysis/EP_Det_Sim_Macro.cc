#include <TSystem>

void EP_Det_Sim_Macro(Long64_t RunEvents, TString out_file_name, Int_t Beamtime, Double_t eta_low, Double_t eta_high)
{
    // Generates simulated flow trees
    // root4star -b -q EP_Det_Sim_Macro.cc\(10000,\"/global/scratch/sd/aschmah/Flow_Sim/EP_Det_Sim_out_200GeV_V2.root\",1\)

    // root4star -b -q EP_Det_Sim_Macro.cc\(10000,\"/project/projectdirs/star/pwg/starjetc/aschmah/Jet/Simulation/flow/Jet_sim_flow_200GeV.root\",1,-1.0,1.0\)
    // root4star -b -q EP_Det_Sim_Macro.cc\(10000,\"/project/projectdirs/star/rnc/aschmah/brandon/Jet_sim_flow_200GeV.root\",1,-1.0,1.0\)

    // Beamtime
    // 0 = 19.6 GeV
    // 1 = 200 GeV

    cout << "EP_Det_Sim_Macro started with " << RunEvents << " events" << endl;

    gSystem ->Load("./.sl64_gcc447/lib/EP_Det_Sim_Event.so");
    gSystem ->Load("./.sl64_gcc447/lib/EP_Det_Sim.so");

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

    EP_Det_Sim *EP_Det_Sim_Ana = new EP_Det_Sim();

    EP_Det_Sim_Ana->SetOutFileName(out_file_name.Data()); // suffix
    EP_Det_Sim_Ana->SetNEvents(RunEvents);
    EP_Det_Sim_Ana->SetBeamtime(Beamtime);
    EP_Det_Sim_Ana->SetEtaRange(eta_low,eta_high);
    EP_Det_Sim_Ana->Init();
    EP_Det_Sim_Ana->Make();
    EP_Det_Sim_Ana->Finish();

    cout << "End of analysis" << endl;
}
