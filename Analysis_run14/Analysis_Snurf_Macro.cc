#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


void Analysis_Snurf_Macro(Int_t StartEvent, Int_t NumOfEvents, Int_t fSE_ME_Flag, TString inputFile, TString suffix,
                          Int_t Analysis_num, Int_t BeamTime_num)
{

    // V1.0
    // SYSTEM:
    // 0 = PDSF
    // 1 = CARVER
    // 2 = PDSF/xrootd
    // 3 = CARVER/xrootd
    // 4 = CARVER/project
    // 5 = PDSF/project
    // 6 = RCF

    Int_t SYSTEM = 5;
    if(SYSTEM == 0) cout << "PDSF environment selected" << endl;
    if(SYSTEM == 1) cout << "CARVER environment selected" << endl;
    if(SYSTEM == 2) cout << "PDSF/xrootd environment selected" << endl;
    if(SYSTEM == 3) cout << "CARVER/xrootd environment selected" << endl;
    if(SYSTEM == 4) cout << "CARVER/project environment selected" << endl;
    if(SYSTEM == 5) cout << "PDSF/project environment selected" << endl;
    if(SYSTEM == 6) cout << "RCF environment selected" << endl;

    // Based on the LBNL pico DSTs, updated 10/19/11
    // Examples to run the analysis
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,9000000,0,\"Split_picoDst_7GeV_6501_6550\",\"V67\",112,0\)   // 7.7 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_11GeV_1751_1800\",\"V67\",1,1\)    // 11.5 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,10000,0,\"Split_picoDst_19GeV_3801_3850\",\"V67\",112,4\)    // 19.6 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_27GeV_8851_8900\",\"V67\",25,5\)    // 27 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,10000,0,\"Split_picoDst_39GeV_8801_8850\",\"V67\",112,2\)    // 39 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_62GeV_24401_24450\",\"V67\",201,3\)  // 62.4 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_UU193GeV_test\",\"V67\",9,3\)  // UU 193 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_xrootd_picoDst_200GeV_run11_104301_104350\",\"V67\",302,6\)  // 200 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_39GeV_19451_19500\",\"V67\",201,2\)  // 39 GeV project

    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_project_picoDst_200GeV_run11_24251_24300\",\"V67\",302,6\)

    // root4star -b -q Analysis_Snurf_Macro.cc\(0,100,0,\"Split_picoDst_14_5GeV_13176_13200\",\"V67\",9,7\)  // 14.5 GeV test

    // project
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_7GeV_13201_13400\",\"V67\",9,0\)   // 7.7 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_11GeV_3151_3200\",\"V67\",9,1\)    // 11.5 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_15GeV_23001_23050\",\"V67\",9,7\)  // 14.5 GeV test
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_19GeV_9001_9050\",\"V67\",9,4\)    // 19.6 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_27GeV_12801_12850\",\"V67\",9,5\)    // 27 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_39GeV_7901_7950\",\"V67\",9,2\)    // 39 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_62GeV_24551_24600\",\"V67\",302,3\)  // 62.4 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_project_picoDst_200GeV_run11_24251_24300\",\"V67\",9,6\) // 200 GeV

    // project pp200 test
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_pp200GeV_12001_12050\",\"V67\",9,6\) // pp200 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_pp510GeV_551_600\",\"V67\",9,6\) // pp510 GeV
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_200GeV_run10_29001_29050\",\"V67\",9,6\) // AuAu 200 GeV run 10 with BEMC


    // D0 analysis on RCF
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_run14_200GeV_9501_9550\",\"V67\",400,8\)
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,10,0,\"test\",\"V67\",400,8\)

    // D0 analysis on PDSF
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_200GeV_run14_preview2_15501_15550\",\"V67\",400,8\)  // preview2
    // root4star -b -q Analysis_Snurf_Macro.cc\(0,1000,0,\"Split_picoDst_200GeV_run14_100001_100050\",\"V67\",400,8\)  // final data

    // Analysis_num:
    // 1   = phi
    // 2   = Lambda
    // 102 = Sigma0 -> Lambda + e+ + e-
    // 3   = K0s
    // 4   = Xi
    // 5   = Omega
    // 6   = D0 + KStar
    // 7   = KStar +/-
    // 8   = hadron v2
    // 88  = hadron spectra
    // 9   = vertex + event display
    // 10  = phi correction
    // 11  = event plane analysis
    // 111 = event plane analysis with different harmonics + eta gap values
    // 112 = BBC analysis
    // 12  = phi same event K+K+ and K-K-
    // 13  = phi V0
    // 14  = anti-Lambda
    // 15  = rho
    // 16  = anti-Xi (Xi+)
    // 17  = Xi(1530) -> Xi- + pi+
    // 18  = Xi(1530) -> Xi+ + pi-
    // 19  = anti-Omega (Omega+)
    // 20  = Merged Omega + Xi analysis + anti-particles
    // 200 = Omega(2250) analysis
    // 201 = d4s analysis
    // 202 = Lambda event analysis
    // 21  = Sigma(1385) analysis
    // 22  = RefMult QA analysis
    // 23  = K->pi pi pi analysis
    // 24  = Theta+ -> K0S + p analysis
    // 25  = J/Psi -> e- + e+
    // 26  = Patrick's e+ e- analysis
    // 30  = proton-proton correlation analysis
    // 31  = fill phi histograms analysis
    // 124 = Patrick's single tree analysis
    // 125 = Event analysis
    // 130 = m2 vs nSigmaP analysis
    // 300 = Jet analysis
    // 301 = Jet analysis - with randomized tracks
    // 302 = Fill Jet tree
    // 35  = LambdaC+ analysis
    // 351 = LambdaC+ V0 analysis
    // 360 = Ach, mean pt analysis
    // 400 = D0 analysis with HFT

    // BeamTime_num:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200
    // 7 = 14.5
    // 8 = 200 run 14

    cout << "Analysis_Snurf_Macro started with " << NumOfEvents << " events" << endl;

    // LBNL
    gROOT->Macro("loadMuDst.C");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/obj/StRoot/StPicoDstMaker/StPicoDstMaker.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/obj/StRoot/StAlexEvent/StAlexEvent.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/obj/StRoot/StAlexV0Event/StAlexV0Event.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/obj/StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/obj/StRoot/StD0Event/StD0Event.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StPicoAlexEvent.so");
    //gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StRefMultCorr.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StStarBbcUtilities.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StParticleEvent.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/Std4sEvent.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StJetTrackEvent.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StLambdaEvent.so");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/StCombPID.so");
    //if(SYSTEM == 0 || SYSTEM == 2 || SYSTEM == 5) gSystem ->Load("/common/star/star53/packages/SL12c/.sl64_gcc447/lib/StBTofUtil.so");
    //if(SYSTEM == 0 || SYSTEM == 2 || SYSTEM == 5) gSystem ->Load("/common/star/star64/packages/SL14a/.sl64_gcc447/lib/StBTofUtil.so");
    //if(SYSTEM == 1 || SYSTEM == 3 || SYSTEM == 4) gSystem ->Load("StBTofUtil");
    gSystem ->Load("../Utils_run14/.sl64_gcc447/lib/MyEventEmbData.so");
    //gSystem ->Load("/global/homes/a/aschmah/local/fastjet_3_0_6/lib/libfastjet.so");
    //gSystem ->Load("/global/homes/a/aschmah/local/fastjet_3_0_6/lib/libfastjettools.so");
    gSystem ->Load("./.sl64_gcc447/lib/StV0TofCorrection.so");
    gSystem ->Load("./.sl64_gcc447/lib/Analysis_new_dev2.so");
    gSystem->Load("StRefMultCorr");
    //gSystem ->AddIncludePath("/global/homes/a/aschmah/local/include/");
    cout << "Shared libraries loaded" << endl;

    chain_input = new StChain();

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
    if(Analysis_num == 9) gStyle->SetCanvasPreferGL(kTRUE);
    //**************************************************************************************

    cout << "Start of analysis" << endl;

    Analysis_Snurf *Snurf_Ana = new Analysis_Snurf();

    TString out_suffix = "_out_";
    out_suffix+=suffix;
    out_suffix+=".root";
    Snurf_Ana->SetOutFileName(out_suffix.Data()); // suffix



    //--------------------------------------------------------------------
    cout << "Define re-centering and shift correction files" << endl;
    TString energy_name[9] = {"7","11","39","62","19","27","200","15","200run14"};
    TString astring;
    TString home_dir;
    if(SYSTEM == 0) home_dir = "/u";
    if(SYSTEM == 1 || SYSTEM == 3 || SYSTEM == 4 || SYSTEM == 5) home_dir = "/global/homes/a";
    if(SYSTEM == 2) home_dir = "/u";
    if(SYSTEM == 6) home_dir = "/star/u";
    astring = home_dir;
    astring += "/aschmah/STAR/Analysis/Corrections/";
    astring += "Recentering_FitParams_AuAu";
    astring += energy_name[BeamTime_num];
    astring += ".root";
    Snurf_Ana->SetReCenteringFile(astring.Data());

    astring = home_dir;
    astring += "/aschmah/STAR/Analysis/Corrections/";
    astring += "Shift_FitParams_AuAu";
    astring += energy_name[BeamTime_num];
    astring += ".root";
    Snurf_Ana->SetShiftFile(astring.Data());
    //--------------------------------------------------------------------



    if(BeamTime_num == 0) // 7.7 GeV
    {
        cout << "Au+Au @ 7.7 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu7/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu7/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu7/"); // location of file list
        }
    }
    if(BeamTime_num == 1) // 11.5 GeV
    {
        cout << "Au+Au @ 11.5 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu11/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu11/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/global/homes/b/bsmckinz/STAR/Analysis_run14/output/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu11/"); // location of file list
        }
    }
    if(BeamTime_num == 2) // 39 GeV
    {
        cout << "Au+Au @ 39 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu39/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu39/"); // location of file list
        }
        if(SYSTEM == 4) // CARVER/project
        {
            Snurf_Ana->SetOutputDir("/global/scratch/sd/aschmah/AuAu39/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu39/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu39/"); // location of file list
        }
        //if(SYSTEM == 5) // PDSF/project
        //{
        //    Snurf_Ana->SetOutputDir("/project/projectdirs/star/starprod/picodsts/Run10/AuAu/aschmah/39GeV/d4s/");
        //    Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu39/"); // location of file list
        //}
    }
    if(BeamTime_num == 3) // 62.4 GeV
    {
        cout << "Au+Au @ 62.4 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu62/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu62/"); // location of file list
        }
        if(SYSTEM == 1) // carver
        {
            Snurf_Ana->SetOutputDir("/global/scratch/sd/aschmah/AuAu62/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu62/"); // location of file list
        }
        if(SYSTEM == 2) // pdsf/xrootd
        {
            Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu62/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/xrootd/AuAu62/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu62/"); // location of file list
        }
    }
    if(BeamTime_num == 4) // 19.6 GeV
    {
        cout << "Au+Au @ 19.6 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu19/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu19/"); // location of file list
        }
        if(SYSTEM == 1) // carver
        {
            Snurf_Ana->SetOutputDir("/global/scratch/sd/aschmah/AuAu19/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu19/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu19/"); // location of file list
        }
    }
    if(BeamTime_num == 5) // 27 GeV
    {
        cout << "Au+Au @ 27 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            //Snurf_Ana->SetOutputDir("/eliza9/starprod/aschmah/AuAu27/AnaOut/");
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/starprod/picodsts/Run10/AuAu/aschmah/27GeV/d4s/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu27/"); // location of file list
        }
        if(SYSTEM == 1) // carver
        {
            Snurf_Ana->SetOutputDir("/global/scratch/sd/aschmah/AuAu27/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu27/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu27/"); // location of file list
        }
    }
    if(BeamTime_num == 6) // 200 GeV
    {
        cout << "Au+Au @ 200 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza17/star/pwg/starlfs/aschmah/AuAu200/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu200/"); // location of file list
        }
        if(SYSTEM == 2) // pdsf/xrootd
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/xrootd/AuAu200_run11/"); // location of file list
        }
        if(SYSTEM == 3) // carver/xrootd
        {
            Snurf_Ana->SetOutputDir("/global/scratch/sd/aschmah/AuAu200/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/xrootd/AuAu200_run11/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu200_run11/"); // location of file list
        }
    }
    if(BeamTime_num == 7) // 14.5 GeV
    {
        cout << "Au+Au @ 14.5 GeV was selected" << endl;
        if(SYSTEM == 0) // pdsf
        {
            Snurf_Ana->SetOutputDir("/eliza17/star/pwg/starlfs/aschmah/AuAu200/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu200/"); // location of file list
        }
        if(SYSTEM == 2) // pdsf/xrootd
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/xrootd/AuAu200_run11/"); // location of file list
        }
        if(SYSTEM == 3) // carver/xrootd
        {
            Snurf_Ana->SetOutputDir("/global/scratch/sd/aschmah/AuAu200/AnaOut/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/xrootd/AuAu200_run11/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/project/AuAu15/"); // location of file list
        }
    }
    if(BeamTime_num == 8) // 200 GeV run 14
    {
        cout << "Au+Au @ 200 GeV run14 was selected" << endl;
        if(SYSTEM == 6) // RCF
        {
            Snurf_Ana->SetOutputDir("/star/institutions/lbl/aschmah/AuAu200_run14/");
            Snurf_Ana->SetInDirList("/star/u/aschmah/STAR/File_lists/AuAu200_run14/"); // location of file list
        }
        if(SYSTEM == 5) // PDSF/project
        {
            Snurf_Ana->SetOutputDir("/project/projectdirs/star/rnc/aschmah/brandon/");
            Snurf_Ana->SetInDirList("/global/homes/a/aschmah/STAR/File_lists/AuAu200_run14/"); // location of file list
        }
    }


    time_t start,end;
    double dif;

    Snurf_Ana->SetInputFile(inputFile);
    Snurf_Ana->SetStartEvent(StartEvent);
    Snurf_Ana->SetNEvents(NumOfEvents);
    Snurf_Ana->SetSE_ME(fSE_ME_Flag);
    Snurf_Ana->SetAnalysisNum(Analysis_num);
    Snurf_Ana->SetBeamTimeNum(BeamTime_num);
    Snurf_Ana->Init();
    time (&start);
    Snurf_Ana->Make();
    time (&end);
    dif = difftime (end,start);
    cout << "Time for make loop = " << dif << " seconds" << endl;
    Snurf_Ana->Finish();



    cout << "End of analysis" << endl;



}


