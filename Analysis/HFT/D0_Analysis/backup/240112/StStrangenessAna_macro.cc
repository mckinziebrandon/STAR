
void StStrangenessAna_macro(Int_t iBeamTimeNum, Int_t iAnalysis, Long64_t iStartEvent_SE, Long64_t iStopEvent_SE, Long64_t iStartEvent_ME, Long64_t iStopEvent_ME, TString file_ext)
{
    // iAnalysis
    // 0 = phi meson
    // 1 = Xi-
    // 2 = Xi+
    // 3 = Omega-
    // 4 = Omega+
    // 5 = Lambda
    // 6 = K0S

    // root4star -b -q StStrangenessAna_macro.cc\(0,0,0,90000000,0,90000000,\"V9\"\)

    // iBeamTimeNum:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200

    Int_t iTree = 0;
    if(iAnalysis == 0) iTree = 0;
    if(iAnalysis == 1) iTree = 1;
    if(iAnalysis == 2) iTree = 1;
    if(iAnalysis == 3) iTree = 1;
    if(iAnalysis == 4) iTree = 1;
    if(iAnalysis == 5) iTree = 2;
    if(iAnalysis == 6) iTree = 3;

    Int_t iEliza = 0;
    if(iBeamTimeNum == 0) // 7.7 GeV
    {
        if(iAnalysis == 0) iEliza = 0; // phi meson
        if(iAnalysis >= 1 && iAnalysis <= 4) iEliza = 0; // Xi+/-,Omega+/-
        if(iAnalysis == 5) iEliza = 1; // Lambda
        if(iAnalysis == 6) iEliza = 1; // K0S
    }
    if(iBeamTimeNum == 1) // 11.5 GeV
    {
        if(iAnalysis == 0) iEliza = 0; // phi meson
        if(iAnalysis >= 1 && iAnalysis <= 4) iEliza = 0; // Xi+/-,Omega+/-
        if(iAnalysis == 5) iEliza = 1; // Lambda
        if(iAnalysis == 6) iEliza = 1; // K0S
    }
    if(iBeamTimeNum == 2) // 39 GeV
    {
        if(iAnalysis == 0) iEliza = 3; // phi meson
        if(iAnalysis >= 1 && iAnalysis <= 4) iEliza = 0; // Xi+/-,Omega+/-
        if(iAnalysis == 5) iEliza = 2; // Lambda
        if(iAnalysis == 6) iEliza = 0; // K0S
    }
    if(iBeamTimeNum == 3) // 62.4 GeV
    {
        if(iAnalysis == 0) iEliza = 3; // phi meson
        if(iAnalysis >= 1 && iAnalysis <= 4) iEliza = 3; // Xi+/-,Omega+/-
        if(iAnalysis == 5) iEliza = 3; // Lambda
        if(iAnalysis == 6) iEliza = 3; // K0S
    }
    if(iBeamTimeNum == 4) // 19.6 GeV
    {
        if(iAnalysis == 0) iEliza = 0; // phi meson
        if(iAnalysis >= 1 && iAnalysis <= 4) iEliza = 0; // Xi+/-,Omega+/-
        if(iAnalysis == 5) iEliza = 1; // Lambda
        if(iAnalysis == 6) iEliza = 0; // K0S
    }
    if(iBeamTimeNum == 5) // 27 GeV
    {
        if(iAnalysis == 0) iEliza = 3; // phi meson
        if(iAnalysis >= 1 && iAnalysis <= 4) iEliza = 3; // Xi+/-,Omega+/-
        if(iAnalysis == 5) iEliza = 3; // Lambda
        if(iAnalysis == 6) iEliza = 3; // K0S
    }

    cout << "StStrangenessAna started" << endl;
    cout << "Loading libraries ..." << endl;
    gSystem->Load("St_base");
    gSystem->Load("StUtilities");        // new addition 22jul99
    gSystem ->Load("../../../Utils/.sl53_gcc432/obj/StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.so");
    gSystem ->Load("../../../Utils/.sl53_gcc432/obj/StRoot/StAlexV0Event/StAlexV0Event.so");
    gSystem ->Load("../../../Utils/.sl53_gcc432/lib/StRefMultCorr.so");
    gSystem ->Load("./.sl53_gcc432/lib/StStrangenessAna.so");

    cout << "Constructor ..." << endl;
    StStrangenessAna* StrangenessAna = new StStrangenessAna();

    TString energy_name[7]   = {"7","11","39","62","19","27","200"};
    TString particle_name[4] = {"Phi","XiOmega","Lambda","K0S"};
    TString dir_name[4]      = {"phi","xiomega","lambda","k0s"};
    TString disk_name[4]     = {"/eliza17/rnc/aschmah/AuAu","/eliza17/star/pwg/starlfs/aschmah/AuAu","/eliza14/star/pwg/starlfs/aschmah/AuAu","/eliza9/starprod/aschmah/AuAu"};
    TString out_file_name    = disk_name[iEliza];;
    out_file_name += energy_name[iBeamTimeNum];
    out_file_name += "/AnaOut/";
    out_file_name += dir_name[iTree];
    out_file_name += "/histo_out/";
    out_file_name += particle_name[iTree];
    out_file_name += "_v2_";
    out_file_name += energy_name[iBeamTimeNum];
    out_file_name += "GeV_histo_";
    out_file_name += file_ext;
    out_file_name += ".root";
    cout << "out_file_name = " << out_file_name.Data() << endl;

    TString astring;

    astring = "Beamenergy: ";
    astring += energy_name[iBeamTimeNum];
    astring += " GeV was selected";

    astring = "./File_Lists/";
    astring += particle_name[iTree];
    astring += "_SE_List_";
    astring += energy_name[iBeamTimeNum];
    astring += "GeV.txt";
    StrangenessAna->setSEList(astring.Data());

    astring = "./File_Lists/";
    astring += particle_name[iTree];
    astring += "_ME_List_";
    astring += energy_name[iBeamTimeNum];
    astring += "GeV.txt";
    StrangenessAna->setMEList(astring.Data());

    astrubg = disk_name[iEliza];
    astring += energy_name[iBeamTimeNum];
    astring += "/AnaOut/phi/";
    StrangenessAna->setInputDir(astring.Data());

    astring = "/u/aschmah/STAR/Analysis/Corrections/";
    astring += "Recentering_FitParams_AuAu";
    astring += energy_name[iBeamTimeNum];
    astring += ".root";
    StrangenessAna->setReCenteringFile(astring.Data());

    astring = "/u/aschmah/STAR/Analysis/Corrections/";
    astring += "Shift_FitParams_AuAu";
    astring += energy_name[iBeamTimeNum];
    astring += ".root";
    StrangenessAna->setShiftFile(astring.Data());

    StrangenessAna->setAnalysis(iAnalysis);
    StrangenessAna->setOutputfile(out_file_name.Data());
    StrangenessAna->setStopEvent_SE(iStopEvent_SE);
    StrangenessAna->setStartEvent_SE(iStartEvent_SE);
    StrangenessAna->setStopEvent_ME(iStopEvent_ME);
    StrangenessAna->setStartEvent_ME(iStartEvent_ME);
    StrangenessAna->setBeamTimeNum(iBeamTimeNum);
    StrangenessAna->init();
    StrangenessAna->loop();
    StrangenessAna->finalize();

}