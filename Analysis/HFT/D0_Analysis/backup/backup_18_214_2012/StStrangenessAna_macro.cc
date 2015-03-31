
void StStrangenessAna_macro(Int_t iBeamTimeNum, Int_t iAnalysis, Long64_t iStartEvent_SE, Long64_t iStopEvent_SE, Long64_t iStartEvent_ME, Long64_t iStopEvent_ME, TString file_ext, TString SE_in_list, TString ME_in_list)
{
    // root4star -b -q StStrangenessAna_macro.cc\(0,2,0,90000000,0,90000000,\"SE\",\"\",\"\"\)    here all files are used, automatic file list input
    // root4star -b -q StStrangenessAna_macro.cc\(4,1,0,10000,0,10000,\"SE\",\"XiOmega_SE_List_19GeV_A.txt\",\"\"\)  here one specific file is used

    // root4star -b -q StStrangenessAna_macro.cc\(5,8,0,10000,0,10000,\"SE\",\"JPsi_SE_List_27GeV_split_1_1\",\"\"\)  here one specific file is used

    // iAnalysis
    // 0 = phi meson
    // 1 = Xi-
    // 2 = Xi+
    // 3 = Omega-
    // 4 = Omega+
    // 5 = Lambda
    // 6 = antiLambda
    // 7 = K0S
    // 8 = JPsi/pi0


    // iBeamTimeNum:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200

    Int_t iTree = 0;
    if(iAnalysis == 0) iTree = 0; // phi
    if(iAnalysis == 1) iTree = 1; // Xi-
    if(iAnalysis == 2) iTree = 1; // Xi+
    if(iAnalysis == 3) iTree = 1; // Omega-
    if(iAnalysis == 4) iTree = 1; // Omega+
    if(iAnalysis == 5) iTree = 2; // Lambda
    if(iAnalysis == 6) iTree = 2; // anti-Lambda
    if(iAnalysis == 7) iTree = 3; // K0S
    if(iAnalysis == 8) iTree = 4; // jpsi

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

    cout << "iEliza = " << iEliza << endl;

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

    TString energy_name[7]       = {"7","11","39","62","19","27","200"};
    TString particle_name[5]     = {"Phi","XiOmega","Lambda","K0S","JPsi"};
    TString particle_ana_name[9] = {"Phi","XiM","XiP","OmegaM","OmegaP","Lambda","antiLambda","K0S","JPsi"};
    TString dir_name[5]          = {"phi","xiomega","lambda","k0s","jpsi"};
    TString disk_name[4]         = {"/eliza17/rnc/aschmah/AuAu","/eliza17/star/pwg/starlfs/aschmah/AuAu","/eliza14/star/pwg/starlfs/aschmah/AuAu","/eliza9/starprod/aschmah/AuAu"};
    //TString out_file_name      = disk_name[iEliza];
    //TString out_file_name        = "/global/scratch/sd/aschmah/AuAu";
    TString out_file_name        = "/eliza9/starprod/aschmah/BES/AuAu";
    out_file_name += energy_name[iBeamTimeNum];
    out_file_name += "/AnaOut/";
    out_file_name += dir_name[iTree];
    out_file_name += "/histo_out/";
    out_file_name += particle_ana_name[iAnalysis];
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

    astring = "./File_Lists/AuAu";
    astring += energy_name[iBeamTimeNum];
    astring += "/";
    if(SE_in_list != "")
    {
        astring += SE_in_list;
        StrangenessAna->setSEList(astring.Data());
    }
    else
    {
        if(ME_in_list == "")
        {
            astring += particle_name[iTree];
            astring += "_SE_List_";
            astring += energy_name[iBeamTimeNum];
            astring += "GeV.txt";
            StrangenessAna->setSEList(astring.Data());
        }
        else StrangenessAna->setSEList(" ");
    }

    astring = "./File_Lists/AuAu";
    astring += energy_name[iBeamTimeNum];
    astring += "/";
    if(ME_in_list != "")
    {
        astring += ME_in_list;
        StrangenessAna->setMEList(astring.Data());
    }
    else
    {
        if(SE_in_list == "")
        {
            astring += particle_name[iTree];
            astring += "_ME_List_";
            astring += energy_name[iBeamTimeNum];
            astring += "GeV.txt";
            StrangenessAna->setMEList(astring.Data());
        }
        else StrangenessAna->setMEList(" ");
    }

    //astring = disk_name[iEliza];
    //astring = "/global/scratch/sd/aschmah/AuAu";
    astring = "/eliza9/starprod/aschmah/BES/AuAu";
    astring += energy_name[iBeamTimeNum];
    astring += "/AnaOut/";
    astring += dir_name[iTree];
    astring += "/";
    StrangenessAna->setInputDir(astring.Data());

    astring = "../../Corrections/";
    astring += "Recentering_FitParams_AuAu";
    astring += energy_name[iBeamTimeNum];
    astring += ".root";
    StrangenessAna->setReCenteringFile(astring.Data());

    astring = "../../Corrections/";
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