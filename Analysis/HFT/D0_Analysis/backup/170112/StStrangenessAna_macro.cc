
void StStrangenessAna_macro(Int_t iBeamTimeNum, Int_t iAnalysis, Long64_t iStartEvent_SE, Long64_t iStopEvent_SE, Long64_t iStartEvent_ME, Long64_t iStopEvent_ME, TString file_ext)
{
    // iAnalysis
    // 0 = phi meson

    // root4star -b -q StStrangenessAna_macro.cc\(0,0,0,90000000,0,90000000,\"V3\"\)

    // iBeamTimeNum:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200

    cout << "StStrangenessAna started" << endl;
    cout << "Loading libraries ..." << endl;
    gSystem->Load("St_base");
    gSystem->Load("StUtilities");        // new addition 22jul99
    gSystem ->Load("../../../Utils/.sl53_gcc432/obj/StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.so");
    gSystem ->Load("../../../Utils/.sl53_gcc432/lib/StRefMultCorr.so");
    gSystem ->Load("./.sl53_gcc432/lib/StStrangenessAna.so");

    cout << "Constructor ..." << endl;
    StStrangenessAna* StrangenessAna = new StStrangenessAna();

    TString energy_name[7] = {"7","11","39","62","19","27","200"};
    TString out_file_name = "/eliza17/rnc/aschmah/AuAu";
    out_file_name += energy_name[iBeamTimeNum];
    out_file_name += "/AnaOut/phi/histo_out/Phi_v2_";
    out_file_name += energy_name[iBeamTimeNum];
    out_file_name += "GeV_histo_";
    out_file_name += file_ext;
    out_file_name += ".root";
    cout << "out_file_name = " << out_file_name.Data() << endl;

    TString astring;

    astring = "Beamenergy: ";
    astring += energy_name[iBeamTimeNum];
    astring += " GeV was selected";

    astring = "./File_Lists/Phi_SE_List_";
    astring += energy_name[iBeamTimeNum];
    astring += "GeV.txt";
    StrangenessAna->setSEList(astring.Data());

    astring = "./File_Lists/Phi_ME_List_";
    astring += energy_name[iBeamTimeNum];
    astring += "GeV.txt";
    StrangenessAna->setMEList(astring.Data());

    astring = "/eliza17/rnc/aschmah/AuAu";
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