
//#include <assert.h>
#include <fstream>
#include <string>
#include "StMessMgr.h"
#include "StStrangenessAna.h"
#include "TString.h"
#include "TLatex.h"
#include "TPad.h"
#include "../event_plane/StRoot/event_plane_analysis/event_plane_analysis_Func.h"

ClassImp(StStrangenessAna)

using std::ifstream ;
using std::string ;
using std::vector ;

static const Int_t N_Analysis = 9;

static const TString V0_TREE[N_Analysis]   = {"AlexPhiMesonEvents","XiM_v2_tree","XiP_v2_tree","OmegaM_v2_tree","OmegaP_v2_tree","Lambda_X_NT","Lambda_X_NT","Lambda_X_NT","AlexPhiMesonEvents"};
static const TString V0_BRANCH[N_Analysis] = {"Events","XiM_v2_branch","XiP_v2_branch","OmegaM_v2_branch","OmegaP_v2_branch","Lambda_X_NT","Lambda_X_NT","Lambda_X_NT","Events"};

static char* ALEXV0_EVENT_TREE;
static char* ALEXV0_EVENT_BRANCH;

//static char* ALEXV0_EVENT_TREE       = V0_TREE[1].Data();
//static char* ALEXV0_EVENT_BRANCH     = V0_BRANCH[1].Data();

static Int_t SE_input_flag = 1;
static Int_t ME_input_flag = 1;

//------------------------------------------------------------------------------------------------------------------------------------
static const Int_t N_Beamtime = 7;
static const Int_t n_bad_run_numbers[N_Beamtime] = {328,27,38,105,35,34,1};
static const Int_t bad_run_list_7GeV[328]  = {11114084,11114085,11114086,11114088,11114089,11114094,11114095,11114100,11114109,11115005,11115007,11115013,11115019,11115025,11115027,11115028,11115030,11115032,11115051,11115062,11115064,11115069,11115072,11115078,11115079,11115080,11115086,11115088,11115094,11116001,11116002,11116005,11116006,11116010,11116014,11116020,11116023,11116028,11116060,11116061,11116062,11116064,11116068,11116070,11116072,11116073,11116075,11117002,11117006,11117031,11117033,11117034,11117036,11117039,11117044,11117045,11117046,11117052,11117055,11117063,11117064,11117071,11117075,11117085,11117088,11117089,11117090,11117093,11117094,11117095,11117098,11117100,11117103,11117104,11117107,11118007,11118008,11118016,11118024,11118025,11118026,11118039,11118044,11119001,11119003,11119006,11119007,11119009,11119012,11119013,11119015,11119016,11119017,11119022,11119024,11119026,11119029,11119030,11119056,11119057,11119060,11119062,11119067,11119069,11119070,11119071,11119074,11119075,11119077,11119079,11119081,11119090,11119091,11119100,11119101,11120003,11120006,11120008,11120011,11120014,11120019,11120023,11120030,11120034,11120037,11120039,11120040,11120045,11120052,11120057,11120062,11120063,11120069,11120070,11120071,11120074,11120077,11120078,11120084,11120092,11121006,11121014,11121015,11121019,11121029,11121030,11121034,11121035,11121043,11121044,11121054,11121058,11121066,11121067,11121070,11121075,11121082,11122001,11122007,11122008,11122010,11122017,11122024,11122037,11122038,11122047,11122048,11122049,11122050,11122053,11122058,11122062,11122069,11122073,11122078,11122085,11122097,11123003,11123004,11123015,11123026,11123028,11123040,11123044,11123055,11123057,11123058,11123059,11123067,11123075,11123076,11123077,11123079,11123081,11123084,11123086,11123088,11123089,11123093,11123094,11123095,11123100,11123101,11123102,11123104,11124001,11124005,11124007,11124008,11124015,11124016,11124018,11124041,11124046,11124050,11124051,11124052,11124053,11124058,11124060,11124061,11124062,11124063,11124064,11124065,11124066,11124069,11125002,11125003,11125004,11125005,11125006,11125008,11125012,11125013,11125014,11125015,11125016,11125017,11125020,11125021,11125022,11125023,11125073,11125081,11125089,11125090,11125096,11125097,11126005,11126006,11126007,11126016,11126018,11126022,11126023,11127001,11127002,11127043,11128005,11128012,11128018,11128050,11128056,11128072,11129018,11129022,11129028,11129051,11130027,11130034,11130057,11131038,11131062,11132013,11132070,11133006,11133019,11134053,11134060,11134067,11134076,11135068,11136003,11136005,11136006,11136007,11136008,11136012,11136013,11136014,11136061,11136076,11136101,11136130,11136160,11136163,11137019,11138027,11138049,11138086,11138124,11139014,11140076,11140086,11141063,11142117,11143026,11143028,11144001,11144009,11144031,11144033,11144040,11144043,11144052,11145008,11145028,11145035,11146061,11146076,11146079,11147004,11147006,11147014,11147017,11147021,11147023};
static const Int_t bad_run_list_11GeV[27]  = {11148039,11148045,11149001,11149008,11149010,11149011,11149015,11149047,11150016,11150025,11150028,11151036,11151040,11151050,11152016,11152036,11152078,11153032,11153042,11155001,11155009,11156003,11156009,11157012,11158006,11158022,11158024};
static const Int_t bad_run_list_19GeV[35]  = {12113091,12114007,12114035,12114078,12114092,12114116,12115009,12115014,12115015,12115016,12115018,12115019,12115020,12115022,12115023,12115062,12115073,12115093,12115094,12116012,12116054,12117010,12117016,12117020,12117065,12119040,12119042,12120017,12120026,12121017,12121022,12121034,12121050,12121067,12122019};
static const Int_t bad_run_list_27GeV[34]  = {12172050,12172051,12172055,12173030,12173031,12173032,12173033,12173034,12174067,12174085,12175062,12175087,12175113,12175114,12175115,12176001,12176044,12176054,12176071,12177015,12177061,12177092,12177099,12177101,12177106,12177107,12177108,12178003,12178004,12178005,12178006,12178013,12178099,12178120};
static const Int_t bad_run_list_39GeV[38]  = {11199124,11100002,11100045,11101046,11102012,11102051,11102052,11102053,11102054,11102055,11102058,11103035,11103056,11103058,11103092,11103093,11105052,11105053,11105054,11105055,11107007,11107042,11107057,11107061,11107065,11107074,11108101,11109013,11109077,11109088,11109090,11109127,11110013,11110034,11110073,11110076,11111084,11111085};
static const Int_t bad_run_list_62GeV[105] = {11080072,11081023,11081025,11082012,11082013,11082046,11082056,11082057,11084009,11084011,11084012,11084013,11084020,11084021,11084035,11084044,11084064,11085015,11085025,11085030,11085046,11085055,11085056,11085057,11086005,11086007,11087001,11087002,11087003,11087004,11088013,11089026,11089028,11089029,11089055,11089068,11089072,11091007,11091015,11091021,11091078,11092010,11092011,11092012,11092032,11092033,11092034,11092067,11092096,11093001,11094016,11094017,11094018,11094019,11094020,11094021,11094022,11094023,11094024,11094027,11094028,11094042,11094044,11094045,11094046,11094047,11094048,11094050,11094051,11094052,11094053,11094054,11094055,11094074,11094075,11094077,11095001,11095002,11095003,11095004,11095005,11095006,11095009,11095010,11095011,11095012,11095013,11095014,11095015,11095022,11095040,11095048,11095050,11095051,11095061,11095062,11095063,11095064,11095082,11095087,11096024,11096039,11096043,11096044,11097093};
static const Int_t bad_run_list_200GeV[1]  = {1};
//------------------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------------------
static const Int_t Mass_electron = 0.00051099892;
//------------------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------------------
// EP method -> resolution
// 0 = full TPC event plane, phi-weights
// 1 = eta sub event plane, phi-weights
// 2 = full TPC event plane, re-centering
// 3 = eta sub event plane, re-centering

// Cuts for event plane calculation
static const Float_t nHitsFitA_EP_cut      = 14;
static const Float_t nHitsPossA_EP_cut     = 0;
static const Float_t nHits_ratio_EP_cut    = 0.52;
static const Float_t MomentumA_EP_low_cut  = 0.15;
static const Float_t MomentumA_EP_high_cut = 5.0;
static const Float_t dcaAB_EP_cut          = 2.0;
static const Float_t eta_EP_cut            = 1.0;
static const Float_t Event_radius_cut      = 2.0;
static const Float_t eta_gap               = 0.05; // 0.05

// 70-80, 60-70, 50-60, 40-50, 30-40, 20-30, 10-20, 5-10, 0-5

//Event plane resolution for system: Au+Au @ 7.7 GeV
static const Double_t Resolution_EP0_7GeV[9] = {14.9913,18.2929,25.6835,35.1828,44.0758,47.8632,43.9594,31.9779,15.5305};
static const Double_t Resolution_EP1_7GeV[9] = {9.98479,11.4681,17.6196,24.2844,31.4525,34.8621,32.922,25.0144,15.4289};
static const Double_t Resolution_EP2_7GeV[9] = {15.0106,18.487,26.2469,36.0243,44.9701,48.856,45.1212,32.753,15.8519};
static const Double_t Resolution_EP3_7GeV[9] = {9.69313,11.5798,17.8237,24.8746,32.0848,35.6442,33.6759,25.871,15.8069};

//Event plane resolution for system: Au+Au @ 11.5 GeV
static const Double_t Resolution_EP0_11GeV[9] = {17.4605,22.3764,31.6758,42.0158,51.0439,54.9112,51.9088,41.1008,23.9172};
static const Double_t Resolution_EP1_11GeV[9] = {10.9215,14.2854,21.1834,29.3017,36.7378,40.4501,38.4486,30.8478,21.1104};
static const Double_t Resolution_EP2_11GeV[9] = {18.051,22.927,32.3827,43.045,52.2768,56.2791,53.551,43.1764,27.7126};
static const Double_t Resolution_EP3_11GeV[9] = {11.1566,14.554,21.5474,29.9147,37.6087,41.4212,39.5566,32.1547,22.8879};

//Event plane resolution for system: Au+Au @ 19.6 GeV
static const Double_t Resolution_EP0_19GeV[9] = {21.098,27.6305,37.6858,48.9311,57.864,62.1655,59.2439,48.58,33.8356};
static const Double_t Resolution_EP1_19GeV[9] = {12.3314,17.4169,25.2445,34.2429,42.0549,46.1318,44.1126,36.2036,26.4456};
static const Double_t Resolution_EP2_19GeV[9] = {21.2273,27.986,38.2177,49.7074,58.7821,63.195,60.2781,49.5695,34.2717};
static const Double_t Resolution_EP3_19GeV[9] = {12.423,17.584,25.59,34.8278,42.79,47.0165,44.9983,36.9999,27.005};

//Event plane resolution for system: Au+Au @ 27 GeV
static const Double_t Resolution_EP0_27GeV[9] = {22.6523,30.4555,41.3613,52.8548,61.5597,65.361,62.2394,51.3692,36.9666};
static const Double_t Resolution_EP1_27GeV[9] = {13.3048,19.1344,27.8721,37.3431,45.0804,48.9593,46.6782,38.3997,28.5737};
static const Double_t Resolution_EP2_27GeV[9] = {22.9087,30.9228,42.0945,53.8355,62.6966,66.5373,63.4645,52.5435,37.6455};
static const Double_t Resolution_EP3_27GeV[9] = {13.4224,19.4013,28.3912,38.0913,46.0539,50.0186,47.7698,39.3833,29.2801};

//Event plane resolution for system: Au+Au @ 39 GeV
static const Double_t Resolution_EP0_39GeV[9] = {25.7775,33.7427,45.13,56.8738,65.384,69.2147,66.0776,55.374,41.4523};
static const Double_t Resolution_EP1_39GeV[9] = {15.0429,21.3115,30.5863,40.5551,48.5286,52.59,50.197,41.6513,31.8176};
static const Double_t Resolution_EP2_39GeV[9] = {25.9956,34.1819,45.8817,57.8561,66.4951,70.363,67.2629,56.386,41.5277};
static const Double_t Resolution_EP3_39GeV[9] = {15.1681,21.5897,31.1031,41.3207,49.5047,53.6854,51.3163,42.5969,32.1888};

//Event plane resolution for system: Au+Au @ 62.4 GeV
static const Double_t Resolution_EP0_62GeV[9] = {30.3177,37.7198,48.7711,60.3508,68.5697,72.2626,69.4578,59.8213,47.0744};
static const Double_t Resolution_EP1_62GeV[9] = {15.0469,22.4996,32.1578,42.6574,50.6306,54.4605,51.6422,42.2871,30.3959};
static const Double_t Resolution_EP2_62GeV[9] = {28.5739,37.8401,50.0268,62.0979,70.3985,74.0446,71.0833,60.6593,45.9837};
static const Double_t Resolution_EP3_62GeV[9] = {16.3416,23.8706,33.878,44.757,53.0614,57.2351,54.8345,46.0863,35.4889};
//------------------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------------------
static const Double_t z_acceptance[N_Beamtime]   = {70.0,50.0,40.0,40.0,70.0,70.0,40.0};
//------------------------------------------------------------------------------------------------------------------------------------

// nBeamTimeNum:
// 0 = 7.7
// 1 = 11.5
// 2 = 39
// 3 = 62.4
// 4 = 19.6
// 5 = 27
// 6 = 200

//____________________________________________________________________________________________________
// Default constructor
StStrangenessAna::StStrangenessAna()
{
  clear() ;
}

//____________________________________________________________________________________________________
// Default destructor
StStrangenessAna::~StStrangenessAna()
{
}

//____________________________________________________________________________________________________
void StStrangenessAna::clear()
{

}

void StStrangenessAna::setInputDir(const TString inputdir)
{
    pinputdir = inputdir.Copy();
    cout << "Input directory was set to: " << pinputdir.Data() << endl;
}
void StStrangenessAna::setOutputfile(const TString outputfile)
{
    poutputfile = outputfile.Copy();
    cout << "Output file was set to: " << poutputfile.Data() << endl;
}
void StStrangenessAna::setSEList(const TString iSEList)
{
    SEList = iSEList.Copy();
    cout << "Same event list was set to: " << SEList.Data() << endl;
}
void StStrangenessAna::setMEList(const TString iMEList)
{
    MEList = iMEList.Copy();
    cout << "Mixed event list was set to: " << MEList.Data() << endl;
}
void StStrangenessAna::setStopEvent_SE(const Long64_t StopEvent_SE)
{
    nStopEvent_SE = StopEvent_SE;
    cout << "nStopEvent_SE = " << nStopEvent_SE << endl;
}
void StStrangenessAna::setStartEvent_SE(const Long64_t StartEvent_SE)
{
    nStartEvent_SE = StartEvent_SE;
    cout << "nStartEvent_SE = " << nStartEvent_SE << endl;
}
void StStrangenessAna::setStopEvent_ME(const Long64_t StopEvent_ME)
{
    nStopEvent_ME = StopEvent_ME;
    cout << "nStopEvent_ME = " << nStopEvent_ME << endl;
}
void StStrangenessAna::setStartEvent_ME(const Long64_t StartEvent_ME)
{
    nStartEvent_ME = StartEvent_ME;
    cout << "nStartEvent_ME = " << nStartEvent_ME << endl;
}

void StStrangenessAna::init()
{

    cout << "Initializing parameters and input/output" << endl;
    Outputfile = new TFile(poutputfile.Data(),"RECREATE");

    refmultCorr = new StRefMultCorr();

    ALEXV0_EVENT_TREE       = (char*)V0_TREE[eAnalysis].Data();
    ALEXV0_EVENT_BRANCH     = (char*)V0_BRANCH[eAnalysis].Data();

    //----------------------------------------------------------------------------------------------------
    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( ALEXV0_EVENT_TREE, ALEXV0_EVENT_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = pinputdir+addfile;
                    input_SE ->AddFile(addfile.Data(),-1, ALEXV0_EVENT_TREE );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
            SE_input_flag = 0;
        }
    }
    // Set the input tree
    if (SE_input_flag == 1 && (eAnalysis <= 4 || eAnalysis == 8) && !input_SE->GetBranch( ALEXV0_EVENT_BRANCH ))
    {
        cerr << "ERROR: Could not find branch '"
            << ALEXV0_EVENT_BRANCH << "'in tree!" << endl;
    }

    if(eAnalysis == 0 || eAnalysis == 8)
    {
        alexPhiMeson_event = new StAlexPhiMesonEvent();
    }
    if(eAnalysis >= 1 && eAnalysis <= 4)
    {
        alexV0_event       = new StAlexV0Event();
    }
    if(SE_input_flag == 1)
    {
        {
            if(eAnalysis == 0 || eAnalysis == 8)
            {
                input_SE  ->SetBranchAddress( ALEXV0_EVENT_BRANCH, &alexPhiMeson_event );
            }
            if(eAnalysis >= 1 && eAnalysis <= 4)
            {
                input_SE  ->SetBranchAddress( ALEXV0_EVENT_BRANCH, &alexV0_event );
            }

            num_events_SE = input_SE->GetEntriesFast();
            cout << "Number of events in file(s) = " << num_events_SE << endl;
            if(nStartEvent_SE > num_events_SE) nStartEvent_SE = num_events_SE;
            if(nStopEvent_SE  > num_events_SE) nStopEvent_SE  = num_events_SE;

            cout << "New nStartEvent_SE = " << nStartEvent_SE << ", new nStopEvent_SE = " << nStopEvent_SE << endl;
        }
    }

    // Mixed event input
    if (!MEList.IsNull())   // if input file is ok
    {
        cout << "Open mixed event file list " << MEList << endl;
        ifstream in(MEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_ME  = new TChain( ALEXV0_EVENT_TREE, ALEXV0_EVENT_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = pinputdir+addfile;
                    input_ME ->AddFile(addfile.Data(),-1, ALEXV0_EVENT_TREE );
                    Long64_t file_entries = input_ME->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
        }
        else
        {
            cout << "WARNING: ME file input is problemtic" << endl;
            ME_input_flag = 0;
        }
    }
    // Set the input tree
    cout << "Read ME branch" << endl;
    if (ME_input_flag == 1 && (eAnalysis <= 4 || eAnalysis == 8)  && !input_ME->GetBranch( ALEXV0_EVENT_BRANCH ))
    {
        cerr << "ERROR: Could not find branch '"
            << ALEXV0_EVENT_BRANCH << "'in tree!" << endl;
    }

    if(ME_input_flag == 1)
    {
        num_events_ME = input_ME->GetEntriesFast();
        cout << "Number of mixed events in file(s) = " << num_events_ME << endl;
        if(nStartEvent_ME > num_events_ME) nStartEvent_ME = num_events_ME;
        if(nStopEvent_ME  > num_events_ME) nStopEvent_ME  = num_events_ME;
    }
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    cout << "Initialize re-centering correction" << endl;
    re_centering = TFile::Open(recentfile.Data());  // open the file
    h_runId_index_rc = (TH1F*)re_centering->FindObjectAny("h_runId_index");
    h_runId_index_rc->SetName("h_runId_index_rc");

    Double_t h_min_rc   = h_runId_index_rc->GetMinimum();
    Double_t h_max_rc   = h_runId_index_rc->GetMaximum();
    const Int_t  h_Nbins_rc = h_runId_index_rc->GetNbinsX();
    h_Nbins_rc_B = h_Nbins_rc+1;
    Int_t    h_Nbins_rc_inv = (Int_t)(h_max_rc-h_min_rc)+1;
    cout << "Min file id = "  << h_min_rc << ", max file id = " << h_max_rc << ", number of indices = " << h_Nbins_rc << endl;
    h_runId_index_rc_inv = new TH1F("h_runId_index_rc_inv","h_runId_index_rc_inv",h_Nbins_rc_inv,h_min_rc,h_max_rc+1.0);

    for(Int_t i_rc = 1; i_rc < h_Nbins_rc+1; i_rc++)
    {
        h_runId_index_rc_inv->SetBinContent(h_runId_index_rc_inv->FindBin(h_runId_index_rc->GetBinContent(i_rc)),(Float_t)i_rc);
    }


    // Open the mean Qx, Qy histograms for re-centering correction
    // [Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full][z-vertex bin][polynomial fit parameters] vs. index
    // h_rc_QxQy_etapm_z_vs_index_[6][n_z_vertex_bins][n_params]; // stores the fit parameters as a function of the index
    cout << "Open re-centering correction histograms" << endl;
    for(Int_t n_qxy = 0; n_qxy < 6; n_qxy++)
    {
        for(Int_t n_z = 0; n_z < n_z_vertex_bins; n_z++)
        {
            for(Int_t n_par = 0; n_par < n_params; n_par++)
            {
                HistName = "h_rc_QxQy_etapm_z_vs_index_";
                HistName += n_qxy;
                HistName += "_";
                HistName += n_z;
                HistName += "_";
                HistName += n_par;
                h_rc_QxQy_etapm_z_vs_index_[n_qxy][n_z][n_par] = (TH1F*)re_centering->FindObjectAny(HistName.Data());
            }
        }
    }

    cout << "Initialize shift correction" << endl;
    shift_method = TFile::Open(shiftfile.Data());  // open the file
    cout << "Read shift parameter histograms" << endl;
    for(Int_t jPsi_EP = 0; jPsi_EP < nPsi_EP; jPsi_EP++)
    {
        for(Int_t n_spar = 0; n_spar < n_shift_par; n_spar++)
        {
            HistName = "hEP_shift_params_";
            HistName += jPsi_EP;
            HistName += "_";
            HistName += n_spar;
            hEP_shift_params[jPsi_EP][n_spar] = (TH1F*)shift_method->FindObjectAny(HistName.Data());
        }
    }
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    cout << "Initialize histograms" << endl;
    Int_t   nAna_InvMass_bins[N_Analysis]       = {200,1250,1250,900,900,200,200,200,200};
    Float_t nAna_InvMass_start[N_Analysis]      = {0.98,1.25,1.25,1.6,1.6,1.07,1.07,0.4,0.0};
    Float_t nAna_InvMass_stop[N_Analysis]       = {1.16,2.5,2.5,2.5,2.5,1.19,1.19,0.6,4.0};
    Int_t   nAna_InvMass_miss_bins[N_Analysis]  = {200,900,900,1250,1250,200,200,200,200};
    Float_t nAna_InvMass_miss_start[N_Analysis] = {0.98,1.6,1.6,1.25,1.25,0.4,0.4,1.07,0.0};
    Float_t nAna_InvMass_miss_stop[N_Analysis]  = {1.16,2.5,2.5,2.5,2.5,0.6,0.6,1.19,1.0};

    for(Int_t SE_ME_loop = 0; SE_ME_loop < 2; SE_ME_loop++) // 0 = same event, 1 = mixed event
    {
        for(Int_t pt_bin = 0; pt_bin < npt_bins; pt_bin++)
        {
            HistName = "hInvMass_pt_SE_ME_";
            HistName += SE_ME_loop;
            HistName += "_pt_";
            HistName += pt_bin;
            hInvMass_pt[SE_ME_loop][pt_bin] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_bins[eAnalysis],nAna_InvMass_start[eAnalysis],nAna_InvMass_stop[eAnalysis]);
            hInvMass_pt[SE_ME_loop][pt_bin] ->Sumw2();

            HistName = "hInvMass_pt_miss_SE_ME_";
            HistName += SE_ME_loop;
            HistName += "_pt_";
            HistName += pt_bin;
            hInvMass_pt_miss[SE_ME_loop][pt_bin] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_miss_bins[eAnalysis],nAna_InvMass_miss_start[eAnalysis],nAna_InvMass_miss_stop[eAnalysis]);
            hInvMass_pt_miss[SE_ME_loop][pt_bin] ->Sumw2();
            for(Int_t cut = 0; cut < n_cuts; cut++)
            {
                // For systematic error calculation
                HistName = "hInvMass_pt_syst_SE_ME_";
                HistName += SE_ME_loop;
                HistName += "_pt_";
                HistName += pt_bin;
                HistName += "_cut_";
                HistName += cut;
                hInvMass_pt_syst[SE_ME_loop*n_cuts + cut][pt_bin] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_bins[eAnalysis],nAna_InvMass_start[eAnalysis],nAna_InvMass_stop[eAnalysis]);
                hInvMass_pt_syst[SE_ME_loop*n_cuts + cut][pt_bin] ->Sumw2();

                HistName = "hInvMass_pt_miss_syst_SE_ME_";
                HistName += SE_ME_loop;
                HistName += "_pt_";
                HistName += pt_bin;
                HistName += "_cut_";
                HistName += cut;
                hInvMass_pt_miss_syst[SE_ME_loop*n_cuts + cut][pt_bin] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_miss_bins[eAnalysis],nAna_InvMass_miss_start[eAnalysis],nAna_InvMass_miss_stop[eAnalysis]);
                hInvMass_pt_miss_syst[SE_ME_loop*n_cuts + cut][pt_bin] ->Sumw2();
            }
        }
    }

    for(Int_t mult_bin = 0; mult_bin < nmult_bins; mult_bin++)
    {
        for(Int_t pt_bin = 0; pt_bin < npt_bins; pt_bin++)
        {
            for(Int_t SE_ME_loop = 0; SE_ME_loop < 2; SE_ME_loop++) // 0 = same event, 1 = mixed event
            {
                HistName = "hInvMass_mult_pt_seme_";
                HistName += mult_bin;
                HistName += "_";
                HistName += pt_bin;
                HistName += "_";
                HistName += SE_ME_loop;
                hInvMass_mult_pt_seme[mult_bin][pt_bin][SE_ME_loop] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_bins[eAnalysis],nAna_InvMass_start[eAnalysis],nAna_InvMass_stop[eAnalysis]);
                hInvMass_mult_pt_seme[mult_bin][pt_bin][SE_ME_loop] ->Sumw2();

                // For systematic error calculation
                for(Int_t cut = 0; cut < n_cuts; cut++)
                {
                    HistName = "hInvMass_mult_pt_seme_syst_";
                    HistName += mult_bin;
                    HistName += "_";
                    HistName += pt_bin;
                    HistName += "_";
                    HistName += SE_ME_loop;
                    HistName += "_cut_";
                    HistName += cut;
                    hInvMass_mult_pt_seme_syst[mult_bin][pt_bin][SE_ME_loop*n_cuts + cut] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_bins[eAnalysis],nAna_InvMass_start[eAnalysis],nAna_InvMass_stop[eAnalysis]);
                    hInvMass_mult_pt_seme_syst[mult_bin][pt_bin][SE_ME_loop*n_cuts + cut] ->Sumw2();
                }

                //cout << "mult_bin = " << mult_bin << ", pt_bin = " << pt_bin << ", SE_ME_loop = " << SE_ME_loop << endl;
                if(mult_bin < nmult_bins4)
                {
                    for(Int_t phi_bin = 0; phi_bin < nphi_bins; phi_bin++)
                    {
                        for(Int_t m_ep = 0; m_ep < nep_methods; m_ep++)
                        {

                            HistName = "hInvMass_mult_pt_phi_ep_seme_";
                            HistName += mult_bin;
                            HistName += "_";
                            HistName += pt_bin;
                            HistName += "_";
                            HistName += phi_bin;
                            HistName += "_";
                            HistName += m_ep;
                            HistName += "_";
                            HistName += SE_ME_loop;
                            hInvMass_mult_pt_phi_ep_seme[mult_bin][pt_bin][phi_bin][m_ep][SE_ME_loop] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_bins[eAnalysis],nAna_InvMass_start[eAnalysis],nAna_InvMass_stop[eAnalysis]);
                            hInvMass_mult_pt_phi_ep_seme[mult_bin][pt_bin][phi_bin][m_ep][SE_ME_loop] ->Sumw2();

                            // For systematic error calculation
                            for(Int_t cut = 0; cut < n_cuts; cut++)
                            {
                                HistName = "hInvMass_mult_pt_phi_ep_seme_syst_";
                                HistName += mult_bin;
                                HistName += "_";
                                HistName += pt_bin;
                                HistName += "_";
                                HistName += phi_bin;
                                HistName += "_";
                                HistName += m_ep;
                                HistName += "_";
                                HistName += SE_ME_loop;
                                HistName += "_cut_";
                                HistName += cut;
                                //cout << "phi_bin = " << phi_bin << ", m_ep =" << m_ep << ", cut = " << cut << endl;
                                hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][pt_bin][phi_bin][m_ep][SE_ME_loop*n_cuts + cut] = new TH1F(HistName.Data(),HistName.Data(),nAna_InvMass_bins[eAnalysis],nAna_InvMass_start[eAnalysis],nAna_InvMass_stop[eAnalysis]);
                                hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][pt_bin][phi_bin][m_ep][SE_ME_loop*n_cuts + cut] ->Sumw2();
                            }

                        }
                    }
                }
            }
        }
    }

    if(eAnalysis == 8)
    {
        NTuple_Four_leptons      = new TNtuple("NTuple_Four_leptons","NTuple_Four_leptons Ntuple","InvABCD:pt:rap:mult:z:m2A:m2B:nsA:nsB:InvAB:qpA:qpB:m2C:m2D:nsC:nsD:InvCD:qpC:qpD");
        NTuple_Four_leptons      ->SetAutoSave( 5000000 );
    }
 



    //----------------------------------------------------------------------------------------------------



    //*********************************************************************************************************
    // Event plane resolutions
    cout << "Initialize event plane resolutions" << endl;
    if(eBeamTimeNum == 0) // 7.7 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_7GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_7GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_7GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_7GeV[e]*0.01;
        }
    }
    if(eBeamTimeNum == 1) // 11.5 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_11GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_11GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_11GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_11GeV[e]*0.01;
        }
    }
    if(eBeamTimeNum == 2) // 39 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_39GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_39GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_39GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_39GeV[e]*0.01;
        }
    }
    if(eBeamTimeNum == 3) // 62.4 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_62GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_62GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_62GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_62GeV[e]*0.01;
        }
    }
    if(eBeamTimeNum == 4) // 19.6 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_19GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_19GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_19GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_19GeV[e]*0.01;
        }
    }
    if(eBeamTimeNum == 5) // 27 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_27GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_27GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_27GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_27GeV[e]*0.01;
        }
    }
    if(eBeamTimeNum == 6) // 200 GeV
    {
        for(Int_t e = 0; e < 9; e++)
        {
            Event_Plane_resolution[e][0] = Resolution_EP0_7GeV[e]*0.01;
            Event_Plane_resolution[e][1] = Resolution_EP1_7GeV[e]*0.01;
            Event_Plane_resolution[e][2] = Resolution_EP2_7GeV[e]*0.01;
            Event_Plane_resolution[e][3] = Resolution_EP3_7GeV[e]*0.01;
        }
    }
    //*********************************************************************************************************



}

void StStrangenessAna::loop()
{

    cout << "Define bad run number array" << endl;
    Int_t bad_run_numbers[n_bad_run_numbers[eBeamTimeNum]];
    for(Int_t i = 0; i < n_bad_run_numbers[eBeamTimeNum]; i++)
    {
        if(eBeamTimeNum == 0) bad_run_numbers[i] = bad_run_list_7GeV[i];
        if(eBeamTimeNum == 1) bad_run_numbers[i] = bad_run_list_11GeV[i];
        if(eBeamTimeNum == 2) bad_run_numbers[i] = bad_run_list_39GeV[i];
        if(eBeamTimeNum == 3) bad_run_numbers[i] = bad_run_list_62GeV[i];
        if(eBeamTimeNum == 4) bad_run_numbers[i] = bad_run_list_19GeV[i];
        if(eBeamTimeNum == 5) bad_run_numbers[i] = bad_run_list_27GeV[i];
        if(eBeamTimeNum == 6) bad_run_numbers[i] = bad_run_list_200GeV[i];
    }

    Long64_t start_event_use = 0;
    Long64_t stop_event_use  = 0;

    Long64_t accepted_SE_ME_counter[2] = {0,0};

    //----------------------------------------------------------------------------------------------------
    // Determine nSigma scaling factor
    Double_t scale_nSigma_fac = 1.0;
    if(eBeamTimeNum == 5) scale_nSigma_fac = 1.9; // 27 GeV has a wrong nSigma calibration
    //----------------------------------------------------------------------------------------------------



    //******************** Cuts *****************************************
    Double_t dcaBC_cut           = 1.0; // 1.0
    Double_t dcaAB_cut           = 1.0; // 1.0
    Double_t lambda_low_cut      = 1.1157-3.0*0.00129;
    Double_t lambda_high_cut     = 1.1157+3.0*0.00129;
    Double_t pion_low_cut        = -0.02; //
    Double_t pion_high_cut       = 0.05;   //
    Double_t proton_low_cut      = 0.6;   //
    Double_t proton_high_cut     = 1.4;   //
    Double_t kaon_low_cut        = 0.125;   //
    Double_t kaon_high_cut       = 0.36;   //
    //*******************************************************************



    //----------------------------------------------------------------------------------------------------

    //--------------------------------
    // Event based variables
    Float_t   EventVertexX        = 0;
    Float_t   EventVertexY        = 0;
    Float_t   EventVertexZ        = 0;
    Int_t     RunId               = 0;
    Float_t   refMult             = 0;
    Int_t     n_prim              = 0;
    Int_t     n_non_prim          = 0;
    Int_t     n_tof_prim          = 0;
    Float_t   EP_Qx_eta_pos_ptw   = 0;
    Float_t   EP_Qy_eta_pos_ptw   = 0;
    Float_t   EP_Qx_eta_neg_ptw   = 0;
    Float_t   EP_Qy_eta_neg_ptw   = 0;
    Float_t   EP_Qx_ptw           = 0;
    Float_t   EP_Qy_ptw           = 0;
    Int_t     Qtracks_eta_pos     = 0;
    Int_t     Qtracks_eta_neg     = 0;
    Int_t     Qtracks_full        = 0;
    Float_t   ZDCx                = 0;  // ZDC coincidence rate
    Float_t   BBCx                = 0;
    Float_t   vzVpd               = 0;
    Int_t   erun_nevents          = 0;
    //--------------------------------


    //--------------------------------
    // Track based variables
    Float_t m2A, m2B,nsA, nsB, dcaA, dcaB, iQxA, iQyA, iQxB, iQyB, etaA, etaB, InvAB, p_t, rap, phi, theta, phi_event_plane;
    Float_t phi_event_plane_eta_gap, Psi_diff_ME, Delta_dip_angle, qpA, qpB, qpC;
    Float_t m2C, nsC, dcaC, iQxC, iQyC, etaC, InvABC, InvAB_miss, InvABC_miss, dcaAB, dcaBC, VerdistX, VerdistY;
    Float_t VerdistX2, VerdistY2, scal_prod, scal_prod2, InvMass_Fill, InvMass_Fill_miss;
    //--------------------------------

    // Variables for Lambdas and K0S Ntuples
    //Float_t pA,pB;
    Float_t RunId_f,n_prim_f,n_non_prim_f,n_tof_prim_f,Qtracks_eta_pos_f,Qtracks_eta_neg_f,Qtracks_full_f;
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    // Event loop
    for(Int_t SE_ME_loop = 0; SE_ME_loop < 2; SE_ME_loop++) // 0 = same event, 1 = mixed event
    {
        start_event_use = 0;
        stop_event_use  = 0;
        if(SE_ME_loop == 0 && SE_input_flag == 1)
        {
            if(eAnalysis == 0 || eAnalysis == 8) // phi meson analysis
            {
                input_SE  ->SetBranchAddress( ALEXV0_EVENT_BRANCH, &alexPhiMeson_event );
            }
            if(eAnalysis >= 1 && eAnalysis <= 4) // Xi-/+, Omega-/+ analysis
            {
                input_SE  ->SetBranchAddress( ALEXV0_EVENT_BRANCH, &alexV0_event );
            }
            if(eAnalysis >= 5 && eAnalysis < 8) // Lambda, anti-Lambda, K0S analysis
            {
                input_SE->SetBranchAddress("InvM",&InvAB);
                input_SE->SetBranchAddress("InvM_K0S",&InvAB_miss);
                input_SE->SetBranchAddress("m2A",&m2A);
                input_SE->SetBranchAddress("m2B",&m2B);
                input_SE->SetBranchAddress("nSP",&nsA); // proton
                input_SE->SetBranchAddress("nSPi",&nsB); // pion
                input_SE->SetBranchAddress("pA",&qpA);
                input_SE->SetBranchAddress("pB",&qpB);
                input_SE->SetBranchAddress("dcaA",&dcaA);
                input_SE->SetBranchAddress("dcaB",&dcaB);
                input_SE->SetBranchAddress("iQxA",&iQxA);
                input_SE->SetBranchAddress("iQyA",&iQyA);
                input_SE->SetBranchAddress("iQxB",&iQxB);
                input_SE->SetBranchAddress("iQyB",&iQyB);
                input_SE->SetBranchAddress("etaA",&etaA);
                input_SE->SetBranchAddress("etaB",&etaB);
                input_SE->SetBranchAddress("mult",&refMult);
                input_SE->SetBranchAddress("dcaAB",&dcaAB);
                input_SE->SetBranchAddress("VerdistX",&VerdistX);
                input_SE->SetBranchAddress("VerdistY",&VerdistY);
                input_SE->SetBranchAddress("pt",&p_t);
                input_SE->SetBranchAddress("rap",&rap);
                input_SE->SetBranchAddress("phiAB",&phi);
                input_SE->SetBranchAddress("thetaAB",&theta);
                input_SE->SetBranchAddress("phi_EP",&phi_event_plane);
                input_SE->SetBranchAddress("phi_EP_eta_gap",&phi_event_plane_eta_gap);
                input_SE->SetBranchAddress("delta_phi_ME_AB",&Psi_diff_ME);
                input_SE->SetBranchAddress("RunId",&RunId_f);
                input_SE->SetBranchAddress("n_prim",&n_prim_f);
                input_SE->SetBranchAddress("n_non_prim",&n_non_prim_f);
                input_SE->SetBranchAddress("n_tof_prim",&n_tof_prim_f);
                input_SE->SetBranchAddress("scalarProduct",&scal_prod);
                input_SE->SetBranchAddress("x",&EventVertexX);
                input_SE->SetBranchAddress("y",&EventVertexY);
                input_SE->SetBranchAddress("z",&EventVertexZ);
                input_SE->SetBranchAddress("EP_Qx_eta_pos_ptw",&EP_Qx_eta_pos_ptw);
                input_SE->SetBranchAddress("EP_Qy_eta_pos_ptw",&EP_Qy_eta_pos_ptw);
                input_SE->SetBranchAddress("EP_Qx_eta_neg_ptw",&EP_Qx_eta_neg_ptw);
                input_SE->SetBranchAddress("EP_Qy_eta_neg_ptw",&EP_Qy_eta_neg_ptw);
                input_SE->SetBranchAddress("EP_Qx_ptw",&EP_Qx_ptw);
                input_SE->SetBranchAddress("EP_Qy_ptw",&EP_Qy_ptw);
                input_SE->SetBranchAddress("Qtracks_eta_pos",&Qtracks_eta_pos_f);
                input_SE->SetBranchAddress("Qtracks_eta_neg",&Qtracks_eta_neg_f);
                input_SE->SetBranchAddress("Qtracks_all",&Qtracks_full_f);
                input_SE->SetBranchAddress("ZDCx",&ZDCx);
                input_SE->SetBranchAddress("vzVpd",&vzVpd);
            }
            cout << "" << endl;
            cout << "------------------- Start looping: Same event -------------------" << endl;
            start_event_use = nStartEvent_SE;
            stop_event_use  = nStopEvent_SE;
            cout << "start_event_use = " << start_event_use << ", stop_event_use = " << stop_event_use << endl;
            input_SE->GetEntry( 0 ); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

            RunId           = (Int_t)RunId_f;
            n_prim          = (Int_t)n_prim_f;
            n_non_prim      = (Int_t)n_non_prim_f;
            n_tof_prim      = (Int_t)n_tof_prim_f;
            Qtracks_eta_pos = (Int_t)Qtracks_eta_pos_f;
            Qtracks_eta_neg = (Int_t)Qtracks_eta_neg_f;
            Qtracks_full    = (Int_t)Qtracks_full_f;
        }
        if(SE_ME_loop == 1 && ME_input_flag == 1)
        {
            if(eAnalysis == 0 || eAnalysis == 8) // phi meson analysis
            {
                input_ME  ->SetBranchAddress( ALEXV0_EVENT_BRANCH, &alexPhiMeson_event );
            }
            if(eAnalysis >= 1 && eAnalysis <= 4) // Xi-/+, Omega-/+ analysis
            {
                input_ME  ->SetBranchAddress( ALEXV0_EVENT_BRANCH, &alexV0_event );
            }
            if(eAnalysis >= 5 && eAnalysis < 8) // Lambda, anti-Lambda, K0S analysis
            {
                input_ME->SetBranchAddress("InvM",&InvAB);
                input_ME->SetBranchAddress("InvM_K0S",&InvAB_miss);
                input_ME->SetBranchAddress("m2A",&m2A);
                input_ME->SetBranchAddress("m2B",&m2B);
                input_ME->SetBranchAddress("nSP",&nsA); // proton
                input_ME->SetBranchAddress("nSPi",&nsB); // pion
                input_ME->SetBranchAddress("pA",&qpA);
                input_ME->SetBranchAddress("pB",&qpB);
                input_ME->SetBranchAddress("dcaA",&dcaA);
                input_ME->SetBranchAddress("dcaB",&dcaB);
                input_ME->SetBranchAddress("iQxA",&iQxA);
                input_ME->SetBranchAddress("iQyA",&iQyA);
                input_ME->SetBranchAddress("iQxB",&iQxB);
                input_ME->SetBranchAddress("iQyB",&iQyB);
                input_ME->SetBranchAddress("etaA",&etaA);
                input_ME->SetBranchAddress("etaB",&etaB);
                input_ME->SetBranchAddress("mult",&refMult);
                input_ME->SetBranchAddress("dcaAB",&dcaAB);
                input_ME->SetBranchAddress("VerdistX",&VerdistX);
                input_ME->SetBranchAddress("VerdistY",&VerdistY);
                input_ME->SetBranchAddress("pt",&p_t);
                input_ME->SetBranchAddress("rap",&rap);
                input_ME->SetBranchAddress("phiAB",&phi);
                input_ME->SetBranchAddress("thetaAB",&theta);
                input_ME->SetBranchAddress("phi_EP",&phi_event_plane);
                input_ME->SetBranchAddress("phi_EP_eta_gap",&phi_event_plane_eta_gap);
                input_ME->SetBranchAddress("delta_phi_ME_AB",&Psi_diff_ME);
                input_ME->SetBranchAddress("RunId",&RunId_f);
                input_ME->SetBranchAddress("n_prim",&n_prim_f);
                input_ME->SetBranchAddress("n_non_prim",&n_non_prim_f);
                input_ME->SetBranchAddress("n_tof_prim",&n_tof_prim_f);
                input_ME->SetBranchAddress("scalarProduct",&scal_prod);
                input_ME->SetBranchAddress("x",&EventVertexX);
                input_ME->SetBranchAddress("y",&EventVertexY);
                input_ME->SetBranchAddress("z",&EventVertexZ);
                input_ME->SetBranchAddress("EP_Qx_eta_pos_ptw",&EP_Qx_eta_pos_ptw);
                input_ME->SetBranchAddress("EP_Qy_eta_pos_ptw",&EP_Qy_eta_pos_ptw);
                input_ME->SetBranchAddress("EP_Qx_eta_neg_ptw",&EP_Qx_eta_neg_ptw);
                input_ME->SetBranchAddress("EP_Qy_eta_neg_ptw",&EP_Qy_eta_neg_ptw);
                input_ME->SetBranchAddress("EP_Qx_ptw",&EP_Qx_ptw);
                input_ME->SetBranchAddress("EP_Qy_ptw",&EP_Qy_ptw);
                input_ME->SetBranchAddress("Qtracks_eta_pos",&Qtracks_eta_pos_f);
                input_ME->SetBranchAddress("Qtracks_eta_neg",&Qtracks_eta_neg_f);
                input_ME->SetBranchAddress("Qtracks_all",&Qtracks_full_f);
                input_ME->SetBranchAddress("ZDCx",&ZDCx);
                input_ME->SetBranchAddress("vzVpd",&vzVpd);   
            }

            cout << "" << endl;
            cout << "------------------- Start looping: Mixed event -------------------" << endl;
            start_event_use = nStartEvent_ME;
            stop_event_use  = nStopEvent_ME;
            input_ME->GetEntry( 0 ); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

            RunId           = (Int_t)RunId_f;
            n_prim          = (Int_t)n_prim_f;
            n_non_prim      = (Int_t)n_non_prim_f;
            n_tof_prim      = (Int_t)n_tof_prim_f;
            Qtracks_eta_pos = (Int_t)Qtracks_eta_pos_f;
            Qtracks_eta_neg = (Int_t)Qtracks_eta_neg_f;
            Qtracks_full    = (Int_t)Qtracks_full_f;
        }

        cout << "SE_ME_loop = " << SE_ME_loop << ", start_event_use = " << start_event_use << ", stop_event_use = " << stop_event_use << ", SE_input_flag = " << SE_input_flag << endl;

        for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
        {
            if (counter != 0  &&  counter % 1000 == 0)
                cout << "." << flush;
            if (counter != 0  &&  counter % 10000 == 0)
            {
                if((stop_event_use-start_event_use) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                    cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_v2) " << flush;
                }
            }

            if(SE_ME_loop == 0)
            {
                if (!input_SE->GetEntry( counter )) // take the event -> information is stored in event
                    break;  // end of data chunk
            }
            if(SE_ME_loop == 1)
            {
                if (!input_ME->GetEntry( counter )) // take the event -> information is stored in event
                    break;  // end of data chunk
            }

            RunId           = (Int_t)RunId_f;
            n_prim          = (Int_t)n_prim_f;
            n_non_prim      = (Int_t)n_non_prim_f;
            n_tof_prim      = (Int_t)n_tof_prim_f;
            Qtracks_eta_pos = (Int_t)Qtracks_eta_pos_f;
            Qtracks_eta_neg = (Int_t)Qtracks_eta_neg_f;
            Qtracks_full    = (Int_t)Qtracks_full_f;

            //---------------------------------------------------------------------------
            if(eAnalysis == 0 || eAnalysis == 8) // phi meson analysis
            {
                erun_nevents         = alexPhiMeson_event->getNumTracks(); // number of tracks in this event
                EventVertexX         = alexPhiMeson_event->getx();
                EventVertexY         = alexPhiMeson_event->gety();
                EventVertexZ         = alexPhiMeson_event->getz();
                RunId                = alexPhiMeson_event->getid();
                refMult              = alexPhiMeson_event->getmult();
                n_prim               = alexPhiMeson_event->getn_prim();
                n_non_prim           = alexPhiMeson_event->getn_non_prim();
                n_tof_prim           = alexPhiMeson_event->getn_tof_prim();
                EP_Qx_eta_pos_ptw    = alexPhiMeson_event->getEP_Qx_eta_pos_ptw();
                EP_Qy_eta_pos_ptw    = alexPhiMeson_event->getEP_Qy_eta_pos_ptw();
                EP_Qx_eta_neg_ptw    = alexPhiMeson_event->getEP_Qx_eta_neg_ptw();
                EP_Qy_eta_neg_ptw    = alexPhiMeson_event->getEP_Qy_eta_neg_ptw();
                EP_Qx_ptw            = alexPhiMeson_event->getEP_Qx_ptw();
                EP_Qy_ptw            = alexPhiMeson_event->getEP_Qy_ptw();
                Qtracks_eta_pos      = alexPhiMeson_event->getQtracks_eta_pos();
                Qtracks_eta_neg      = alexPhiMeson_event->getQtracks_eta_neg();
                Qtracks_full         = alexPhiMeson_event->getQtracks_full();
                ZDCx                 = alexPhiMeson_event->getZDCx();
                BBCx                 = alexPhiMeson_event->getBBCx();
                vzVpd                = alexPhiMeson_event->getvzVpd();
            }

            if(eAnalysis >= 1 && eAnalysis <= 4) // Xi-/+, Omega-/+ analysis
            {
                erun_nevents         =  alexV0_event->getNumTracks(); // number of tracks in this event
                EventVertexX         =  alexV0_event->getx();
                EventVertexY         =  alexV0_event->gety();
                EventVertexZ         =  alexV0_event->getz();
                RunId                =  alexV0_event->getid();
                refMult              =  alexV0_event->getmult();
                n_prim               =  alexV0_event->getn_prim();
                n_non_prim           =  alexV0_event->getn_non_prim();
                n_tof_prim           =  alexV0_event->getn_tof_prim();
                EP_Qx_eta_pos_ptw    =  alexV0_event->getEP_Qx_eta_pos_ptw();
                EP_Qy_eta_pos_ptw    =  alexV0_event->getEP_Qy_eta_pos_ptw();
                EP_Qx_eta_neg_ptw    =  alexV0_event->getEP_Qx_eta_neg_ptw();
                EP_Qy_eta_neg_ptw    =  alexV0_event->getEP_Qy_eta_neg_ptw();
                EP_Qx_ptw            =  alexV0_event->getEP_Qx_ptw();
                EP_Qy_ptw            =  alexV0_event->getEP_Qy_ptw();
                Qtracks_eta_pos      =  alexV0_event->getQtracks_eta_pos();
                Qtracks_eta_neg      =  alexV0_event->getQtracks_eta_neg();
                Qtracks_full         =  alexV0_event->getQtracks_full();
                ZDCx                 =  alexV0_event->getZDCx();
                BBCx                 =  alexV0_event->getBBCx();
                vzVpd                =  alexV0_event->getvzVpd();
            }

            if(eAnalysis >= 5 && eAnalysis < 8) // Lambda, anti-Lambda, K0S analysis
            {
                erun_nevents      =  1;
                nsA*=scale_nSigma_fac;
                nsB*=scale_nSigma_fac;
                InvMass_Fill      = InvAB;
                InvMass_Fill_miss = InvAB_miss;
                qpC               = 0.0;
                iQxC              = 0.0;
                iQyC              = 0.0;
                etaC              = 0.0;
            }
            //---------------------------------------------------------------------------


            //---------------------------------------------------------------------------
            // Good/Bad run selection
            Int_t flag_good_run = 1;
            for(Int_t bad_run = 0; bad_run < n_bad_run_numbers[eBeamTimeNum]; bad_run++)
            {
                if(bad_run_numbers[bad_run] == (Int_t)RunId)
                {
                    flag_good_run = 0;
                    break;
                }
            }
            //---------------------------------------------------------------------------


            //---------------------------------------------------------------------------
            // Centrality determination
            refmultCorr->init((Int_t)RunId);
            refmultCorr->initEvent(refMult,EventVertexZ);
            //Double_t mult_corr   = refmultCorr->getRefMultCorr();
            Double_t reweight    = refmultCorr->getWeight();
            erefMult_bin16       = refmultCorr->getCentralityBin16();
            erefMult_bin         = refmultCorr->getCentralityBin9();
            //---------------------------------------------------------------------------


            //---------------------------------------------------------------------------
            // z-vertex bin determination
            Int_t z_bin = -1;
            Float_t start_z = -z_acceptance[eBeamTimeNum];
            Float_t delta_z = 2.0*z_acceptance[eBeamTimeNum]/((Float_t)n_z_vertex_bins);
            z_bin = (Int_t)((EventVertexZ-start_z)/delta_z); // z-vertex bin for re-centering correction
            //---------------------------------------------------------------------------


            Int_t day, file_of_day, file_id;
            Get_day_file_id((Int_t)RunId,day,file_of_day,file_id);
            Int_t file_bin_rc = (Int_t)(h_runId_index_rc_inv->GetBinContent(h_runId_index_rc_inv->FindBin(file_id))); // index bin for re-centering correction


            //----------------------------------------------------------------------------------------------------
            // Determine re-centering correction
            // h_rc_QxQy_etapm_z_vs_index_[q_id][z_bin][p_id]:  [Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full][z-bin][parameter]
            Float_t corr_params_rc[6][3]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full, [p_id]: par0, par1, par2

            for(Int_t q_id = 0; q_id < 6; q_id++)  // q_id: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
            {
                for(Int_t p_id = 0; p_id < 3; p_id++) // p_id: par0, par1, par2
                {
                    corr_params_rc[q_id][p_id] = h_rc_QxQy_etapm_z_vs_index_[q_id][z_bin][p_id]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[q_id][z_bin][p_id]->FindBin(file_bin_rc));
                }
            }


            // Correction parameters (always normalized to number of tracks/event)
            Float_t corr_Q_vectors_rc[6]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
            Float_t EP_Q_vectors[6]          = {EP_Qx_eta_pos_ptw,EP_Qy_eta_pos_ptw,EP_Qx_eta_neg_ptw,EP_Qy_eta_neg_ptw,EP_Qx_ptw,EP_Qy_ptw}; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
            Float_t N_tracks_EP_Q_vectors[6] = {(Float_t)Qtracks_eta_pos,(Float_t)Qtracks_eta_pos,(Float_t)Qtracks_eta_neg,(Float_t)Qtracks_eta_neg,(Float_t)Qtracks_full,(Float_t)Qtracks_full}; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
            // Re-centering corrected Q-vectors
            Float_t EP_Q_vectors_rc[6]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
            for(Int_t q_id = 0; q_id < 6; q_id++)  // q_id: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
            {
                corr_Q_vectors_rc[q_id] = corr_params_rc[q_id][0] + corr_params_rc[q_id][1]*refMult + corr_params_rc[q_id][2]*refMult*refMult;
                if(N_tracks_EP_Q_vectors[q_id] > 0.0)
                {
                    EP_Q_vectors_rc[q_id] = (EP_Q_vectors[q_id]/N_tracks_EP_Q_vectors[q_id]) - corr_Q_vectors_rc[q_id];
                }
                else
                {
                    EP_Q_vectors_rc[q_id] = EP_Q_vectors[q_id];
                }
            }
            //----------------------------------------------------------------------------------------------------



            //----------------------------------------------------------------------------------------------------
            if(
               flag_good_run == 1
               && EventVertexZ > -z_acceptance[eBeamTimeNum]
               && EventVertexZ < z_acceptance[eBeamTimeNum]
               && n_tof_prim >= 2
               && n_non_prim < 3000
               && erefMult_bin >= 0
               && erefMult_bin < 9
              )
            {
                if(eAnalysis == 8) // pi0 analysis
                {
                    for(UShort_t i = 0; i < erun_nevents; ++i) // loop over all tracks of the actual event
                    {
                        alexPhiMeson_track      = alexPhiMeson_event->getTrack( i ); // take the track
                        m2A                     = alexPhiMeson_track->getm2A();
                        m2B                     = alexPhiMeson_track->getm2B();
                        nsA                     = alexPhiMeson_track->getnsA()*scale_nSigma_fac;
                        nsB                     = alexPhiMeson_track->getnsB()*scale_nSigma_fac;
                        dcaA                    = alexPhiMeson_track->getdcaA();
                        dcaB                    = alexPhiMeson_track->getdcaB();
                        iQxA                    = alexPhiMeson_track->getiQxA();
                        iQyA                    = alexPhiMeson_track->getiQyA();
                        iQxB                    = alexPhiMeson_track->getiQxB();
                        iQyB                    = alexPhiMeson_track->getiQyB();
                        etaA                    = alexPhiMeson_track->getetaA();
                        etaB                    = alexPhiMeson_track->getetaB();
                        InvAB                   = alexPhiMeson_track->getInvAB();
                        p_t                     = alexPhiMeson_track->getpt();
                        rap                     = alexPhiMeson_track->getrap(); // true Rapidity, not pseudo-Rapidity!
                        phi                     = alexPhiMeson_track->getphi();
                        theta                   = alexPhiMeson_track->gettheta();
                        phi_event_plane         = alexPhiMeson_track->getPsi_ep();
                        phi_event_plane_eta_gap = alexPhiMeson_track->getPsi_ep_eta();
                        Psi_diff_ME             = alexPhiMeson_track->getPsi_diff_ME();
                        Delta_dip_angle         = alexPhiMeson_track->getDelta_dip_angle();
                        qpA                     = alexPhiMeson_track->getqpA();
                        qpB                     = alexPhiMeson_track->getqpB();
                        qpC                     = 0.0;
                        iQxC                    = 0.0;
                        iQyC                    = 0.0;
                        etaC                    = 0.0;

                        InvMass_Fill      = InvAB;
                        InvMass_Fill_miss = InvAB;
                        Double_t m2A_A    = m2A;
                        Double_t m2B_A    = m2B;
                        Double_t nsA_A    = nsA;
                        Double_t nsB_A    = nsB;
                        Double_t InvAB_A  = InvAB;
                        Double_t qpA_A    = qpA;
                        Double_t qpB_A    = qpB;

                        if(InvAB < 0.07 && InvAB > 0.005)
                        {
                            Double_t gamma_p = TMath::Abs(qpA) + TMath::Abs(qpB); // total momentum of the two electrons is the gamma momentum
                            TLorentzVector trackAB;
                            TVector3 trackAB_vec3;
                            trackAB_vec3.SetMagThetaPhi(gamma_p,theta,phi);
                            Double_t E_tot = TMath::Sqrt((2.0*Mass_electron)*(2.0*Mass_electron) + gamma_p*gamma_p);
                            trackAB.SetPxPyPzE(trackAB_vec3.Px(),trackAB_vec3.Py(),trackAB_vec3.Pz(),E_tot);
                            //cout << "i = " << i << ", px = " << trackAB.E() << endl;

                            for(UShort_t j = i+1; j < erun_nevents; ++j) // loop over all tracks of the actual event
                            {
                                alexPhiMeson_track      = alexPhiMeson_event->getTrack( j ); // take the track
                                m2A                     = alexPhiMeson_track->getm2A();
                                m2B                     = alexPhiMeson_track->getm2B();
                                nsA                     = alexPhiMeson_track->getnsA()*scale_nSigma_fac;
                                nsB                     = alexPhiMeson_track->getnsB()*scale_nSigma_fac;
                                dcaA                    = alexPhiMeson_track->getdcaA();
                                dcaB                    = alexPhiMeson_track->getdcaB();
                                iQxA                    = alexPhiMeson_track->getiQxA();
                                iQyA                    = alexPhiMeson_track->getiQyA();
                                iQxB                    = alexPhiMeson_track->getiQxB();
                                iQyB                    = alexPhiMeson_track->getiQyB();
                                etaA                    = alexPhiMeson_track->getetaA();
                                etaB                    = alexPhiMeson_track->getetaB();
                                InvAB                   = alexPhiMeson_track->getInvAB();
                                p_t                     = alexPhiMeson_track->getpt();
                                rap                     = alexPhiMeson_track->getrap(); // true Rapidity, not pseudo-Rapidity!
                                phi                     = alexPhiMeson_track->getphi();
                                theta                   = alexPhiMeson_track->gettheta();
                                phi_event_plane         = alexPhiMeson_track->getPsi_ep();
                                phi_event_plane_eta_gap = alexPhiMeson_track->getPsi_ep_eta();
                                Psi_diff_ME             = alexPhiMeson_track->getPsi_diff_ME();
                                Delta_dip_angle         = alexPhiMeson_track->getDelta_dip_angle();
                                qpA                     = alexPhiMeson_track->getqpA();
                                qpB                     = alexPhiMeson_track->getqpB();
                                qpC                     = 0.0;
                                iQxC                    = 0.0;
                                iQyC                    = 0.0;
                                etaC                    = 0.0;

                                if(InvAB < 0.07 && InvAB > 0.005)
                                {
                                    Double_t gamma_p = TMath::Abs(qpA) + TMath::Abs(qpB); // total momentum of the two electrons is the gamma momentum
                                    TLorentzVector trackCD;
                                    TVector3 trackCD_vec3;
                                    trackCD_vec3.SetMagThetaPhi(gamma_p,theta,phi);
                                    Double_t E_tot = TMath::Sqrt((2.0*Mass_electron)*(2.0*Mass_electron) + gamma_p*gamma_p);
                                    trackCD.SetPxPyPzE(trackCD_vec3.Px(),trackCD_vec3.Py(),trackCD_vec3.Pz(),E_tot);

                                    TLorentzVector trackABCD = trackAB + trackCD;

                                    //cout << "p1 = {" << trackAB.Px() << ", " << trackAB.Py() << ", " << trackAB.Pz() << "}, E1 = " << trackAB.E() << endl;
                                    //cout << "p2 = {" << trackCD.Px() << ", " << trackCD.Py() << ", " << trackCD.Pz() << "}, E2 = " << trackCD.E() << endl;
                                    //cout << "InvAB = " << InvAB << ", tracCD.M = " << trackCD.M() << ", trackABCD.M = " << trackABCD.M() << ", gamma_p = " << gamma_p << endl;

                                    //cout << "i = " << i << ", j = " << j << ", InvABCD = " << trackABCD.M() << ", theta = " << theta << ", phi = " << phi
                                    //    << ", pA = " << trackAB.P() << ", pB = " << trackCD.P() << endl;

                                    // Filling Ntuple
                                    Four_leptons_NTDataArray[0]      =(Float_t)trackABCD.M();
                                    Four_leptons_NTDataArray[1]      =(Float_t)trackABCD.Pt();
                                    Four_leptons_NTDataArray[2]      =(Float_t)trackABCD.Rapidity();
                                    Four_leptons_NTDataArray[3]      =(Float_t)refMult;
                                    Four_leptons_NTDataArray[4]      =(Float_t)EventVertexZ;
                                    Four_leptons_NTDataArray[5]      =(Float_t)m2A_A;
                                    Four_leptons_NTDataArray[6]      =(Float_t)m2B_A;
                                    Four_leptons_NTDataArray[7]      =(Float_t)nsA_A;
                                    Four_leptons_NTDataArray[8]      =(Float_t)nsB_A;
                                    Four_leptons_NTDataArray[9]      =(Float_t)InvAB_A;
                                    Four_leptons_NTDataArray[10]     =(Float_t)qpA_A;
                                    Four_leptons_NTDataArray[11]     =(Float_t)qpB_A;
                                    Four_leptons_NTDataArray[12]     =(Float_t)m2A;
                                    Four_leptons_NTDataArray[13]     =(Float_t)m2B;
                                    Four_leptons_NTDataArray[14]     =(Float_t)nsA;
                                    Four_leptons_NTDataArray[15]     =(Float_t)nsB;
                                    Four_leptons_NTDataArray[16]     =(Float_t)InvAB;
                                    Four_leptons_NTDataArray[17]     =(Float_t)qpA;
                                    Four_leptons_NTDataArray[18]     =(Float_t)qpB;

                                    NTuple_Four_leptons->Fill(Four_leptons_NTDataArray);
                                }
                            }
                        }

                    }
                }
                else
                {
                    for(UShort_t i = 0; i < erun_nevents; ++i) // loop over all tracks of the actual event
                    {
                        if(eAnalysis == 0) // phi meson analysis
                        {
                            alexPhiMeson_track      = alexPhiMeson_event->getTrack( i ); // take the track
                            m2A                     = alexPhiMeson_track->getm2A();
                            m2B                     = alexPhiMeson_track->getm2B();
                            nsA                     = alexPhiMeson_track->getnsA()*scale_nSigma_fac;
                            nsB                     = alexPhiMeson_track->getnsB()*scale_nSigma_fac;
                            dcaA                    = alexPhiMeson_track->getdcaA();
                            dcaB                    = alexPhiMeson_track->getdcaB();
                            iQxA                    = alexPhiMeson_track->getiQxA();
                            iQyA                    = alexPhiMeson_track->getiQyA();
                            iQxB                    = alexPhiMeson_track->getiQxB();
                            iQyB                    = alexPhiMeson_track->getiQyB();
                            etaA                    = alexPhiMeson_track->getetaA();
                            etaB                    = alexPhiMeson_track->getetaB();
                            InvAB                   = alexPhiMeson_track->getInvAB();
                            p_t                     = alexPhiMeson_track->getpt();
                            rap                     = alexPhiMeson_track->getrap(); // true Rapidity, not pseudo-Rapidity!
                            phi                     = alexPhiMeson_track->getphi();
                            theta                   = alexPhiMeson_track->gettheta();
                            phi_event_plane         = alexPhiMeson_track->getPsi_ep();
                            phi_event_plane_eta_gap = alexPhiMeson_track->getPsi_ep_eta();
                            Psi_diff_ME             = alexPhiMeson_track->getPsi_diff_ME();
                            Delta_dip_angle         = alexPhiMeson_track->getDelta_dip_angle();
                            qpA                     = alexPhiMeson_track->getqpA();
                            qpB                     = alexPhiMeson_track->getqpB();
                            qpC                     = 0.0;
                            iQxC                    = 0.0;
                            iQyC                    = 0.0;
                            etaC                    = 0.0;

                            InvMass_Fill      = InvAB;
                            InvMass_Fill_miss = InvAB;
                        }

                        if(eAnalysis >= 1 && eAnalysis <= 4) // Xi-/+, Omega-/+ analysis
                        {
                            alexV0_track = alexV0_event->getTrack( i ); // take the track
                            m2A                     = alexV0_track->getm2A();
                            m2B                     = alexV0_track->getm2B();
                            m2C                     = alexV0_track->getm2C();
                            nsA                     = alexV0_track->getnsA()*scale_nSigma_fac;
                            nsB                     = alexV0_track->getnsB()*scale_nSigma_fac;
                            nsC                     = alexV0_track->getnsC()*scale_nSigma_fac;
                            dcaA                    = alexV0_track->getdcaA();
                            dcaB                    = alexV0_track->getdcaB();
                            dcaC                    = alexV0_track->getdcaC();
                            iQxA                    = alexV0_track->getiQxA();
                            iQyA                    = alexV0_track->getiQyA();
                            iQxB                    = alexV0_track->getiQxB();
                            iQyB                    = alexV0_track->getiQyB();
                            iQxC                    = alexV0_track->getiQxC();
                            iQyC                    = alexV0_track->getiQyC();
                            etaA                    = alexV0_track->getetaA();
                            etaB                    = alexV0_track->getetaB();
                            etaC                    = alexV0_track->getetaC();
                            InvAB                   = alexV0_track->getInvAB();
                            InvABC                  = alexV0_track->getInvABC();
                            InvAB_miss              = alexV0_track->getInvAB_miss();
                            InvABC_miss             = alexV0_track->getInvABC_miss();
                            dcaAB                   = alexV0_track->getdcaAB();
                            dcaBC                   = alexV0_track->getdcaBC();
                            Delta_dip_angle         = alexV0_track->getdcaABC();
                            VerdistX                = alexV0_track->getVerdistX();
                            VerdistY                = alexV0_track->getVerdistY();
                            VerdistX2               = alexV0_track->getVerdistX2();
                            VerdistY2               = alexV0_track->getVerdistY2();
                            p_t                     = alexV0_track->getpt();
                            rap                     = alexV0_track->getrap();
                            phi                     = alexV0_track->getphi();
                            theta                   = alexV0_track->gettheta();
                            phi_event_plane         = alexV0_track->getPsi_ep();
                            phi_event_plane_eta_gap = alexV0_track->getPsi_ep_eta();
                            Psi_diff_ME             = alexV0_track->getPsi_diff_ME();
                            scal_prod               = alexV0_track->getscal_prod();
                            scal_prod2              = alexV0_track->getscal_prod2();

                            qpA = TMath::Sqrt(iQxA*iQxA+iQyA*iQyA)*TMath::CosH(etaA);
                            qpB = TMath::Sqrt(iQxB*iQxB+iQyB*iQyB)*TMath::CosH(etaB);
                            qpC = TMath::Sqrt(iQxC*iQxC+iQyC*iQyC)*TMath::CosH(etaC);

                            InvMass_Fill      = InvABC;
                            InvMass_Fill_miss = InvABC_miss;
                        }


                        //---------------------------------------------------------------
                        // Determine pt bin
                        Int_t ept_bin = -1;
                        for(Int_t r = 0; r < npt_bins-1; r++)
                        {
                            if(p_t >= pt_bin_ranges[0])
                            {
                                if(p_t >= pt_bin_ranges[r] && p_t < pt_bin_ranges[r+1])
                                {
                                    ept_bin = r;
                                }
                            }
                            else break;
                        }
                        //---------------------------------------------------------------


                        //----------------------------------------------------------------------------------------------------
                        // Event plane corrections
                        // 0 = Full TPC event plane sub A, phi-weight method
                        // 1 = Full TPC event plane sub B, phi-weight method
                        // 2 = Full TPC event plane all, phi-weight method
                        // 3 = Full TPC event plane sub A, re-centering method
                        // 4 = Full TPC event plane sub B, re-centering method
                        // 5 = Full TPC event plane all, re-centering method
                        // 6 = Eta sub event plane eta plus, phi-weight method
                        // 7 = Eta sub event plane eta neg, phi-weight method
                        // 8 = Eta sub event plane eta plus, re-centering method
                        // 9 = Eta sub event plane eta neg, re-centering method
                        // used for example here: hEP_shift_params



                        // Apply re-centering correction and calculate Psi angles
                        Float_t phi_event_plane_rc             = -100.0;
                        Float_t phi_event_plane_eta_pos_gap_rc = -100.0;
                        Float_t phi_event_plane_eta_neg_gap_rc = -100.0;

                        // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
                        Float_t flag_weight_part_A_EP[3] = {0.0,0.0,0.0}; // [eta pos, eta neg, full TPC]
                        Float_t flag_weight_part_B_EP[3] = {0.0,0.0,0.0}; // [eta pos, eta neg, full TPC]
                        Float_t flag_weight_part_C_EP[3] = {0.0,0.0,0.0}; // [eta pos, eta neg, full TPC]
                        if(fabs(etaA) < eta_EP_cut && dcaA < dcaAB_EP_cut &&
                           fabs(qpA) > MomentumA_EP_low_cut && fabs(qpA) < MomentumA_EP_high_cut)
                        {
                            if(etaA > 0 && fabs(etaA) > eta_gap) flag_weight_part_A_EP[0] = 1.0; // track was used for EP calculation, eta pos
                            if(etaA < 0 && fabs(etaA) > eta_gap) flag_weight_part_A_EP[1] = 1.0; // track was used for EP calculation, eta neg
                            flag_weight_part_A_EP[2] = 1.0;
                        }
                        if(fabs(etaB) < eta_EP_cut && dcaB < dcaAB_EP_cut &&
                           fabs(qpB) > MomentumA_EP_low_cut && fabs(qpB) < MomentumA_EP_high_cut && SE_ME_loop == 0)
                        {
                            if(etaB > 0 && fabs(etaB) > eta_gap) flag_weight_part_B_EP[0] = 1.0; // track was used for EP calculation, eta pos
                            if(etaB < 0 && fabs(etaB) > eta_gap) flag_weight_part_B_EP[1] = 1.0; // track was used for EP calculation, eta neg
                            flag_weight_part_B_EP[2] = 1.0;
                        }
                        if(fabs(etaC) < eta_EP_cut && dcaC < dcaAB_EP_cut &&
                           fabs(qpC) > MomentumA_EP_low_cut && fabs(qpC) < MomentumA_EP_high_cut && SE_ME_loop == 0)
                        {
                            if(etaC > 0 && fabs(etaC) > eta_gap) flag_weight_part_C_EP[0] = 1.0; // track was used for EP calculation, eta pos
                            if(etaC < 0 && fabs(etaC) > eta_gap) flag_weight_part_C_EP[1] = 1.0; // track was used for EP calculation, eta neg
                            flag_weight_part_C_EP[2] = 1.0;
                        }

                        // --> apply auto correlation correction!!!
                        if(N_tracks_EP_Q_vectors[0] > 0.0 && N_tracks_EP_Q_vectors[2] > 0.0)
                        {
                            phi_event_plane_eta_pos_gap_rc = calc_phi_event_plane_2nd(EP_Q_vectors_rc[0]-(flag_weight_part_A_EP[0]*iQxA+flag_weight_part_B_EP[0]*iQxB+flag_weight_part_C_EP[0]*iQxC)/N_tracks_EP_Q_vectors[0],EP_Q_vectors_rc[1]-(flag_weight_part_A_EP[0]*iQyA+flag_weight_part_B_EP[0]*iQyB+flag_weight_part_C_EP[0]*iQyC)/N_tracks_EP_Q_vectors[0]); //
                            phi_event_plane_eta_neg_gap_rc = calc_phi_event_plane_2nd(EP_Q_vectors_rc[2]-(flag_weight_part_A_EP[1]*iQxA+flag_weight_part_B_EP[1]*iQxB+flag_weight_part_C_EP[1]*iQxC)/N_tracks_EP_Q_vectors[2],EP_Q_vectors_rc[3]-(flag_weight_part_A_EP[1]*iQyA+flag_weight_part_B_EP[1]*iQyB+flag_weight_part_C_EP[1]*iQyC)/N_tracks_EP_Q_vectors[2]); //
                        }
                        if(N_tracks_EP_Q_vectors[4] > 0.0 && N_tracks_EP_Q_vectors[5] > 0.0)
                        {
                            phi_event_plane_rc             = calc_phi_event_plane_2nd(EP_Q_vectors_rc[4]-(flag_weight_part_A_EP[2]*iQxA+flag_weight_part_B_EP[2]*iQxB+flag_weight_part_C_EP[2]*iQxC)/N_tracks_EP_Q_vectors[4],EP_Q_vectors_rc[5]-(flag_weight_part_A_EP[2]*iQyA+flag_weight_part_B_EP[2]*iQyB+flag_weight_part_C_EP[2]*iQyC)/N_tracks_EP_Q_vectors[5]); // subtract the Q-vector of the track to avoid auto correlations
                        }

                        Int_t eta_index = 6;
                        Float_t phi_event_plane_eta_gap_rc = phi_event_plane_eta_pos_gap_rc;
                        if(rap > 0.0)
                        {
                            eta_index = 7; // eta particle is positive -> use negative eta sub event plane
                            phi_event_plane_eta_gap_rc = phi_event_plane_eta_neg_gap_rc;
                        }

                        // Determine shift correction
                        Float_t shift_phi_event_plane            = Get_shift_Psi(phi_event_plane,2,hEP_shift_params,file_bin_rc); // Full TPC event plane all, phi-weight method
                        Float_t shift_phi_event_plane_eta_gap    = Get_shift_Psi(phi_event_plane_eta_gap,eta_index,hEP_shift_params,file_bin_rc); // Eta sub event plane eta plus/neg, phi-weight method
                        Float_t shift_phi_event_plane_rc         = Get_shift_Psi(phi_event_plane_rc,5,hEP_shift_params,file_bin_rc); // Full TPC event plane all, re-centering method
                        Float_t shift_phi_event_plane_eta_gap_rc = Get_shift_Psi(phi_event_plane_eta_gap_rc,eta_index+2,hEP_shift_params,file_bin_rc); // Eta sub event plane eta plus/neg, re-centering method

                        //Double_t phi_event_plane_all_shift[nep_methods] = {shift_phi_event_plane,shift_phi_event_plane_eta_gap,shift_phi_event_plane_rc,shift_phi_event_plane_eta_gap_rc};

                        // Apply shift correction
                        phi_event_plane            += shift_phi_event_plane;
                        phi_event_plane_eta_gap    += shift_phi_event_plane_eta_gap;
                        phi_event_plane_rc         += shift_phi_event_plane_rc;
                        phi_event_plane_eta_gap_rc += shift_phi_event_plane_eta_gap_rc;

                        Double_t phi_event_plane_all[nep_methods] = {phi_event_plane,phi_event_plane_eta_gap,phi_event_plane_rc,phi_event_plane_eta_gap_rc};
                        Double_t delta_phi_event_plane_all[nep_methods];
                        Int_t phi_bin_all[nep_methods];

                        Float_t phi_event_plane_use = 0.0;
                        if(eEP_method == 0) // full TPC event plane, phi-weights
                        {
                            phi_event_plane_use = phi_event_plane;
                        }
                        if(eEP_method == 1) // eta sub event plane, phi-weights
                        {
                            phi_event_plane_use = phi_event_plane_eta_gap;
                        }
                        if(eEP_method == 2) // full TPC event plane, re-centering
                        {
                            phi_event_plane_use = phi_event_plane_rc;
                        }
                        if(eEP_method == 3) // eta sub event plane, re-centering
                        {
                            phi_event_plane_use = phi_event_plane_eta_gap_rc;
                        }


                        //Float_t phi_orig = phi;


                        //if(eBeamTimeNum == 0 && Momentum < 0.0) // 7.7 GeV and anti-protons -> use different bins
                        //{
                        //    ept_bin = (Int_t)(ept_bin/2.0);
                        //}

                        // -pi/2..pi/2 -> 0..pi
                        if(phi_event_plane_use < 0.0) phi_event_plane_use += TMath::Pi();

                        // -pi..pi -> 0..pi
                        if(phi < 0.0) phi = phi + TMath::Pi();

                        // -pi..pi, delta_phi_angle: 0..pi
                        Double_t delta_phi_angle = phi - phi_event_plane_use;
                        if(phi >= 0.0 && delta_phi_angle >= 0.0) delta_phi_angle = delta_phi_angle;
                        if(phi >= 0.0 && delta_phi_angle < 0.0)  delta_phi_angle += TMath::Pi();
                        if(phi < 0.0  && delta_phi_angle >= -TMath::Pi()) delta_phi_angle += TMath::Pi();
                        if(phi < 0.0  && delta_phi_angle < -TMath::Pi())  delta_phi_angle += 2.0*TMath::Pi();


                        for(Int_t m = 0; m < nep_methods; m++)
                        {
                            // -pi/2..pi/2 -> 0..pi
                            if(phi_event_plane_all[m] < 0.0) phi_event_plane_all[m] += TMath::Pi();

                            // -pi..pi, delta_phi_angle: 0..pi
                            delta_phi_event_plane_all[m] = phi - phi_event_plane_all[m];
                            if(phi >= 0.0 && delta_phi_event_plane_all[m] >= 0.0) delta_phi_event_plane_all[m] = delta_phi_event_plane_all[m];
                            if(phi >= 0.0 && delta_phi_event_plane_all[m] < 0.0)  delta_phi_event_plane_all[m] += TMath::Pi();
                            if(phi < 0.0  && delta_phi_event_plane_all[m] >= -TMath::Pi()) delta_phi_event_plane_all[m] += TMath::Pi();
                            if(phi < 0.0  && delta_phi_event_plane_all[m] < -TMath::Pi())  delta_phi_event_plane_all[m] += 2.0*TMath::Pi();

                            phi_bin_all[m] = (Int_t)((delta_phi_event_plane_all[m]-phi_start)/delta_phi);

                            if(phi_bin_all[m] == 13) phi_bin_all[m] = 0;
                            if(phi_bin_all[m] == 12) phi_bin_all[m] = 1;
                            if(phi_bin_all[m] == 11) phi_bin_all[m] = 2;
                            if(phi_bin_all[m] == 10) phi_bin_all[m] = 3;
                            if(phi_bin_all[m] == 9)  phi_bin_all[m] = 4;
                            if(phi_bin_all[m] == 8)  phi_bin_all[m] = 5;
                            if(phi_bin_all[m] == 7)  phi_bin_all[m] = 6;
                        }
                        //----------------------------------------------------------------------------------------------------



                        //cout << "InvAB = " << InvAB << ", ept_bin = " << ept_bin << ", rap = " << rap << endl;
                        //----------------------------------------------------------------------------------------------------

                        Float_t dcaA_add[n_cuts]      = {0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.05,-0.1,0.0,0.0,0.0,0.1,0.0,0.0,0.0,-0.05,0.05,-0.05,0.05};
                        Float_t dcaB_add[n_cuts]      = {0.0,-0.1,0.0,0.0,0.0,0.0,0.0,0.1,0.0,0.1,-0.1,0.0,0.0,-0.1,0.0,0.0,0.0,0.0,-0.1,0.1};
                        Float_t dcaC_add[n_cuts]      = {0.0,-0.1,0.0,0.0,0.0,0.0,0.0,0.1,0.0,0.1,-0.1,0.0,0.0,-0.1,0.0,0.0,0.0,0.0,-0.1,0.1};
                        Float_t VerdistY_add[n_cuts]  = {0.0,0.05,0.0,0.0,0.15,0.1,0.05,-0.1,0.0,-0.5,0.5,0.05,0.0,0.0,-0.1,0.0,0.0,0.0,-0.05,0.05};
                        Float_t VerdistX_add[n_cuts]  = {0.0,-0.5,0.0,0.0,0.0,0.0,0.0,-0.5,0.0,-0.5,0.5,0.0,-0.5,0.0,1.5,1.0,0.0,0.0,-0.5,0.5};
                        Float_t VerdistX2_add[n_cuts] = {0.0,-0.5,1.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.5,1.0,0.0,1.0,0.0,0.0,-0.5,0.5};
                        Float_t VerdistY2_add[n_cuts] = {0.0,-0.1,0.0,0.0,0.0,0.0,0.0,0.1,-0.2,0.0,0.0,0.2,0.0,0.0,0.0,0.0,-0.1,0.2,0.1,-0.1};

                        Float_t p_t_add[n_cuts]       = {0.0,0.3,-0.2,0.2,0.2,0.1,-0.2,0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.1,0.1,-0.05,0.05};
                        Float_t p_add[n_cuts]         = {0.0,0.15,-0.2,0.2,-0.2,0.2,-0.2,0.2,0.0,0.0,0.4,-0.1,0.2,-0.2,0.0,0.0,-0.1,0.1,-0.05,0.05};
                        Float_t nsLow_add[n_cuts]     = {0.0,0.15,-0.2,0.2,0.1,-0.1,0.0,0.0,-0.1,-0.2,0.0,0.1,0.0,0.0,-0.4,0.4,0.0,0.0,0.0,0.0};
                        Float_t nsHigh_add[n_cuts]    = {0.0,-0.15,-0.2,0.2,0.0,0.0,0.0,0.0,0.1,0.2,0.0,0.1,0.0,0.0,-0.4,0.4,0.0,0.0,0.0,0.0};
                        Float_t t_add[n_cuts]         = {0.0,0.2,0.5,-0.5,-1.0,1.0,0.0,0.0,0.0,0.0,0.5,-0.5,-0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0};
                        Float_t m_add[n_cuts]         = {0.0,-0.5,-0.5,0.5,0.0,0.0,-1.0,-1.0,1.0,0.0,-0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

                        Double_t K0S_low_cut          = 0.483;
                        Double_t K0S_high_cut         = 0.513;
                        Double_t Lambda_low_cut       = 1.108;
                        Double_t Lambda_high_cut      = 1.122;

                        if(
                           fabs(rap) < 1.0
                           && ept_bin >= 0 && ept_bin < npt_bins
                           && erefMult_bin >= 0 && erefMult_bin < nmult_bins_all
                          )
                        {
                            for(Int_t icut = 0; icut < n_cuts; icut++)
                            {
                                if(
                                   (
                                    //----------------------------------------------------------------------------------------------------------------
                                    eAnalysis == 0 && // phi meson analysis
                                    (
                                     (
                                      (eBeamTimeNum == 0 || eBeamTimeNum == 1) &&
                                      (p_t < 1.2+p_t_add[icut] || (p_t >= 1.2+p_t_add[icut] && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36) ) )) &&
                                      ((nsA < ((11.0+m_add[icut])*fabs(qpA)-(7.86+t_add[icut])) && m2A > 0.16 && m2A < 0.36) || (nsA > ((11.0+m_add[icut])*fabs(qpA)-(7.86+t_add[icut])))) &&
                                      ((nsB < ((11.0+m_add[icut])*fabs(qpB)-(7.86+t_add[icut])) && m2B > 0.16 && m2B < 0.36) || (nsB > ((11.0+m_add[icut])*fabs(qpB)-(7.86+t_add[icut])))) &&
                                      ((m2A < -10 && nsA < 2.5+nsHigh_add[icut] && nsA > -(1.2+nsLow_add[icut])) || (m2A > 0.16 && m2A < 0.36)) &&
                                      ((m2B < -10 && nsB < 2.5+nsHigh_add[icut] && nsB > -(1.2+nsLow_add[icut])) || (m2B > 0.16 && m2B < 0.36)))
                                     ||
                                     (
                                      (eBeamTimeNum >= 2) &&
                                      ((fabs(qpA) <= 0.65+p_add[icut] && m2A < -10) || (m2A > 0 && ((fabs(qpA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(qpA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
                                      ((fabs(qpB) <= 0.65+p_add[icut] && m2B < -10) || (m2B > 0 && ((fabs(qpB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(qpB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
                                      ((p_t+p_t_add[icut]) < 0.8 || ((p_t+p_t_add[icut]) >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
                                      (((m2A < -10 && nsA < 2.5+nsHigh_add[icut] && nsA > -(1.5+nsLow_add[icut])) || (m2A > 0.16 && m2A < 0.36)) &&
                                       ((m2B < -10 && nsB < 2.5+nsHigh_add[icut] && nsB > -(1.5+nsLow_add[icut])) || (m2B > 0.16 && m2B < 0.36))
                                      )
                                     )
                                    )
                                    //----------------------------------------------------------------------------------------------------------------
                                   ) ||
                                   (
                                    //----------------------------------------------------------------------------------------------------------------
                                    (eAnalysis == 1 || eAnalysis == 2) && // Xi-/+ analysis
                                    (
                                     (
                                      dcaAB        < dcaAB_cut &&
                                      dcaBC        < dcaBC_cut &&
                                      InvAB        > lambda_low_cut &&
                                      InvAB        < lambda_high_cut
                                     ) &&
                                     (
                                      ( // All particles have a mass
                                       (m2A > proton_low_cut && m2A < proton_high_cut)
                                       && (m2B > pion_low_cut   && m2B  < pion_high_cut)
                                       && (m2C > pion_low_cut   && m2C  < pion_high_cut)
                                       &&
                                       (
                                        dcaA            > 0.0+dcaA_add[icut]
                                        && dcaB         > 1.0+dcaB_add[icut]
                                        && dcaC         > 0.6+dcaC_add[icut]
                                        && VerdistY     > 0.2+VerdistY_add[icut]
                                        && VerdistY2    < 0.6+VerdistY2_add[icut]
                                        && VerdistX     > 4.0+VerdistX_add[icut]
                                        && VerdistX2    > 3.0+VerdistX2_add[icut]
                                       )
                                      )
                                      ||
                                      ( // No particle has a mass
                                       (m2A < -10 || m2A > 5)
                                       && (m2B < -10 || m2B > 5)
                                       && (m2C < -10 || m2C > 5)
                                       &&
                                       (
                                        dcaA            > 0.2+dcaA_add[icut]
                                        && dcaB         > 1.0+dcaB_add[icut]
                                        && dcaC         > 0.8+dcaC_add[icut]
                                        && VerdistY     > 0.3+VerdistY_add[icut]
                                        && VerdistY2    < 0.6+VerdistY2_add[icut]
                                        && VerdistX     > 4.0+VerdistX_add[icut]
                                        && VerdistX2    > 3.5+VerdistX2_add[icut]
                                       )
                                      )
                                      ||
                                      ( // Only 2nd pion has a mass
                                       (m2A < -10 || m2A > 5)
                                       && (m2B < -10 || m2B > 5)
                                       && (m2C > pion_low_cut && m2C  < pion_high_cut)
                                       &&
                                       (
                                        dcaA            > 0.1+dcaA_add[icut]
                                        && dcaB         > 1.5+dcaB_add[icut]
                                        && dcaC         > 1.2+dcaC_add[icut]
                                        && VerdistY     > 0.1+VerdistY_add[icut]
                                        && VerdistY2    < 0.5+VerdistY2_add[icut]
                                        && VerdistX     > 5.0+VerdistX_add[icut]
                                        && VerdistX2    > 3.5+VerdistX2_add[icut]
                                       )
                                      )
                                      ||
                                      ( // All the rest
                                       dcaA            > 0.2 +dcaA_add[icut]
                                       && dcaB         > 1.25+dcaB_add[icut]
                                       && dcaC         > 1.0 +dcaC_add[icut]
                                       && VerdistY     > 0.3 +VerdistY_add[icut]
                                       && VerdistY2    < 0.5 +VerdistY2_add[icut]
                                       && VerdistX     > 5.0 +VerdistX_add[icut]
                                       && VerdistX2    > 3.0 +VerdistX2_add[icut]
                                      )
                                      //||
                                      //(
                                      // eBeamTimeNum == 0 &&
                                      // (
                                      //  erefMult_bin <= 3
                                      // )
                                      //)
                                     )
                                    )
                                    //----------------------------------------------------------------------------------------------------------------
                                   ) ||
                                   (
                                    //----------------------------------------------------------------------------------------------------------------
                                    (eAnalysis == 3 || eAnalysis == 4) && // Omega-/+ analysis
                                    (
                                     (
                                      dcaAB        < dcaAB_cut &&
                                      dcaBC        < dcaBC_cut &&
                                      InvAB        > lambda_low_cut &&
                                      InvAB        < lambda_high_cut
                                     ) &&
                                     (
                                      ( // All particles have a mass
                                       (m2A > proton_low_cut && m2A < proton_high_cut)
                                       && (m2B > pion_low_cut   && m2B  < pion_high_cut)
                                       && (m2C > kaon_low_cut   && m2C  < kaon_high_cut)
                                       &&
                                       (
                                        dcaA            > 0.2+dcaA_add[icut]
                                        && dcaB         > 1.5+dcaB_add[icut]
                                        && dcaC         > 0.6+dcaC_add[icut]
                                        && VerdistY     > 0.3+VerdistY_add[icut]
                                        && VerdistY2    < 0.4+VerdistY2_add[icut]
                                        && VerdistX     > 4.0+VerdistX_add[icut]
                                        && VerdistX2    > 3.0+VerdistX2_add[icut]
                                       )
                                      )
                                      ||
                                      ( // No particle has a mass
                                       (m2A < -10 || m2A > 5)
                                       && (m2B < -10 || m2B > 5)
                                       && (m2C < -10 || m2C > 5)
                                       &&
                                       (
                                        dcaA            > 0.85+dcaA_add[icut]
                                        && dcaB         > 2.5 +dcaB_add[icut]
                                        && dcaC         > 1.05+dcaC_add[icut]
                                        && VerdistY     > 0.65+VerdistY_add[icut]
                                        && VerdistY2    < 0.25+VerdistY2_add[icut]
                                        && VerdistX     > 8.5 +VerdistX_add[icut]
                                        && VerdistX2    > 4.2 +VerdistX2_add[icut]
                                       )
                                      )
                                      ||
                                      ( // Only kaon has a mass
                                       (m2A < -10 || m2A > 5)
                                       && (m2B < -10 || m2B > 5)
                                       && (m2C > kaon_low_cut && m2C  < kaon_high_cut)
                                       &&
                                       (
                                        dcaA            > 0.3+dcaA_add[icut]
                                        && dcaB         > 1.5+dcaB_add[icut]
                                        && dcaC         > 0.7+dcaC_add[icut]
                                        && VerdistY     > 0.4+VerdistY_add[icut]
                                        && VerdistY2    < 0.4+VerdistY2_add[icut]
                                        && VerdistX     > 5.0+VerdistX_add[icut]
                                        && VerdistX2    > 3.6+VerdistX2_add[icut]
                                       )
                                      )
                                      ||
                                      ( // All the rest
                                       dcaA            > 0.4+dcaA_add[icut]
                                       && dcaB         > 1.5+dcaB_add[icut]
                                       && dcaC         > 1.1+dcaC_add[icut]
                                       && VerdistY     > 0.5+VerdistY_add[icut]
                                       && VerdistY2    < 0.4+VerdistY2_add[icut]
                                       && VerdistX     > 6.5+VerdistX_add[icut]
                                       && VerdistX2    > 3.8+VerdistX2_add[icut]
                                      )
                                     )
                                    )
                                    //----------------------------------------------------------------------------------------------------------------
                                   ) ||
                                   (
                                    //----------------------------------------------------------------------------------------------------------------
                                    (eAnalysis == 5 || eAnalysis == 6) && // Lambda & anti-Lambda analysis
                                    (
                                     ( // no particle has a mass
                                      (m2A < -10 || m2A > 5) && (m2B < -10 || m2B > 5)
                                      && dcaA          > 0.6+dcaA_add[icut]
                                      && dcaB          > 1.7+dcaB_add[icut]
                                      && dcaAB         < 1.0
                                      && VerdistX      > 4.0+VerdistX_add[icut]
                                      && VerdistY      < 0.75+VerdistY_add[icut]
                                     )
                                     ||
                                     ( // only pion has a mass
                                      (m2A < -10 || m2A > 5) && (m2B > (-0.013*qpB+0.017) && m2B < 0.04)
                                      && dcaA          > 0.5+dcaA_add[icut]
                                      && dcaB          > 1.5+dcaB_add[icut]
                                      && dcaAB         < 1.0
                                      && VerdistX      > 3.5+VerdistX_add[icut]
                                      && VerdistY      < 0.75+VerdistY_add[icut]
                                     )
                                     ||
                                     ( // only proton has a mass
                                      (m2B < -10 || m2B > 5) && (m2A > 0.5 && m2A < 1.5)
                                      && dcaA          > 0.15+dcaA_add[icut]
                                      && dcaB          > 0.8+dcaB_add[icut]
                                      && dcaAB         < 1.0
                                      && VerdistX      > 2.5+VerdistX_add[icut]
                                      && VerdistY      < 1.2+VerdistY_add[icut]
                                     )
                                     ||
                                     ( // both particles have a mass
                                      (m2A > 0.5 && m2A < 1.5) && (m2B > (-0.013*qpB+0.017) && m2B < 0.04)
                                      && dcaA          > 0.1+dcaA_add[icut]
                                      && dcaB          > 0.7+dcaB_add[icut]
                                      && dcaAB         < 1.0
                                      && VerdistX      > 2.0+VerdistX_add[icut]
                                      && VerdistY      < 1.3+VerdistY_add[icut]
                                     )
                                    )
                                    //----------------------------------------------------------------------------------------------------------------
                                   ) ||
                                   (
                                    //----------------------------------------------------------------------------------------------------------------
                                    (eAnalysis == 7) && // K0S analysis
                                    (
                                     dcaA             > 0.7
                                     && dcaB          > 0.7
                                     && dcaAB         < 1.0
                                     && VerdistX      > 3.0
                                     && VerdistY      < 0.8
                                     && scal_prod     > 0.0
                                     && ((m2A > -0.5 && m2A < 0.15) || (m2A < -10) || (m2A > 5))
                                     && ((m2B > -0.5 && m2B < 0.15) || (m2B < -10) || (m2B > 5))
                                    )
                                    //----------------------------------------------------------------------------------------------------------------
                                   )
                                  )
                                {
                                    hInvMass_pt_miss_syst[SE_ME_loop*n_cuts + icut][ept_bin]    ->Fill(InvMass_Fill_miss,reweight);

                                    if( // Bat-anti missidentification cuts :)
                                       (
                                        (eAnalysis == 3 || eAnalysis == 4) && // Omega-/+ analysis
                                        //(!(InvMass_Fill_miss > 1.313 && InvMass_Fill_miss < 1.328)))
                                        //(InvMass_Fill_miss > 1.328)
                                        (InvMass_Fill_miss > 1.331)
                                       )
                                       //|| // Omega removal not needed for Xi analyis. This will only reduce the Xi signal
                                       //((eAnalysis == 1 || eAnalysis == 2) && // Xi-/+ analysis
                                       // (!(InvMass_Fill_miss > 1.665 && InvMass_Fill_miss < 1.678)))
                                       ||
                                       (eAnalysis == 0) // phi analysis
                                       ||
                                       (eAnalysis == 1 || eAnalysis == 2)// Xi-/+ analysis
                                       ||
                                       (eAnalysis == 5) // Lambda analysis
                                       ||
                                       (eAnalysis == 6 &&
                                        (
                                         ((m2A < -10 || m2A > 5) && !(InvMass_Fill_miss > K0S_low_cut && InvMass_Fill_miss < K0S_high_cut)) ||
                                         (m2A > 0.5 && m2A < 1.5)
                                        )
                                       ) // anti-Lambda analysis
                                       ||
                                       (eAnalysis == 7 && !(InvMass_Fill_miss > Lambda_low_cut && InvMass_Fill_miss < Lambda_high_cut)) // K0S analysis
                                      )
                                    {

                                        if(icut == 0) accepted_SE_ME_counter[SE_ME_loop]++;
                                        //if(icut == 5 && ept_bin == 4) cout << "SE_ME_loop = " << SE_ME_loop << ", icut = " << icut << ", bin = " << SE_ME_loop*n_cuts + icut << ", InvMass = " << InvMass_Fill << ", ept_bin = " << ept_bin << ", reweight = " << reweight << endl;
                                        hInvMass_pt_syst[SE_ME_loop*n_cuts + icut][ept_bin]      ->Fill(InvMass_Fill,reweight);

                                        hInvMass_mult_pt_seme_syst[4+erefMult_bin][ept_bin][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight);

                                        Int_t mult_bin = 0; // 0-80%
                                        hInvMass_mult_pt_seme_syst[mult_bin][ept_bin][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight);
                                        if(erefMult_bin == 7 || erefMult_bin == 8) // 0-10
                                        {
                                            mult_bin = 1;
                                            hInvMass_mult_pt_seme_syst[mult_bin][ept_bin][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight);
                                        }
                                        if(erefMult_bin == 4 || erefMult_bin == 5 || erefMult_bin == 6) // 10-40
                                        {
                                            mult_bin = 2;
                                            hInvMass_mult_pt_seme_syst[mult_bin][ept_bin][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight);
                                        }
                                        if(erefMult_bin == 0 || erefMult_bin == 1 || erefMult_bin == 2 || erefMult_bin == 3) // 40-80
                                        {
                                            mult_bin = 3;
                                            hInvMass_mult_pt_seme_syst[mult_bin][ept_bin][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight);
                                        }

                                        for(Int_t m_ep = 0; m_ep < nep_methods; m_ep++)
                                        {
                                            if(phi_bin_all[m_ep] >= 0 && phi_bin_all[m_ep] < nphi_bins)
                                            {
                                                //hInvMass_mult_pt_phi_ep_seme[4+erefMult_bin][ept_bin][phi_bin_all[m_ep]][m_ep][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight/Event_Plane_resolution[erefMult_bin][m_ep]);

                                                Int_t mult_bin = 0; // 0-80%
                                                hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][ept_bin][phi_bin_all[m_ep]][m_ep][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight/Event_Plane_resolution[erefMult_bin][m_ep]);
                                                if(erefMult_bin == 7 || erefMult_bin == 8) // 0-10
                                                {
                                                    mult_bin = 1;
                                                    hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][ept_bin][phi_bin_all[m_ep]][m_ep][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight/Event_Plane_resolution[erefMult_bin][m_ep]);
                                                }
                                                if(erefMult_bin == 4 || erefMult_bin == 5 || erefMult_bin == 6) // 10-40
                                                {
                                                    mult_bin = 2;
                                                    hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][ept_bin][phi_bin_all[m_ep]][m_ep][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight/Event_Plane_resolution[erefMult_bin][m_ep]);
                                                }
                                                if(erefMult_bin == 0 || erefMult_bin == 1 || erefMult_bin == 2 || erefMult_bin == 3) // 40-80
                                                {
                                                    mult_bin = 3;
                                                    hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][ept_bin][phi_bin_all[m_ep]][m_ep][SE_ME_loop*n_cuts + icut] ->Fill(InvMass_Fill,reweight/Event_Plane_resolution[erefMult_bin][m_ep]);
                                                }
                                            }

                                        }
                                    } // accepted track/event
                                    //--------------------------------------------------------------------------------
                                }
                            } // end loop over cut combinations
                        }
                        //----------------------------------------------------------------------------------------------------
                    } // end of track loop
                }
            }
            //----------------------------------------------------------------------------------------------------

        } // end of event loop
    } // end of same event <-> mixed event loop

    cout << "" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Number of accepted same event entries: "  << accepted_SE_ME_counter[0] << endl;
    cout << "Number of accepted mixed event entries: " << accepted_SE_ME_counter[1] << endl;
    cout << "---------------------------------------------------------------------" << endl;


}

void StStrangenessAna::finalize()
{

    delete refmultCorr;

    cout << "" << endl;
    cout << "Write output" << endl;
    Outputfile      ->cd();

    if(eAnalysis == 8) // pi0 analysis
    {
        NTuple_Four_leptons->Write();
    }
    else
    {

        TDirectory *cut_dir[n_cuts];
        for(Int_t icut = 0; icut < n_cuts; icut++)
        {
            Outputfile      ->cd();
            HistName = "cut";
            HistName += icut;
            cout << "Create directory: " << HistName.Data() << endl;
            cut_dir[icut] = Outputfile      ->mkdir(HistName.Data());
            cut_dir[icut] ->mkdir("hInvMass_pt_cuts");
            cut_dir[icut] ->cd("hInvMass_pt_cuts");

            for(Int_t SE_ME_loop = 0; SE_ME_loop < 2; SE_ME_loop++) // 0 = same event, 1 = mixed event
            {
                for(Int_t pt_bin = 0; pt_bin < npt_bins; pt_bin++)
                {
                    hInvMass_pt_syst[SE_ME_loop*n_cuts + icut][pt_bin]       ->Write();
                    hInvMass_pt_miss_syst[SE_ME_loop*n_cuts + icut][pt_bin]  ->Write();
                }
            }

            cut_dir[icut] ->cd();
            cut_dir[icut] ->mkdir("hInvMass_mult_pt_phi_ep_seme");
            cut_dir[icut] ->cd("hInvMass_mult_pt_phi_ep_seme");

            for(Int_t mult_bin = 0; mult_bin < nmult_bins; mult_bin++)
            {
                for(Int_t pt_bin = 0; pt_bin < npt_bins; pt_bin++)
                {
                    for(Int_t SE_ME_loop = 0; SE_ME_loop < 2; SE_ME_loop++) // 0 = same event, 1 = mixed event
                    {
                        hInvMass_mult_pt_seme_syst[mult_bin][pt_bin][SE_ME_loop*n_cuts + icut] ->Write();
                        if(mult_bin < nmult_bins4)
                        {
                            for(Int_t phi_bin = 0; phi_bin < nphi_bins; phi_bin++)
                            {
                                for(Int_t m_ep = 0; m_ep < nep_methods; m_ep++)
                                {

                                    hInvMass_mult_pt_phi_ep_seme_syst[mult_bin][pt_bin][phi_bin][m_ep][SE_ME_loop*n_cuts + icut] ->Write();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Outputfile      ->cd();
    cout << "Close output file" << endl;
    cout << "" << endl;
    Outputfile      ->Close();

}

