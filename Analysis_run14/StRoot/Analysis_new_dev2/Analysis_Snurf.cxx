#include "Analysis_Snurf.h"


static Int_t Lambda_ME_counter = 0;
static Int_t eEvent_numA = 0;
static Int_t eEvent_numB = 0;
static Int_t total_counter = 0;
static Long64_t n_total_entries = 0;


static const Bool_t debug_flag = kFALSE;

static const Double_t MAGFIELDFACTOR = kilogauss;

//*********************************************************************************************************
// Beam time dependend cuts
static Double_t vertex_z_cut  = 0.0;
static Int_t fBeamTimeNum;
//*********************************************************************************************************



//*********************************************************************************************************
// Helix parameters
static const Double_t curv_to_invpt_ratio = 0.768572/514.164;
//*********************************************************************************************************



//******************* Important constants *****************************************************************
static const Float_t clight = 29.9792458;//[cm/ns]
static const Int_t N_mass = 16;
static const Double_t mass_array[N_mass]     = {0.0,0.0,0.00051099892,0.00051099892,0.0,0.105658369,0.105658369,0.1349766,0.13957018,0.13957018,0.497648,0.493677,0.493677,0.93956536,0.93827203,0.93827203}; // in GeV/c^2  {empty,gamma,e+,e-,neutrino,mu+,mu-,pi0,pi+,pi-,K0long,K+,K-,n,p,anti-p}
static const Double_t Mass2_low_cut[N_mass]  = {-1.0,-1.0,-0.2,-0.2,-1.0,-1.0,-1.0,-1.0,-0.2,-0.2,-1.0,0.125,0.125,-1.0,0.5,0.5};
static const Double_t Mass2_high_cut[N_mass] = {1.0,1.0,0.1,0.1,1.0,1.0,1.0,1.0,0.1,0.1,1.0,0.36,0.36,1.0,1.4,1.4};
static const Double_t Mass2_low_cut_strict[N_mass]  = {-1.0,-1.0,-0.2,-0.2,-1.0,-1.0,-1.0,-1.0,-0.03,-0.03,-1.0,0.125,0.125,-1.0,0.7,0.7};
static const Double_t Mass2_high_cut_strict[N_mass] = {1.0,1.0,0.1,0.1,1.0,1.0,1.0,1.0,0.05,0.05,1.0,0.36,0.36,1.0,1.1,1.1};
//*********************************************************************************************************


//-------------------------------------------------------------
// Parameters for re-centering and shift correction
static TFile *re_centering; // re-centering correction file
static TFile *shift_method; // shift method correction file
static TH1F* h_runId_index_rc; // runId <-> index for re-centering
static TH1F* h_runId_index_rc_inv; // index for re-centering <-> runId
static Int_t h_Nbins_rc_B;
static const Int_t n_params = 3;
static const Int_t n_z_vertex_bins = 10;
static const Int_t nPsi_EP = 10;
static const Int_t n_shift_par = 5;
static TH1F* h_rc_QxQy_etapm_z_vs_index_[6][n_z_vertex_bins][n_params];
static TH1F* hEP_shift_params[nPsi_EP][n_shift_par];
//-------------------------------------------------------------

StV0TofCorrection* V0_tof_corr;



//-------------------------------------------------------------
static TH1F* hBBC_ADC_tiles_Eeast_West[2][24]; // [East, West][tiles]
static TProfile* pBBC_ADC_tiles_vs_RefMult_Eeast_West[2][24]; // [East, West][tiles]
static TH2F* hBBC_ADC_tiles_vs_RefMult_Eeast_West[2][24]; // [East, West][tiles]
static TH2F* hrefMult_BBC_hits_East_West[2];   // [East, West]
static TH1F* heta_EP;
static TH2F* hFTPC_BBC_corr;
static TH1F* h_FTPC_phi;
//-------------------------------------------------------------



static StRefMultCorr* refmultCorrUtil;
static Int_t erefMult_bin;
static Int_t erefMult_bin16;

static Double_t nsigma_scaling_fac = 1.0; // nsigma is wrongly calibrated for 27 GeV -> a factor ~1.9 too small at 27 GeV

// mixed event -> STAR

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
static const Int_t n_harmonics             = 3; // number of harmonics used
static const Int_t n_eta_gap_values        = 6; // number of eta gap values used
static const Double_t eta_gap_array[n_eta_gap_values] = {0.05,0.1,0.2,0.3,0.4,0.5}; // eta gap values

static Int_t vertex_counter = 0;

// From random sub event -> event plane analysis
static TRandom ran_seed_number;
static Int_t seed_number;

static const Int_t N_Ana                = 4;  // number of different analysis
static const Int_t N_max_PIDs           = 16; // maximum number of PIDs
static const Int_t N_max_PIDs_per_track = 4;  // maximum number of simultaneous PIDs for one track
static const Int_t N_refmult_bins       = 9;
static const Int_t N_z_axis_bins        = 14;
static const Int_t N_EP_bins            = 10;

static const Double_t phi_EP_start    = -TMath::Pi()/2.0;
static const Double_t phi_EP_stop     = TMath::Pi()/2.0;
static const Double_t phi_EP_step     = (phi_EP_stop-phi_EP_start)/((Double_t)N_EP_bins);

static const Int_t N_max_sample         = 7; // 7 maximum number of events[z][refmult] until mixing is started
static Int_t N_max_sample_array[N_refmult_bins][N_z_axis_bins][N_EP_bins];


static const Int_t N_max_tracks         = 5000; // maximum number of tracks per event
static Int_t Sample_Counter[N_refmult_bins][N_z_axis_bins][N_EP_bins]; // counts the number of events for each bin
static Int_t Sample_Event_Num[N_refmult_bins][N_z_axis_bins][N_EP_bins][N_max_sample]; // event number for each bin
static const Int_t N_Beamtime = 9; // 0 == 7.7, 1 == 11.5, 2 == 39, 3 == 62.4, 4 == 19.6, 5 == 27, 6 == 200, 7 == 14.5, 8 = 200 run14
static Double_t Z_axis_table[N_Beamtime][N_z_axis_bins];
static Int_t PIDs_track[N_Ana][N_max_PIDs_per_track];

static Double_t non_prim_to_prim_ratio = 0.0;
static Double_t tracks_left_to_right_ratio = 0.0;
static Double_t tracks_left_to_right_diff = 0.0;

static Int_t good_omega_dca_counter = 0;
static Int_t bad_omega_dca_counter  = 0;

static Int_t fAnalysisNum = 0;
static const Int_t    np_bins        = 15;
static TCutG* cut_Kaon[np_bins];
static Float_t p_bin_ranges[np_bins+1] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.5,4.0,5.0};


//------------------------------------------------------------------------------------------------------------------------------------
// bad run lists
static const Int_t n_bad_run_numbers[N_Beamtime] = {328,27,38,105,35,34,179,1,1};
static const Int_t bad_run_list_7GeV[328]   = {11114084,11114085,11114086,11114088,11114089,11114094,11114095,11114100,11114109,11115005,11115007,11115013,11115019,11115025,11115027,11115028,11115030,11115032,11115051,11115062,11115064,11115069,11115072,11115078,11115079,11115080,11115086,11115088,11115094,11116001,11116002,11116005,11116006,11116010,11116014,11116020,11116023,11116028,11116060,11116061,11116062,11116064,11116068,11116070,11116072,11116073,11116075,11117002,11117006,11117031,11117033,11117034,11117036,11117039,11117044,11117045,11117046,11117052,11117055,11117063,11117064,11117071,11117075,11117085,11117088,11117089,11117090,11117093,11117094,11117095,11117098,11117100,11117103,11117104,11117107,11118007,11118008,11118016,11118024,11118025,11118026,11118039,11118044,11119001,11119003,11119006,11119007,11119009,11119012,11119013,11119015,11119016,11119017,11119022,11119024,11119026,11119029,11119030,11119056,11119057,11119060,11119062,11119067,11119069,11119070,11119071,11119074,11119075,11119077,11119079,11119081,11119090,11119091,11119100,11119101,11120003,11120006,11120008,11120011,11120014,11120019,11120023,11120030,11120034,11120037,11120039,11120040,11120045,11120052,11120057,11120062,11120063,11120069,11120070,11120071,11120074,11120077,11120078,11120084,11120092,11121006,11121014,11121015,11121019,11121029,11121030,11121034,11121035,11121043,11121044,11121054,11121058,11121066,11121067,11121070,11121075,11121082,11122001,11122007,11122008,11122010,11122017,11122024,11122037,11122038,11122047,11122048,11122049,11122050,11122053,11122058,11122062,11122069,11122073,11122078,11122085,11122097,11123003,11123004,11123015,11123026,11123028,11123040,11123044,11123055,11123057,11123058,11123059,11123067,11123075,11123076,11123077,11123079,11123081,11123084,11123086,11123088,11123089,11123093,11123094,11123095,11123100,11123101,11123102,11123104,11124001,11124005,11124007,11124008,11124015,11124016,11124018,11124041,11124046,11124050,11124051,11124052,11124053,11124058,11124060,11124061,11124062,11124063,11124064,11124065,11124066,11124069,11125002,11125003,11125004,11125005,11125006,11125008,11125012,11125013,11125014,11125015,11125016,11125017,11125020,11125021,11125022,11125023,11125073,11125081,11125089,11125090,11125096,11125097,11126005,11126006,11126007,11126016,11126018,11126022,11126023,11127001,11127002,11127043,11128005,11128012,11128018,11128050,11128056,11128072,11129018,11129022,11129028,11129051,11130027,11130034,11130057,11131038,11131062,11132013,11132070,11133006,11133019,11134053,11134060,11134067,11134076,11135068,11136003,11136005,11136006,11136007,11136008,11136012,11136013,11136014,11136061,11136076,11136101,11136130,11136160,11136163,11137019,11138027,11138049,11138086,11138124,11139014,11140076,11140086,11141063,11142117,11143026,11143028,11144001,11144009,11144031,11144033,11144040,11144043,11144052,11145008,11145028,11145035,11146061,11146076,11146079,11147004,11147006,11147014,11147017,11147021,11147023};
static const Int_t bad_run_list_11GeV[27]   = {11148039,11148045,11149001,11149008,11149010,11149011,11149015,11149047,11150016,11150025,11150028,11151036,11151040,11151050,11152016,11152036,11152078,11153032,11153042,11155001,11155009,11156003,11156009,11157012,11158006,11158022,11158024};
static const Int_t bad_run_list_19GeV[35]   = {12113091,12114007,12114035,12114078,12114092,12114116,12115009,12115014,12115015,12115016,12115018,12115019,12115020,12115022,12115023,12115062,12115073,12115093,12115094,12116012,12116054,12117010,12117016,12117020,12117065,12119040,12119042,12120017,12120026,12121017,12121022,12121034,12121050,12121067,12122019};
static const Int_t bad_run_list_27GeV[34]   = {12172050,12172051,12172055,12173030,12173031,12173032,12173033,12173034,12174067,12174085,12175062,12175087,12175113,12175114,12175115,12176001,12176044,12176054,12176071,12177015,12177061,12177092,12177099,12177101,12177106,12177107,12177108,12178003,12178004,12178005,12178006,12178013,12178099,12178120};
static const Int_t bad_run_list_39GeV[38]   = {11199124,11100002,11100045,11101046,11102012,11102051,11102052,11102053,11102054,11102055,11102058,11103035,11103056,11103058,11103092,11103093,11105052,11105053,11105054,11105055,11107007,11107042,11107057,11107061,11107065,11107074,11108101,11109013,11109077,11109088,11109090,11109127,11110013,11110034,11110073,11110076,11111084,11111085};
static const Int_t bad_run_list_62GeV[105]  = {11080072,11081023,11081025,11082012,11082013,11082046,11082056,11082057,11084009,11084011,11084012,11084013,11084020,11084021,11084035,11084044,11084064,11085015,11085025,11085030,11085046,11085055,11085056,11085057,11086005,11086007,11087001,11087002,11087003,11087004,11088013,11089026,11089028,11089029,11089055,11089068,11089072,11091007,11091015,11091021,11091078,11092010,11092011,11092012,11092032,11092033,11092034,11092067,11092096,11093001,11094016,11094017,11094018,11094019,11094020,11094021,11094022,11094023,11094024,11094027,11094028,11094042,11094044,11094045,11094046,11094047,11094048,11094050,11094051,11094052,11094053,11094054,11094055,11094074,11094075,11094077,11095001,11095002,11095003,11095004,11095005,11095006,11095009,11095010,11095011,11095012,11095013,11095014,11095015,11095022,11095040,11095048,11095050,11095051,11095061,11095062,11095063,11095064,11095082,11095087,11096024,11096039,11096043,11096044,11097093};
static const Int_t bad_run_list_200GeV[179] = {12126101,12127003,12127017,12127018,12127032,12128025,12132043,12133018,12134023,12136005,12136006,12136014,12136017,12136022,12136023,12136024,12136025,12136027,12136028,12136029,12136030,12136031,12136034,12136054,12138017,12138021,12138081,12138082,12139006,12139007,12139015,12139016,12139028,12139059,12139075,12139076,12139077,12139078,12139079,12139080,12140010,12140011,12140012,12140013,12140014,12140015,12140016,12140018,12140019,12140020,12140021,12140025,12140026,12140027,12140028,12140029,12140042,12140051,12140052,12140053,12140054,12140055,12140056,12140064,12140066,12140067,12141001,12141002,12141003,12141004,12141005,12141006,12141009,12141014,12141015,12141016,12141017,12141018,12141019,12141026,12141027,12141028,12141029,12141030,12141032,12141033,12141034,12141035,12141036,12141041,12141042,12141043,12141044,12141045,12141046,12141048,12141050,12141051,12141052,12141056,12141059,12141060,12141061,12141062,12141063,12141064,12141065,12141066,12141067,12141071,12141072,12142001,12142002,12142003,12142006,12142013,12142014,12142015,12142016,12142017,12142018,12142019,12142020,12142021,12142022,12142023,12142026,12142027,12142032,12142033,12142034,12142046,12142047,12142048,12142049,12142050,12142051,12142061,12142062,12142063,12142076,12142077,12143016,12143018,12143054,12143075,12144001,12144002,12144013,12144014,12144027,12144028,12157038,12157051,12158040,12158041,12158054,12158056,12158057,12162055,12162056,12162057,12162058,12164037,12164078,12164079,12166002,12166003,12167015,12167024,12167052,12168002,12168009,12168022,12168077,12170044,12170045,12170054,12170056};
static const Int_t bad_run_list_15GeV[1]   = {1};
static const Int_t bad_run_list_200GeV_run14[1]   = {1};
//------------------------------------------------------------------------------------------------------------------------------------



//*****************************************************************************
// Vertexing
static Int_t N_primaries_new  = 0;
static Int_t N_primaries_new2 = 0;
static Int_t N_primaries_old  = 0;
static Int_t N_tofmatch_new  = 0;
static Int_t N_tofmatch_new2 = 0;
static Int_t N_tofmatch_old  = 0;
static Double_t vector_dist_prim_second = 0.0;
//*****************************************************************************

static StAlexEvent  alex_event;
static StAlexEvent *alex_event_ptr;
static StAlexTrack *alex_track;
static TTree       *Tree_hadron_v2;
static const char ALEX_EVENT_TREE[]   = "AlexEvents";
static const char ALEX_EVENT_BRANCH[] = "Events";

static StAlexV0Event  alexV0_event;
static StAlexV0Event *alexV0_event_ptr;
static StAlexV0Track *alexV0_track;
static TTree         *Tree_V0_v2;
static const char ALEXV0_EVENT_TREE[]   = "AlexV0Events";
static const char ALEXV0_EVENT_BRANCH[] = "Events";

// Group A is for Omega- or Omega+
static StAlexV0Event  alexV0_event_A;
static StAlexV0Event *alexV0_event_ptr_A;
static StAlexV0Track *alexV0_track_A;

// Group B is for Xi- or Xi+
static StAlexV0Event  alexV0_event_B;
static StAlexV0Event *alexV0_event_ptr_B;
static StAlexV0Track *alexV0_track_B;

static TTree         *Tree_Sigma0_v2;
static TTree         *Tree_OmegaPV0_v2;
static TTree         *Tree_OmegaMV0_v2;
static TTree         *Tree_XiPV0_v2;
static TTree         *Tree_XiMV0_v2;
static TTree         *Tree_KaonV0_v2;

static StAlexPhiMesonEvent  alexPhiMeson_event;
static StAlexPhiMesonEvent *alexPhiMeson_event_ptr;
static StAlexPhiMesonTrack *alexPhiMeson_track;
static TTree         *Tree_PhiMeson_v2;
static const char ALEXPHIMESON_EVENT_TREE[]   = "AlexPhiMesonEvents";
static const char ALEXPHIMESON_EVENT_BRANCH[] = "Events";

static StD0Event  D0_event;
static StD0Event *D0_event_ptr;
static StD0Track *D0_track;
static TTree     *Tree_D0_v2;
static const char D0_EVENT_TREE[]   = "D0Events";
static const char D0_EVENT_BRANCH[] = "Events";



static const Int_t n_poly_marker_track   = 5000;
static const Int_t n_cEventDisplay_array = 10;
static TCanvas* cEventDisplay_array[n_cEventDisplay_array];
static TCanvas* cEventDCA_array[n_cEventDisplay_array];
static TCanvas* cEvent_MeanX_array[n_cEventDisplay_array];
static TCanvas* cEvent_MeanY_array[n_cEventDisplay_array];
static TCanvas* cEvent_MeanZ_array[n_cEventDisplay_array];
static TPolyLine3D   *pTrack[n_poly_marker_track];
static TPolyLine3D   *pTrack_extrapolate[n_poly_marker_track];
static TPolyLine3D   *pTrack_extrapolate_outside[n_poly_marker_track];
static TPolyLine3D   *pTrack_extrapolate_outsideB[n_poly_marker_track];
static TPolyLine3D   *pTrack_extrapolate_to_Tof[n_poly_marker_track];
static TPolyLine3D   *XAxis;
static TPolyLine3D   *YAxis;
static TPolyLine3D   *ZAxis;
static TPolyLine3D   *BeamLine;
static TMarker3DBox  *pPrimaryVertex;
static TMarker3DBox  *pPrimaryVertex2;
static TMarker3DBox  *pPrimaryVertex_old;
static TMarker3DBox  *pTofHit[n_poly_marker_track];
static TPolyLine3D   *TPC_endcaps[4];
static TH1F          *hEventDCA_new;
static TH1F          *hEventDCA_old;
static TH1F          *hEvent_fit_points[n_cEventDisplay_array];
static TH1F          *hEvent_pt[n_cEventDisplay_array];
static TH1F          *hEvent_p[n_cEventDisplay_array];
static TH1F          *hEvent_dca[n_cEventDisplay_array];
static TH1F          *hEvent_MeanX;
static TH1F          *hEvent_MeanY;
static TH1F          *hEvent_MeanZ;
static TH1D          *hnon_prim_to_prim_ratio;
static TH1D          *htracks_left_to_right_ratio;
static TH1D          *htracks_left_to_right_diff;
static TH1F          *hrho_pippim_inv_mass;

static TH1D          *h_dca_diff;
static TH1D          *h_decay_length_diff;

static TH1D          *h_phi_pp_corr[2][3];
static TH1D          *h_cos_phi_pp_corr[2][3];

const Int_t n_pp_harm = 4;
static TProfile2D    *p_pp_corr[N_refmult_bins][2][n_pp_harm][3]; // [centrality bin][w/o weights, with weights][harmonic][p-p,pbar-pbar,p-pbar]

// Phi correction histograms
static const Int_t nPhi_corr_days     = 365; // 365
static const Int_t nPhi_corr_z_bins   = 5;   // 5
static const Int_t nPhi_corr_eta_bins = 6;   // 6
static const Int_t nPhi_corr_pt_bins  = 4;   // 4
static TH1F* hPhi_corr[nPhi_corr_days][nPhi_corr_z_bins][nPhi_corr_eta_bins][2][nPhi_corr_pt_bins];  // [2] for different polarities, 0 = +, 1 = -
static TH1F* hPhi_corr_in[nPhi_corr_days][nPhi_corr_z_bins][nPhi_corr_eta_bins][2][nPhi_corr_pt_bins];
static Double_t hPhi_corr_in_max_entry[nPhi_corr_days][nPhi_corr_z_bins][nPhi_corr_eta_bins][2][nPhi_corr_pt_bins];
static Int_t hPhi_days_use[nPhi_corr_days];
static Int_t hPhi_days_in[nPhi_corr_days];

static const Float_t vertex_z_array[N_Beamtime] = {70,50,40,40,70,70,40,70,13}; // 7.7, 11.5, 39, 62.4, 19.6, 27, 200, 14.5, 200 run14

static const Int_t n_pt_bins_pp_corr = 5;
static Double_t pp_corr_pt_bins[2][n_pt_bins_pp_corr] =
{
    {0.2,0.6,1.0,1.5,2.0},
    {0.6,1.0,1.5,2.0,3.0}
};


static Float_t phi_corr_z_start;
static Float_t phi_corr_z_stop;
static Float_t phi_corr_eta_start;
static Float_t phi_corr_eta_stop;
static Float_t phi_corr_pt_start;
static Float_t phi_corr_pt_stop;
static Float_t phi_corr_delta_z;
static Float_t phi_corr_delta_eta;
static Float_t phi_corr_delta_pt;


static const Int_t nMinDCA = 20;
static TH1F* hMinDCA[nMinDCA];
static TH1F* hMinDCA_X[nMinDCA];
static TH1F* hMinDCA_origX[nMinDCA];
static TCanvas* cMinDCA;
static TCanvas* cEventDisplay;
static TCanvas* cVertex_distr;
static TH1F* heta;
static TH1F* hphi_eta_pos;
static TH1F* hphi_eta_neg;
static TH1F* hphi;
static TH1F* hsign;
static TH2F* hTPC_dEdx_vs_p;
static TH2F* hTOF_vs_p;
static TH1F* hMass;
static TH1F* hx_vertex_distr;
static TH1F* hy_vertex_distr;
static TH1F* hz_vertex_distr;
static TH2F* h_nhitsFit;
static TH2F* h_eta_phi[10];
static TH2F* h_nhitsFit_TOF;

//**************************************************************************************
// K*(892)0 analysis
static TH1F* hKStar_inv_mass;
static TH1F* hKStar_inv_massB;
static const    Int_t    nmult_bins   = 1;
static const    Int_t    npt_bins     = 7;
static const    Int_t    npt_bins_D0         = 110;
static const    Int_t    npt_bins_Sigma1385  = 45;
static const    Int_t    nphi_bins    = 12;
static const    Double_t phi_start    = 0.0;
static const    Double_t phi_stop     = TMath::Pi();
static Double_t delta_phi             = (phi_stop-phi_start)/((Double_t)nphi_bins);
static const Int_t    n_inv_mass_bins  = 8000;
static const Double_t start_inv_mass  = -10.0;
static const Double_t stop_inv_mass   = 10.0;  // maybe something funny at high masses :)

static TH1D* h_inv_mass_mult_pt_phi[nmult_bins][npt_bins][nphi_bins];
static TH1D* h_inv_mass_D0_pt[npt_bins_D0];
static TH1D* h_inv_mass_Sigma1385_pt[npt_bins_Sigma1385];
static Int_t phi_bin_counter[nphi_bins];
static TH1D* hDelta_phi_KStar;
static TH1D* hEP_Psi_KStar;
static TH1D* hPhi_KStar;
static TH2D* hPhi_vs_Psi;
static TH2D* hTheta_vs_phi;
//**************************************************************************************


//**************************************************************************************
// rho0 analysis
static TH1F* pipi_inv_mass;
static const Double_t start_inv_mass_pipi = 0.0;
static const Double_t stop_inv_mass_pipi  = 20.0;
static const Int_t N_pipi_bins            = 10000;
static const Int_t N_pipi_spectra         = 6;
static const Int_t N_pipi_pt_bins         = 10;

static TH1D* h_inv_mass_pipi[N_pipi_spectra][N_pipi_pt_bins];
//**************************************************************************************

static Int_t ED_NT_counter = 0;
static Int_t good_event  = 0;
static Int_t total_event = 0;
static Int_t bad_event   = 0;
static Int_t event_counter = 0;
static Int_t ED_counter = 0;
static Int_t n_poly_track_counter = 0;
static Int_t comb_counter_global = 0;
static Float_t x_mean = -999.0; // new vertex position
static Float_t y_mean = -999.0;
static Float_t z_mean = -999.0;
static Float_t x_meanA = -999.0; // new vertex position
static Float_t y_meanA = -999.0;
static Float_t z_meanA = -999.0;
static Float_t x_meanB = -999.0; // new vertex position
static Float_t y_meanB = -999.0;
static Float_t z_meanB = -999.0;
static Float_t x_mean2 = -999.0; // new vertex position
static Float_t y_mean2 = -999.0;
static Float_t z_mean2 = -999.0;
static Float_t EP_Qx   = 0.0; // Event plane vector
static Float_t EP_Qy   = 0.0;
static Float_t EP_Qx_B = 0.0; // Event plane vector for second event -> event mixing
static Float_t EP_Qy_B = 0.0;
static Float_t EP_Qx_eta_pos = 0.0;   // sub event vectors for event plane resolution with pt and phi weight
static Float_t EP_Qy_eta_pos = 0.0;
static Float_t EP_Qx_eta_neg = 0.0;
static Float_t EP_Qy_eta_neg = 0.0;
static Float_t EP_Qx_eta_pos_ptw = 0.0;   // sub event vectors for event plane resolution only with pt weight
static Float_t EP_Qy_eta_pos_ptw = 0.0;
static Float_t EP_Qx_eta_neg_ptw = 0.0;
static Float_t EP_Qy_eta_neg_ptw = 0.0;
static Float_t EP_Qx1_eta_pos_ptw = 0.0;   // sub event vectors for event plane resolution only with pt weight
static Float_t EP_Qy1_eta_pos_ptw = 0.0;
static Float_t EP_Qx1_eta_neg_ptw = 0.0;
static Float_t EP_Qy1_eta_neg_ptw = 0.0;
static Float_t EP_Qx_ptw         = 0.0;   // event plane vector only with pt weight
static Float_t EP_Qy_ptw         = 0.0;
static Float_t EP_Qx_eta_pos_B   = 0.0;   // sub event vectors for event plane resolution
static Float_t EP_Qy_eta_pos_B   = 0.0;
static Float_t EP_Qx_eta_neg_B   = 0.0;
static Float_t EP_Qy_eta_neg_B   = 0.0;
static Float_t EP_Qx_subA_ptw    = 0.0;   // sub event vectors for event plane resolution
static Float_t EP_Qy_subA_ptw    = 0.0;
static Float_t EP_Qx_subB_ptw    = 0.0;
static Float_t EP_Qy_subB_ptw    = 0.0;
static Int_t Qtracks_used_eta_pos = 0;
static Int_t Qtracks_used_eta_neg = 0;
static Int_t Qtracks_used         = 0;

static Int_t n_primaries     = 0;
static Int_t n_non_primaries = 0;
static Int_t n_tofmatch_prim = 0;

// Ntuple
static TNtuple *PhiKPKM_NT;
static Float_t PhiKPKM_NTDataArray[19];
static TNtuple *PhiKPKM_V0_NT;
static Float_t PhiKPKM_V0_NTDataArray[33];
static TNtuple *D0_NT;
static Float_t D0_NTDataArray[19];
static TNtuple *Phi_K0SK0S_NT;
static Float_t Phi_K0SK0S_NTDataArray[5];
static TNtuple *DiLepton_NT;
static Float_t DiLepton_NTDataArray[15];
static TNtuple *K0SPiPPiM_NT;
static Float_t K0SPiPPiM_NTDataArray[15];
static TNtuple *K0LPiPPiM_NT;
static Float_t K0LPiPPiM_NTDataArray[15];
static TNtuple *Lambda_X_NT;
static Float_t Lambda_X_NTDataArray[46];
static TNtuple *ThetaPlus_NT;
static Float_t ThetaPlus_NTDataArray[25];
static TNtuple *K0S_pippim_NT;
static Float_t K0S_pippim_NTDataArray[28];
static TNtuple *Xi_NT;
static Float_t Xi_NTDataArray[42];
static TNtuple *Omega_NT;
static Float_t Omega_NTDataArray[43];
static TNtuple *Kaon_NT;
static Float_t Kaon_NTDataArray[17];
static TNtuple *KStarPM_NT;
static Float_t KStarPM_NTDataArray[25];
static TNtuple *D0KMPiP_NT;
static Float_t D0KMPiP_NTDataArray[15];
static TNtuple *Xi1530_NT;
static Float_t Xi1530_NTDataArray[29];
static TNtuple *Vertex_NT;
static Float_t Vertex_NTDataArray[34];
static TNtuple *EDisplay_NT;
static Float_t EDisplay_NTDataArray[41];
static TNtuple *EventPlane_NT;
static Float_t EventPlane_NTDataArray[50];
static TNtuple *EventAnalysis_NT;
static Float_t EventAnalysis_NTDataArray[13];
static TNtuple *EventPlane_array_NT[n_harmonics];
static const Int_t n_variables_event_plane_harm = (24 + 10*n_eta_gap_values*3 + 6);
static Float_t EventPlane_harm_NTDataArray[n_variables_event_plane_harm];
//static TNtuple *Hadron_v2_NT;
//static Float_t Hadron_v2_NTDataArray[33];
static TNtuple *Phi_Corr_NT;
static Float_t Phi_Corr_NTDataArray[11];
static TNtuple *Hadron_NT;
static Float_t Hadron_NTDataArray[12];
static TNtuple *RefMult_QA_NT;
static Float_t RefMult_QA_NTDataArray[38];

static TNtuple *Helix_NT;
static Float_t Helix_NTDataArray[12];

static TNtuple *Vertex_lin_NT;
static Float_t Vertex_lin_NTDataArray[13];




static const Int_t N_d4s_Lambda = 300;
static TNtuple *d4s_NT;
static Float_t d4s_NTDataArray[24];
static StPhysicalHelixD d4s_Lambda_helix[2][N_d4s_Lambda];
static StThreeVectorF d4s_Lambda_vertex[2][N_d4s_Lambda];
static TLorentzVector d4s_Lambda_TLV[3][N_d4s_Lambda];
static Float_t d4s_Lambda_daughters[18][N_d4s_Lambda];
static Float_t d4s_Lambda_properties[4][N_d4s_Lambda];
static TH1* h_d4s_InvMassAB[4];

static StPhysicalHelixD d4s_Xi_helix[2][N_d4s_Lambda];
static StThreeVectorF d4s_Xi_vertex[N_d4s_Lambda];
static TLorentzVector d4s_Xi_TLV[2][N_d4s_Lambda];
static Float_t d4s_Xi_daughters[9][N_d4s_Lambda];
static Float_t d4s_Xi_properties[4][N_d4s_Lambda];



//--------------------------------------------------------------
// Only for test purposes
static StThreeVectorF d4s_Lambda_vertex_old[2][N_d4s_Lambda];
static TLorentzVector d4s_Lambda_TLV_old[2][N_d4s_Lambda];
static Float_t d4s_Lambda_daughters_old[2][10][N_d4s_Lambda];
static Float_t d4s_Lambda_properties_old[2][3][N_d4s_Lambda];

static StThreeVectorF d4s_Xi_vertex_old[N_d4s_Lambda];
static TLorentzVector d4s_Xi_TLV_old[N_d4s_Lambda];
static StPhysicalHelixD d4s_Xi_helix_old[N_d4s_Lambda];
static Float_t d4s_Xi_daughters_old[5][N_d4s_Lambda];
static Float_t d4s_Xi_properties_old[4][N_d4s_Lambda];
//--------------------------------------------------------------



static Std4sEvent    d4s_event;
static Std4sMother  *d4s_Mother;
static Std4sEvent   *d4s_event_ptr;
static Std4sXi      *d4s_Xi;
static Std4sLambda  *d4s_Lambda;
static Std4sTrack   *d4s_Track;
static TTree        *Tree_d4s;
static const char D4S_EVENT_TREE[]   = "d4s_Events";
static const char D4S_EVENT_BRANCH[] = "Events";

static StJetTrackEvent JetTrackEvent;

static StJetTrackEvent *JetTrackEvent_ptr;
static StJetTrackParticle  *JetTrackParticle;
static TTree        *Tree_JetTrackEvent;
static const char JETTRACK_EVENT_TREE[]   = "JetTrackEvent";
static const char JETTRACK_EVENT_BRANCH[] = "Events";

static StLambdaEvent    Lambda_event;
static StLambdaEvent   *Lambda_event_ptr;
static StLambdaLambda  *Lambda_Lambda;
static TTree           *Tree_Lambda;
static const char LAMBDA_EVENT_TREE[]   = "Lambda_Events";
static const char LAMBDA_EVENT_BRANCH[] = "Events";



static TCutG* KaonP_dEdx_cut;


static TFile *Outputfile;

Long64_t dummy_counter      = 0;
Long64_t dummy_counter_array[9];
Long64_t dummy_counter_loop = 0;
Long64_t dummy_counter_A    = 0;
time_t start_time_ana,end_time_ana;
double diff_time_ana;

static StPicoAlexEvent* event_A_ana;
static StPicoAlexEvent* event_B_ana;
static StPicoAlexEvent* event_C_ana;
static StPicoAlexEvent* event_D_ana;

static StPicoAlexEvent* event_SE_ME_ana[2];

static StPicoAlexEvent* PicoAlexEventA;
static StPicoAlexEvent* PicoAlexEventB;
static StPicoAlexEvent* PicoAlexEvent_use;

// output tree based on structure of StParticleEvent class
static const char PARTICLE_EVENT_TREE[]   = "ParticleEvents";
static const char PARTICLE_EVENT_BRANCH[] = "Events";
static TTree* my_tree;
static StParticleEvent* my_event;
static Float_t EP_Q_vectors_rc_glob[6]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
static Float_t Qtracks_used_Arr_glob[6];
static TVector2* Arr_iQadd_FullTPC[N_max_tracks];
static TVector2* Arr_iQadd_EtaPos[N_max_tracks]; 
static TVector2* Arr_iQadd_EtaNeg[N_max_tracks];


static TH2F *h_phi_vs_track_number[N_refmult_bins][2]; // [centrality bin][+ charge, - charge]
static TH2F *h_pT_vs_track_number[N_refmult_bins][2];  // [centrality bin][+ charge, - charge]
static TProfile *h_mean_pt_vs_phi[N_refmult_bins][2]; // [centrality bin][+ charge, - charge]
static TProfile *p_mean_pt_vs_charge_asymm[2];

static TH2F* hDecVertex_YX;
static TH2F* hDecVertex_YZ;
static TH1F* hInvMassAB[3];
static TH1F* hDcaDiff[3];
static TH1F* hnSigElDist[4];

static TTree* my_tree_single;
static MyEventEmbData* my_event_single;

static TTree* BBC_tree;
static Int_t    BBC_refMult;
static Int_t    BBC_runId;
static UShort_t BBC_ADC[48];
static Double_t BBC_x, BBC_y, BBC_z;

static TH2D* h2D_Omega2250_m_vs_p_cent[9];
static TH1D* h1D_Omega2250_m_cent[9];
static TH1D* h1D_Lambda;
static TH1D* h1D_Lambda_B;
static TH1D* h1D_Xi;


const static Int_t N_2D_m2_nSigmaP_pT_bins = 15;
const static Float_t low_pT_m2_nSigmaP     = 0.4;
const static Float_t delta_pT_m2_nSigmaP   = 0.2;
static TH2F* h2D_m2_nSigmaP[N_refmult_bins][N_2D_m2_nSigmaP_pT_bins][2]; // [][][+,-]


// Jet analysis
static const Int_t N_jet_histos = 10;
static TH1F* h_jet_pt[N_jet_histos];
static TH1F* h_jet_pt_sub[N_jet_histos];
static TH1F* h_jet_area[N_jet_histos];
static TH1F* h_jet_per_event[N_jet_histos];
static TRandom3 r3;

static const Int_t npt_bins_pt_spectra  = 20;
static const Int_t n_hadrons_pt_spectra = 2; // positive charge, negative charge
static Double_t pt_bin_ranges_spectra[npt_bins_pt_spectra] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0};

static TH2F* h_NewY_vs_NewX_mult_pt[9][npt_bins_pt_spectra][n_hadrons_pt_spectra];
static TH1F* h_m2_mult_pt[9][npt_bins_pt_spectra][n_hadrons_pt_spectra][2];
static TH2F* h_ToF_local_YZ_vs_m2[2];
static StCombPID* combPID;


static const Int_t N_InvMass_lambdaCplus = 12;
static TH2F* h_2D_InvMass_LambdaCplus[N_InvMass_lambdaCplus];

static const Int_t N_InvMass_lambdaCplus_V0        = 100;
static const Int_t N_InvMass_lambdaCplus_V0_Lambda = 8;
static TH1F* h_massA_LambdaCplus;
static TH1F* h_massB_LambdaCplus;
static TH1F* h_massC_LambdaCplus;
static TH1F* h_massA_corr_LambdaCplus;
static TH1F* h_massB_corr_LambdaCplus;
static TH1F* h_InvMass_AB_LambdaCplus[N_InvMass_lambdaCplus_V0_Lambda];
static TH1F* h_InvMass_ABC_LambdaCplus[N_InvMass_lambdaCplus_V0];


static TProfile* p_pt_Ach[2];
static TH1F*     h_Ach;

#include "Analysis_Snurf_Func.h"

ClassImp(Analysis_Snurf)
    //_____________________________________________________________________________
    Analysis_Snurf::Analysis_Snurf() {

    }
//_____________________________________________________________________________
Analysis_Snurf::~Analysis_Snurf() {

}

Int_t Analysis_Snurf::Init()
{
    cout << "***************** Initializing objects *******************" << endl;

    combPID         = new StCombPID();
    V0_tof_corr     = new StV0TofCorrection();

    r3.SetSeed(0);
    gRandom->SetSeed(0);

    //refmultCorrUtil = new StRefMultCorr();
    refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();

    cout << "MAGFIELDFACTOR = " << MAGFIELDFACTOR << endl;

    TString HistName,HistName2,HistName3;

    TString inFile_name =  inDir_asciifile;
    inFile_name         += inFile;
    inFile_name         += ".list";

    cout << "Read data from input file: " << inFile_name.Data() << endl;
    picoMaker = new StPicoDstMaker(0,inFile_name.Data(),"picoDst");
    picoMaker ->Init();
    input     = picoMaker->chain();

    n_total_entries = input->GetEntries();

    cout << "" << endl;
    cout << "************************************************************************************************************" << endl;
    cout << "Number of entries in tree = " << n_total_entries << endl;
    cout << "************************************************************************************************************" << endl;
    cout << "" << endl;

    if(nEvents > n_total_entries) nEvents = n_total_entries;

    // Initialize the counters
    memset(&Sample_Counter[0][0],0,sizeof(Int_t)*N_refmult_bins*N_z_axis_bins*N_EP_bins);
    memset(&Sample_Event_Num[0][0][0],0,sizeof(Int_t)*N_refmult_bins*N_z_axis_bins*N_max_sample*N_EP_bins);

    for(Int_t i = 0; i < N_refmult_bins; i++)
    {
        for(Int_t j = 0; j < N_z_axis_bins; j++)
        {
            for(Int_t k = 0; k < N_EP_bins; k++)
            {
                N_max_sample_array[i][j][k] = N_max_sample;
                //if(j < 2) N_max_sample_array[i][j]   = 5; // number has to be the same for all bins, otherwise the ratio is changed
                //if(j == 13) N_max_sample_array[i][j] = 5;
            }
        }
    }



    phi_corr_z_start   = -vertex_z_array[eBeamTimeNum];
    phi_corr_z_stop    = vertex_z_array[eBeamTimeNum];
    phi_corr_eta_start = -1.0; // -1.3
    phi_corr_eta_stop  = 1.0;  // 1.3
    //if(eBeamTimeNum == 0 || eBeamTimeNum == 1) // old values
    //{
    //    phi_corr_eta_start = -1.3; // -1.3
    //    phi_corr_eta_stop  = 1.3;  // 1.3
    //}
    phi_corr_pt_start  = 0.0;  // 0.0
    phi_corr_pt_stop   = 1.6;  // 1.6
    phi_corr_delta_z   = (phi_corr_z_stop-phi_corr_z_start)/((Float_t)nPhi_corr_z_bins);
    phi_corr_delta_eta = (phi_corr_eta_stop-phi_corr_eta_start)/((Float_t)nPhi_corr_eta_bins);
    phi_corr_delta_pt  = (phi_corr_pt_stop-phi_corr_pt_start)/((Float_t)nPhi_corr_pt_bins);


    //****************** Get the correction histograms for the event plane phi distribution ******************************
    cout << "Open the phi correction histograms for event plane analysis" << endl;

    Double_t MaxEntry    = 0.0;
    Int_t Good_phi_file  = 0;
    Int_t Total_phi_file = 0;
    Int_t Start_day_phi_corr = 0;
    Int_t Stop_day_phi_corr  = 0;

    // new
    TFile *filephi_day_z_eta;
    if(eBeamTimeNum == 0) // 7.7 GeV
    {
        Start_day_phi_corr = 110;
        Stop_day_phi_corr  = 150;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu7.root");  // open the file
    }
    if(eBeamTimeNum == 1) // 11.5 GeV
    {
        Start_day_phi_corr = 146;
        Stop_day_phi_corr  = 161;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu11.root");  // open the file
    }
    if(eBeamTimeNum == 2) // 39 GeV
    {
        Start_day_phi_corr = 97;
        Stop_day_phi_corr  = 115;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu39.root");  // open the file
    }
    if(eBeamTimeNum == 3) // 62.4 GeV 
    {
        Start_day_phi_corr = 77;
        Stop_day_phi_corr  = 100;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu62.root");  // open the file
    }
    if(eBeamTimeNum == 4) // 19.6 GeV  
    {
        Start_day_phi_corr = 110;
        Stop_day_phi_corr  = 123;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu19.root");  // open the file
    }
    if(eBeamTimeNum == 5) // 27 GeV 
    {
        Start_day_phi_corr = 170; 
        Stop_day_phi_corr  = 181;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu27.root");  // open the file
    }
    if(eBeamTimeNum == 6) // 200 GeV
    {
        Start_day_phi_corr = 170; 
        Stop_day_phi_corr  = 181;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu27.root");  // open the file
    }
    if(eBeamTimeNum == 7) // 14.5 GeV
    {
        Start_day_phi_corr = 170; 
        Stop_day_phi_corr  = 181;
        filephi_day_z_eta = TFile::Open("../Analysis/Corrections/Phi_weights_AuAu27.root");  // open the file
    }
    if(eBeamTimeNum == 8) // 200 GeV run14
    {
        Start_day_phi_corr = 170; 
        Stop_day_phi_corr  = 181;
        filephi_day_z_eta = TFile::Open("/global/homes/a/aschmah/STAR/Analysis/Corrections/Phi_weights_AuAu27.root");  // open the file
    }

    for(Int_t i = Start_day_phi_corr; i < Stop_day_phi_corr; i++)
    {
        hPhi_days_use[i] = 0;
        hPhi_days_in[i]  = 0;
        for(Int_t j = 0; j < nPhi_corr_z_bins; j++)
        {
            for(Int_t k = 0; k < nPhi_corr_eta_bins; k++)
            {
                for(Int_t l = 0; l < 2; l++)
                {
                    for(Int_t m = 0; m < nPhi_corr_pt_bins; m++)
                    {
                        HistName = "hPhi_corr_";
                        HistName += i;
                        HistName += "_";
                        HistName += j;
                        HistName += "_";
                        HistName += k;
                        HistName += "_";
                        HistName += l;
                        HistName += "_";
                        HistName += m;
                        if(
                           filephi_day_z_eta->FindObjectAny(HistName.Data()) != 0
                          )
                        {
                            hPhi_corr_in[i][j][k][l][m] = (TH1F*)filephi_day_z_eta->FindObjectAny(HistName.Data());
                            HistName = "hPhi_corr_in_";
                            HistName += i;
                            HistName += "_";
                            HistName += j;
                            HistName += "_";
                            HistName += k;
                            HistName += "_";
                            HistName += l;
                            HistName += "_";
                            HistName += m;
                            hPhi_corr_in[i][j][k][l][m]->SetName(HistName.Data());
                            MaxEntry = hPhi_corr_in[i][j][k][l][m]->GetBinContent(hPhi_corr_in[i][j][k][l][m]->GetMaximumBin());

                            // Calculate the mean of the 10 max entries to reduce fluctuations
                            Double_t Max_Mean_Val_array[10];
                            for(Int_t q = 0; q < 10; q++)
                            {
                                Max_Mean_Val_array[q] = 0.0;
                            }
                            for(Int_t p = 0; p < hPhi_corr_in[i][j][k][l][m]->GetNbinsX(); p++)
                            {
                                Double_t Phi_corr_in_val = hPhi_corr_in[i][j][k][l][m]->GetBinContent(p);
                                if(Phi_corr_in_val > Max_Mean_Val_array[9])
                                {
                                    Max_Mean_Val_array[9] = Phi_corr_in_val;
                                    for(Int_t q = 0; q < 9; q++)
                                    {
                                        if(Phi_corr_in_val > Max_Mean_Val_array[8-q])
                                        {
                                            Max_Mean_Val_array[8-q+1] = Max_Mean_Val_array[8-q];
                                            Max_Mean_Val_array[8-q]   = Phi_corr_in_val;
                                        }
                                    }
                                }
                            }
                            Double_t MaxEntry_mean = 0.0;
                            Int_t    MaxEntry_mean_counter = 0;
                            for(Int_t q = 0; q < 10; q++)
                            {
                                //cout << "Max_Mean_Val_array[" << q << "] = " << Max_Mean_Val_array[q] << endl;
                                if(Max_Mean_Val_array[q] > 0.0)
                                {
                                    MaxEntry_mean = MaxEntry_mean + Max_Mean_Val_array[q];
                                    MaxEntry_mean_counter++;
                                }
                            }
                            if(MaxEntry_mean_counter > 0)
                            {
                                MaxEntry_mean = MaxEntry_mean/((Double_t)MaxEntry_mean_counter);
                            }
                            //cout << "MaxEntry = " << MaxEntry << ", MaxEntry_mean = " << MaxEntry_mean << endl;

                            hPhi_corr_in_max_entry[i][j][k][l][m] = MaxEntry;

                            if(MaxEntry_mean > 0.0)
                            {
                                hPhi_corr_in[i][j][k][l][m]->Scale(1.0/MaxEntry_mean);
                                Total_phi_file++;
                            }
                            if(MaxEntry >= 100)
                            {
                                Good_phi_file++;
                            }

                            hPhi_days_in[i] = 1;
                            //if(MaxEntry < 100)
                            //cout << "day_bin = " << i << ", Phi correction file opened: " << HistName.Data() << ", MaxEntry = " << MaxEntry << ", hPhi_days_in = " << hPhi_days_in[i] << endl;
                        }
                    }
                }
            }
        }
    }

    for(Int_t i = 0; i < 365; i++)
    {
        for(Int_t j = 0; j < nPhi_corr_z_bins; j++)
        {
            for(Int_t k = 0; k < nPhi_corr_eta_bins; k++)
            {
                for(Int_t l = 0; l < 2; l++)
                {
                    for(Int_t m = 0; m < nPhi_corr_pt_bins; m++)
                    {
                        HistName = "hPhi_corr_";
                        HistName += i;
                        HistName += "_";
                        HistName += j;
                        HistName += "_";
                        HistName += k;
                        HistName += "_";
                        HistName += l;
                        HistName += "_";
                        HistName += m;
                        hPhi_corr[i][j][k][l][m] = new TH1F(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
                    }
                }
            }
        }
    }



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



    for(Int_t xa = 0; xa < 10; xa++)
    {
        HistName = "h_eta_phi_";
        HistName += xa;
        h_eta_phi[xa] = new TH2F(HistName.Data(),HistName.Data(),200,-TMath::Pi(),TMath::Pi(),200,-1.2,1.2);
    }
    h_nhitsFit = new TH2F("h_nhitsFit","h_nhitsFit",200,0,5,50,0,50);
    h_nhitsFit_TOF = new TH2F("h_nhitsFit_TOF","h_nhitsFit_TOF",200,0,5,50,0,50);


    cout << "" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << "Total number of phi correction files = " << Total_phi_file << endl;
    cout << "Files with max entry > 100 = " << Good_phi_file << endl;
    cout << "-------------------------------------------------------------------------" << endl;


    //filephi_day_z_eta->Close();

    //filephi->Close();
    //********************************************************************************************************************



    //****************** Get the cut files *******************************************************************************
    cout << "Open the cut files" << endl;

    TFile *file_kaon_cut;
    if(eBeamTimeNum == 0) // 7.7 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 1) // 11.5 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 2) // 39 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 3) // 62.4 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 4) // 19.6 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 5) // 27 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 6) // 27 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 7) // 14.5 GeV
    {
        file_kaon_cut = TFile::Open("../Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    if(eBeamTimeNum == 8) // 200 GeV run14
    {
        file_kaon_cut = TFile::Open("/global/homes/a/aschmah/STAR/Analysis/cuts/Kaon_cut_39_GeV_5.root");  // open the file
    }
    for(Int_t i = 0; i < np_bins; i++)
    {
        HistName = "cut_Kaon_";
        HistName += i;
        cut_Kaon[i] = (TCutG*)file_kaon_cut->FindObjectAny(HistName.Data()); // Kaon cut with a S/B > 5.0, starting from 1.2 GeV/c to 3.5 GeV/c
    }
    cout << "Cut files opened" << endl;
    //********************************************************************************************************************


    cout << "Defining output" << endl;
    TString fulloutFileName =  outputDir;
    fulloutFileName         += inFile.Data();
    fulloutFileName         += outFileName.Data();

    Outputfile = new TFile(fulloutFileName,"RECREATE");
    cout << "Output file: " << fulloutFileName.Data() << endl;



    cout << "*******************************************************************************************************" << endl;
    cout << "Defining ntuples" << endl;
    if(eAnalysisNum == 1 || eAnalysisNum == 12 || eAnalysisNum == 25)
    {
        Outputfile->cd();
        Tree_PhiMeson_v2 = NULL;
        alexPhiMeson_event_ptr = &alexPhiMeson_event;
        Tree_PhiMeson_v2 = new TTree( ALEXPHIMESON_EVENT_TREE, ALEXPHIMESON_EVENT_TREE );
        Tree_PhiMeson_v2->Branch( ALEXPHIMESON_EVENT_BRANCH, "StAlexPhiMesonEvent", &alexPhiMeson_event_ptr );
        Tree_PhiMeson_v2->SetAutoSave( 5000000 );

        //PhiKPKM_NT      = new TNtuple("PhiKPKM_NT","PhiKPKM_NT Ntuple","InvMassAB:Mass2A:Mass2B:nSigmaKaonA:nSigmaKaonB:MomentumA:MomentumB:BetaA:BetaB:refMult:dcaA:dcaB:pt:rap:DeltaDipAngle:phi_event_plane:phiAB:TriggerId:EventVertexZ");
        //PhiKPKM_NT      ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 400)
    {
        Outputfile->cd();
        Tree_D0_v2 = NULL;
        D0_event_ptr = &D0_event;
        Tree_D0_v2 = new TTree( D0_EVENT_TREE, D0_EVENT_TREE );
        Tree_D0_v2->Branch( D0_EVENT_BRANCH, "StD0Event", &D0_event_ptr );
        Tree_D0_v2->SetAutoSave( 5000000 );

        h_dca_diff           = new TH1D("h_dca_diff","h_dca_diff",500,-0.1,0.1);
        h_decay_length_diff  = new TH1D("h_decay_length_diff","h_decay_length_diff",500,-0.1,0.1);
    }
    if(eAnalysisNum == 2 || eAnalysisNum == 14 || eAnalysisNum == 3)
    {
        //Outputfile->cd();
        //Tree_V0_v2 = NULL;
        //alexV0_event_ptr = &alexV0_event;
        //Tree_V0_v2 = new TTree( ALEXV0_EVENT_TREE, ALEXV0_EVENT_TREE );
        //Tree_V0_v2->Branch( ALEXV0_EVENT_BRANCH, "StAlexV0Event", &alexV0_event_ptr );
        //Tree_V0_v2->Write();
        //Tree_V0_v2->SetAutoSave( 5000000 );

        Lambda_X_NT     = new TNtuple("Lambda_X_NT","Lambda_X_NT Ntuple","InvM:InvM_K0S:m2A:m2B:nSP:nSPi:pA:pB:dcaA:dcaB:iQxA:iQyA:iQxB:iQyB:etaA:etaB:mult:dcaAB:VerdistX:VerdistY:pt:rap:phiAB:thetaAB:phi_EP:phi_EP_eta_gap:delta_phi_ME_AB:RunId:n_prim:n_non_prim:n_tof_prim:scalarProduct:x:y:z:EP_Qx_eta_pos_ptw:EP_Qy_eta_pos_ptw:EP_Qx_eta_neg_ptw:EP_Qy_eta_neg_ptw:EP_Qx_ptw:EP_Qy_ptw:Qtracks_eta_pos:Qtracks_eta_neg:Qtracks_all:ZDCx:vzVpd");
        Lambda_X_NT     ->SetAutoSave( 5000000 );

        Vertex_lin_NT     = new TNtuple("Vertex_lin_NT","Vertex_lin_NT Ntuple","x1:y1:z1:x2:y2:z2:x3:y3:z3:flag:pA:pB:dcaAB");
        //Vertex_lin_NT     ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 24)
    {
        ThetaPlus_NT     = new TNtuple("ThetaPlus_NT","ThetaPlus_NT Ntuple","InvMAB:InvMABC:m2A:m2B:m2C:nsPP:nsPM:nsP:pA:pB:pC:dcaA:dcaB:dcaC:mult:dcaAB:VerdistX:VerdistY:pt:rap:id:scalprod:X:Y:Z");
        ThetaPlus_NT     ->SetAutoSave( 5000000 );
    }

    if(eAnalysisNum == 35)
    {
        for(Int_t i = 0; i < N_InvMass_lambdaCplus; i++)
        {
            HistName = "h_2D_InvMass_LambdaCplus_";
            HistName += i;
            h_2D_InvMass_LambdaCplus[i] = new TH2F(HistName.Data(),HistName.Data(),100,0,5,1000,1.45,2.45);
        }
    }
    if(eAnalysisNum == 351)
    {
        for(Int_t i = 0; i < N_InvMass_lambdaCplus_V0_Lambda; i++)
        {
            HistName = "h_InvMass_AB_LambdaCplus_";
            HistName += i;
            h_InvMass_AB_LambdaCplus[i] = new TH1F(HistName.Data(),HistName.Data(),2000,0.23,1.23);  // 0.23..1.23
        }
        for(Int_t i = 0; i < N_InvMass_lambdaCplus_V0; i++)
        {
            HistName = "h_InvMass_ABC_LambdaCplus_";
            HistName += i;
            h_InvMass_ABC_LambdaCplus[i] = new TH1F(HistName.Data(),HistName.Data(),10000,1.23,6.23);
        }

        h_massA_LambdaCplus      = new TH1F("h_massA_LambdaCplus","h_massA_LambdaCplus",500,-0.3,2.3);
        h_massB_LambdaCplus      = new TH1F("h_massB_LambdaCplus","h_massB_LambdaCplus",500,-0.3,2.3);
        h_massC_LambdaCplus      = new TH1F("h_massC_LambdaCplus","h_massC_LambdaCplus",500,-0.3,2.3);
        h_massA_corr_LambdaCplus = new TH1F("h_massA_corr_LambdaCplus","h_massA_corr_LambdaCplus",500,-0.3,2.3);
        h_massB_corr_LambdaCplus = new TH1F("h_massB_corr_LambdaCplus","h_massB_corr_LambdaCplus",500,-0.3,2.3);

    }

    if(eAnalysisNum == 360)
    {
        for(Int_t i_charge = 0; i_charge < 2; i_charge++)
        {
            HistName = "p_pt_Ach_";
            HistName += i_charge;
            p_pt_Ach[i_charge] = new TProfile(HistName.Data(),HistName.Data(),15,-0.25,0.25);
        }

        h_Ach           = new TH1F("h_Ach","h_Ach",60,-0.3,0.3);
    }

    //if(eAnalysisNum == 4 || eAnalysisNum == 16  || eAnalysisNum == 20)
    //{
    //    Xi_NT        = new TNtuple("Xi_NT","Xi_NT Ntuple","InvMassAB:InvMassABC:MomentumA:MomentumB:MomentumB2:refMultA:dcaA:dcaB:dcaAB:dcaB2:dcaBB2:VerdistX:VerdistY:VerdistX2:VerdistY2:pt2:rap2:phiABC:phi_event_plane:TriggerId:EventVertexZ:Mass2B2:nSigmaPionB2:TPCdEdxB2:BetaB2:PolarityB2:Mass2A:Mass2B:Mass2Ao:Mass2Bo:Mass2B2o:BetaA:BetaB:phi_event_plane_eta_gap:dca_xi:fDCA_Helix_out:RunIdA:n_primaries:n_non_primaries:n_tofmatch_prim:scalarProduct:scalarProduct2");
    //    Xi_NT        ->SetAutoSave( 5000000 );
    //}
    //if(eAnalysisNum == 5 || eAnalysisNum == 20)
    //{
    //    Omega_NT        = new TNtuple("Omega_NT","Omega_NT Ntuple","InvMassAB:InvMassABC:MomentumA:MomentumB:MomentumB2:refMultA:dcaA:dcaB:dcaAB:dcaB2:dcaBB2:VerdistX:VerdistY:VerdistX2:VerdistY2:pt2:rap2:phiABC:phi_event_plane:TriggerId:EventVertexZ:Mass2B2:nSigmaKaonB2:TPCdEdxB2:BetaB2:PolarityB2:Mass2A:Mass2B:Mass2Ao:Mass2Bo:Mass2B2o:BetaA:BetaB:phi_event_plane_eta_gap:dca_xi:fDCA_Helix_out:RunIdA:n_primaries:n_non_primaries:n_tofmatch_prim:scalarProduct:scalarProduct2:InvMassABC_Xi");
    //    Omega_NT        ->SetAutoSave( 5000000 );
    //}
    if(eAnalysisNum == 20)
    {
        Outputfile->cd();
        Tree_OmegaPV0_v2   = NULL;
        Tree_OmegaMV0_v2   = NULL;
        Tree_XiPV0_v2      = NULL;
        Tree_XiMV0_v2      = NULL;
        alexV0_event_ptr_A = &alexV0_event_A;
        alexV0_event_ptr_B = &alexV0_event_B;
        Tree_OmegaPV0_v2   = new TTree("OmegaP_v2_tree" , ALEXV0_EVENT_TREE );
        Tree_OmegaMV0_v2   = new TTree("OmegaM_v2_tree" , ALEXV0_EVENT_TREE );
        Tree_XiPV0_v2      = new TTree("XiP_v2_tree"    , ALEXV0_EVENT_TREE );
        Tree_XiMV0_v2      = new TTree("XiM_v2_tree"    , ALEXV0_EVENT_TREE );
        Tree_OmegaPV0_v2   ->Branch("OmegaP_v2_branch"  , "StAlexV0Event", &alexV0_event_ptr_A );
        Tree_OmegaMV0_v2   ->Branch("OmegaM_v2_branch"  , "StAlexV0Event", &alexV0_event_ptr_A );
        Tree_XiPV0_v2      ->Branch("XiP_v2_branch"     , "StAlexV0Event", &alexV0_event_ptr_B );
        Tree_XiMV0_v2      ->Branch("XiM_v2_branch"     , "StAlexV0Event", &alexV0_event_ptr_B );
    }
    if(eAnalysisNum == 102)
    {
        Outputfile->cd();
        Tree_Sigma0_v2   = NULL;
        alexV0_event_ptr_A = &alexV0_event_A;
        Tree_Sigma0_v2   = new TTree("Tree_Sigma0_v2_tree" , ALEXV0_EVENT_TREE );
        Tree_Sigma0_v2   ->Branch("Tree_Sigma0_v2_branch"  , "StAlexV0Event", &alexV0_event_ptr_A );
    }
    if(eAnalysisNum == 8)
    {
        Outputfile->cd();
        Tree_hadron_v2 = NULL;
        alex_event_ptr = &alex_event;
        Tree_hadron_v2 = new TTree( ALEX_EVENT_TREE, ALEX_EVENT_TREE );
        Tree_hadron_v2->Branch( ALEX_EVENT_BRANCH, "StAlexEvent", &alex_event_ptr );
        Tree_hadron_v2->Write();
        Tree_hadron_v2->SetAutoSave( 5000000 );

        //Hadron_v2_NT    = new TNtuple("Hadron_v2_NT","Hadron_v2_NT Ntuple","phi:p_t:eta:dca:pq:m2:nsK:nsPi:nsP:nsEl:dEdx:iQx:iQy:x:y:z:id:mult:phi_EP:phi_EP_raw:phi_EP_eta_gap:n_prim:n_non_prim:n_tofmatch_prim:EP_Qx_eta_pos_ptw:EP_Qy_eta_pos_ptw:EP_Qx_eta_neg_ptw:EP_Qy_eta_neg_ptw:EP_Qx_ptw:EP_Qy_ptw:Qtracks_eta_pos:Qtracks_eta_neg:Qtracks_full");
        //Hadron_v2_NT    ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 88)
    {
        const Int_t n_bins_X_pt[npt_bins_pt_spectra] = {500,500,500,400,400,300,300,300,300,200,200,200,200,100,100,100,100,100,100,100};
        const Int_t n_bins_Y_pt[npt_bins_pt_spectra] = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100};

        Double_t start_stop_X_Y_pt[2][2][npt_bins_pt_spectra]; // [start,stop][X,Y][pt_bin]
        for(Int_t j = 0; j < 4; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.4;
            start_stop_X_Y_pt[1][0][j] = 1.2;
            start_stop_X_Y_pt[0][1][j] = -0.2;
            start_stop_X_Y_pt[1][1][j] = 0.2;
        }
        for(Int_t j = 4; j < 8; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.4;
            start_stop_X_Y_pt[1][0][j] = 1.3;
            start_stop_X_Y_pt[0][1][j] = -0.3;
            start_stop_X_Y_pt[1][1][j] = 0.2;
        }
        for(Int_t j = 8; j < 10; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.4;
            start_stop_X_Y_pt[1][0][j] = 1.3;
            start_stop_X_Y_pt[0][1][j] = -0.6;
            start_stop_X_Y_pt[1][1][j] = 0.3;
        }
        for(Int_t j = 10; j < 12; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.4;
            start_stop_X_Y_pt[1][0][j] = 1.3;
            start_stop_X_Y_pt[0][1][j] = -1.0;
            start_stop_X_Y_pt[1][1][j] = 0.6;
        }
        for(Int_t j = 12; j < 14; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.5;
            start_stop_X_Y_pt[1][0][j] = 1.4;
            start_stop_X_Y_pt[0][1][j] = -1.2;
            start_stop_X_Y_pt[1][1][j] = 0.7;
        }
        for(Int_t j = 14; j < 16; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.6;
            start_stop_X_Y_pt[1][0][j] = 1.4;
            start_stop_X_Y_pt[0][1][j] = -1.4;
            start_stop_X_Y_pt[1][1][j] = 0.8;
        }
        for(Int_t j = 16; j < 18; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.7;
            start_stop_X_Y_pt[1][0][j] = 1.6;
            start_stop_X_Y_pt[0][1][j] = -1.5;
            start_stop_X_Y_pt[1][1][j] = 0.9;
        }
        for(Int_t j = 18; j < npt_bins; j++)
        {
            start_stop_X_Y_pt[0][0][j] = -0.8;
            start_stop_X_Y_pt[1][0][j] = 1.8;
            start_stop_X_Y_pt[0][1][j] = -1.6;
            start_stop_X_Y_pt[1][1][j] = 1.0;
        }


        for(Int_t i = 0; i < 2; i++)
        {
            HistName = "h_ToF_local_YZ_vs_m2_";
            HistName += i;
            h_ToF_local_YZ_vs_m2[i] = new TH2F(HistName.Data(),HistName.Data(),500,-0.4,1.6,200,-7,7);
        }
        for(Int_t i = 0; i < 9; i++) // 8
        {
            for(Int_t j = 0; j < npt_bins_pt_spectra; j++) // 20
            {
                for(Int_t k = 0; k < n_hadrons_pt_spectra; k++) // 2
                {
                    HistName = "h_NewY_vs_NewX_mult_pt_m";
                    HistName += i;
                    HistName += "_pt_";
                    HistName += j;
                    HistName += "_pid_";
                    HistName += k;
                    h_NewY_vs_NewX_mult_pt[i][j][k] = new TH2F(HistName.Data(),HistName.Data(),n_bins_X_pt[j],start_stop_X_Y_pt[0][0][j],start_stop_X_Y_pt[1][0][j],n_bins_Y_pt[j],start_stop_X_Y_pt[0][1][j],start_stop_X_Y_pt[1][1][j]);
                    h_NewY_vs_NewX_mult_pt[i][j][k] ->Sumw2();

                    for(Int_t l = 0; l < 2; l++)
                    {
                        HistName = "h_m2_mult_pt_m";
                        HistName += i;
                        HistName += "_pt_";
                        HistName += j;
                        HistName += "_pid_";
                        HistName += k;
                        HistName += "_cut_";
                        HistName += l;
                        h_m2_mult_pt[i][j][k][l] = new TH1F(HistName.Data(),HistName.Data(),800,-0.4,1.6);
                        h_m2_mult_pt[i][j][k][l] ->Sumw2();
                    }
                }
            }
        }
    }
    if(eAnalysisNum == 130)
    {
        Outputfile->cd();
        Float_t low_nSigma_range_pT[N_2D_m2_nSigmaP_pT_bins] = {-24,-16,-11,-9,-7,-5,-4,-4,-4,-4,-4,-4,-4,-4,-4};
        Float_t N_bins_nSigma[N_2D_m2_nSigmaP_pT_bins]       = {300,300,300,300,300,300,300,300,300,300,300,300,300,300,300};
        Float_t N_bins_m2[N_2D_m2_nSigmaP_pT_bins]           = {600,600,500,500,400,300,300,300,300,300,300,300,300,300,300};
        for(Int_t nrefmult = 0; nrefmult < N_refmult_bins; nrefmult++)
        {
            for(Int_t npt = 0; npt < N_2D_m2_nSigmaP_pT_bins; npt++)
            {
                for(Int_t ncharge = 0; ncharge < 2; ncharge++)
                {
                    HistName = "h2D_m2_nSigmaP_refmult";
                    HistName += nrefmult;
                    HistName += "_pt";
                    HistName += npt;
                    HistName += "_charge";
                    HistName += ncharge;
                    h2D_m2_nSigmaP[nrefmult][npt][ncharge] = new TH2F(HistName.Data(),HistName.Data(),N_bins_nSigma[npt],low_nSigma_range_pT[npt],5,N_bins_m2[npt],-1,2);
                }
            }
        }
    }
    if(eAnalysisNum == 112)
    {
        BBC_tree = new TTree("BBC_tree","BBC_tree");
        BBC_tree->Branch("refMult",&BBC_refMult,"BBC_refMult/I");
        BBC_tree->Branch("x",&BBC_x,"BBC_x/D");
        BBC_tree->Branch("y",&BBC_y,"BBC_y/D");
        BBC_tree->Branch("z",&BBC_z,"BBC_z/D");
        BBC_tree->Branch("runId",&BBC_runId,"BBC_runId/I");
        BBC_tree->Branch("TPC_Qx1_eta_pos",&EP_Qx1_eta_pos_ptw,"TPC_Qx1_eta_pos/F");
        BBC_tree->Branch("TPC_Qy1_eta_pos",&EP_Qy1_eta_pos_ptw,"TPC_Qy1_eta_pos/F");
        BBC_tree->Branch("TPC_Qx1_eta_neg",&EP_Qx1_eta_neg_ptw,"TPC_Qx1_eta_neg/F");
        BBC_tree->Branch("TPC_Qy1_eta_neg",&EP_Qy1_eta_neg_ptw,"TPC_Qy1_eta_neg/F");
        BBC_tree->Branch("TPC_Qx2_eta_pos",&EP_Qx_eta_pos_ptw,"TPC_Qx2_eta_pos/F");
        BBC_tree->Branch("TPC_Qy2_eta_pos",&EP_Qy_eta_pos_ptw,"TPC_Qy2_eta_pos/F");
        BBC_tree->Branch("TPC_Qx2_eta_neg",&EP_Qx_eta_neg_ptw,"TPC_Qx2_eta_neg/F");
        BBC_tree->Branch("TPC_Qy2_eta_neg",&EP_Qy_eta_neg_ptw,"TPC_Qy2_eta_neg/F");
        BBC_tree->Branch("TPC_N_tracks_eta_pos",&Qtracks_used_eta_pos,"TPC_N_tracks_eta_pos/I");
        BBC_tree->Branch("TPC_N_tracks_eta_neg",&Qtracks_used_eta_neg,"TPC_N_tracks_eta_neg/I");
        BBC_tree->Branch("ADC",BBC_ADC,"BBC_ADC[48]/s");
        //BBC_tree->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 11)
    {
        EventPlane_NT   = new TNtuple("EventPlane_NT","EventPlane_NT Ntuple","Qx:Qy:EventVertexX:EventVertexY:EventVertexZ:refMult:event_number:Qtracks_used:phi_mean:EP_Qx_subA:EP_Qy_subA:EP_Qx_subB:EP_Qy_subB:mean_px:mean_py:mean_pz:mean_pxw:mean_pyw:mean_pzw:Psi:Psi_phi_weight:Psi_no_weight:EP_Qx_eta_pos:EP_Qy_eta_pos:EP_Qx_eta_neg:EP_Qy_eta_neg:Psi_eta_pos:Psi_eta_neg:Qtracks_used_eta_pos:Qtracks_used_eta_neg:mean_py_pos:mean_py_neg:mean_px_pos:mean_px_neg:EP_Qx_eta_pos_ptw:EP_Qy_eta_pos_ptw:EP_Qx_eta_neg_ptw:EP_Qy_eta_neg_ptw:EP_Qx_ptw:EP_Qy_ptw:n_tofmatch_prim:n_non_primaries:RunId:ZDCx:BBCx:vzVpd:EP_Qx_subA_ptw:EP_Qy_subA_ptw:EP_Qx_subB_ptw:EP_Qy_subB_ptw");
        EventPlane_NT   ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 200 || eAnalysisNum == 201 || eAnalysisNum == 203)
    {
        for(Int_t m = 0; m < 9; m++)
        {
            HistName = "h2D_Omega2250_m_vs_p_cent";
            HistName += m;
            h2D_Omega2250_m_vs_p_cent[m] = new TH2D(HistName.Data(),HistName.Data(),100,0,6,2800,1.9,8.9);

            HistName = "h1D_Omega2250_m_cent";
            HistName += m;
            h1D_Omega2250_m_cent[m] = new TH1D(HistName.Data(),HistName.Data(),2800,1.9,8.9);
        }
        h1D_Lambda   = new TH1D("h1D_Lambda","h1D_Lambda",200,1.0,1.2);
        h1D_Lambda_B = new TH1D("h1D_Lambda_B","h1D_Lambda_B",200,1.0,1.2);
        h1D_Xi       = new TH1D("h1D_Xi","h1D_Xi",200,1.2,1.4);

        d4s_NT        = new TNtuple("d4s_NT","d4s_NT Ntuple","InvABCDE:InvAB:InvABC:InvDE:VDX_AB:VDX_ABC:VDX_DE:VDX_ABCDE:refMult:P:VDY_AB:VDY_ABC:VDY_DE:dcaA:dcaB:dcaC:dcaD:dcaE:dcaABC_DE:dcaDE:dcaAB:dcaC_AB:dcaABCDE:refMultB");
        d4s_NT        ->SetAutoSave( 5000000 );

        Tree_d4s = NULL;
        d4s_event_ptr = &d4s_event;
        Tree_d4s = new TTree(D4S_EVENT_TREE ,  D4S_EVENT_TREE);
        Tree_d4s->Branch(D4S_EVENT_BRANCH , "Std4sEvent", &d4s_event_ptr );
        Tree_d4s->SetAutoSave( 5000000 );

        for(Int_t i = 0; i < 4; i++)
        {
            HistName = "h_d4s_InvMassAB_";
            HistName += i;
            h_d4s_InvMassAB[i] = new TH1F(HistName.Data(),HistName.Data(),800,0.4,2.0);
        }
    }
    if(eAnalysisNum == 202)
    {
        Tree_Lambda = NULL;
        Lambda_event_ptr = &Lambda_event;
        Tree_Lambda = new TTree(LAMBDA_EVENT_TREE ,  LAMBDA_EVENT_TREE);
        Tree_Lambda->Branch(LAMBDA_EVENT_BRANCH , "StLambdaEvent", &Lambda_event_ptr );
        Tree_Lambda->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 302)
    {
        Tree_JetTrackEvent = NULL;
        JetTrackEvent_ptr = &JetTrackEvent;
        Tree_JetTrackEvent = new TTree(JETTRACK_EVENT_TREE ,  JETTRACK_EVENT_TREE);
        Tree_JetTrackEvent->Branch(JETTRACK_EVENT_BRANCH , "StJetTrackEvent", &JetTrackEvent_ptr );
        Tree_JetTrackEvent->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 300 || eAnalysisNum == 301)
    {
        for(Int_t i = 0; i < N_jet_histos; i++)
        {
            HistName = "h_jet_pt_";
            HistName += i;
            h_jet_pt[i]        = new TH1F(HistName.Data(),HistName.Data(),700,-20,50.0);

            HistName = "h_jet_pt_sub_";
            HistName += i;
            h_jet_pt_sub[i]    = new TH1F(HistName.Data(),HistName.Data(),700,-20,50.0);

            HistName = "h_jet_area_";
            HistName += i;
            h_jet_area[i]      = new TH1F(HistName.Data(),HistName.Data(),500,0.0,10.0);

            HistName = "h_jet_per_event_";
            HistName += i;
            h_jet_per_event[i] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,100.0);
        }
    }
    if(eAnalysisNum == 125)
    {
        EventAnalysis_NT   = new TNtuple("EventAnalysis_NT","EventAnalysis_NT Ntuple","X:Y:Z:refmult:Npos:Nneg:mpt_pos:mpt_neg:NposB:NnegB:mpt_posB:mpt_negB:Apm");
        EventAnalysis_NT   ->SetAutoSave( 5000000 );

        for(Int_t k = 0; k < 2; k++)
        {
            HistName = "p_mean_pt_vs_charge_asymm_";
            HistName += k;
            p_mean_pt_vs_charge_asymm[k] = new TProfile(HistName.Data(),HistName.Data(),500,-1,1);
        }
    }
    if(eAnalysisNum == 111)
    {
        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics (5)
        {
            HistName = "EventPlane_array_NT_";
            HistName += n_harm;
            HistName2 = "EventPlane_array_NT_";
            HistName2 += n_harm;
            HistName2 += " Ntuple";

            HistName3 = "X:Y:Z:mult:ev:evid:nTOFprim:nnonPrim:RunId:ZDCx:BBCx:vzVpd:Qx_full:Qy_full:Qx_full_pt:Qy_full_pt:Qx_full_subA:Qy_full_subA:Qx_full_subB:Qy_full_subB:Qx_full_pt_subA:Qy_full_pt_subA:Qx_full_pt_subB:Qy_full_pt_subB:N_tracks_full";
            HistName3 += ":QxE_BBC:QyE_BBC:QxW_BBC:QyW_BBC:N_BBCE:N_BBCW";
            for(Int_t n_eta_gap = 0; n_eta_gap < n_eta_gap_values; n_eta_gap++)  // loop over different eta gap values (6)
            {
                for(Int_t n_charge = 0; n_charge < 3; n_charge++) // loop over charges: all, +, -
                {
                    HistName3 += ":Qx_eta_pos_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":Qy_eta_pos_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":Qx_eta_neg_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":Qy_eta_neg_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;

                    HistName3 += ":Qx_eta_pt_pos_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":Qy_eta_pt_pos_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":Qx_eta_pt_neg_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":Qy_eta_pt_neg_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;

                    HistName3 += ":N_tracks_eta_pos_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                    HistName3 += ":N_tracks_eta_neg_";
                    HistName3 += n_eta_gap;
                    HistName3 += "_c";
                    HistName3 += n_charge;
                }
            }

            EventPlane_array_NT[n_harm]   = new TNtuple(HistName.Data(),HistName2.Data(),HistName3.Data());
            EventPlane_array_NT[n_harm]   ->SetAutoSave( 5000000 );
        }
    }
    if(eAnalysisNum == 13)
    {
        PhiKPKM_V0_NT      = new TNtuple("PhiKPKM_V0_NT","PhiKPKM_V0_NT Ntuple","InvMassAB:Mass2A:Mass2B:nSigmaKaonA:nSigmaKaonB:MomentumA:MomentumB:BetaA:BetaB:refMult:dcaA:dcaB:pt:rap:DeltaDipAngle:phi_event_plane:phiAB:TriggerId:EventVertexZ:VerdistX:InvMassAB_K0S:InvMassAB_Lambda:InvMassAB_KStar:InvMassAB_AntiKStar:InvMassABPrim:RunId:phi_event_plane_eta_gap:n_primaries:n_non_primaries:n_tofmatch_prim:scalarProduct:VerdistY:dcaAB");
        PhiKPKM_V0_NT      ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 17 || eAnalysisNum == 18)
    {
        Xi1530_NT        = new TNtuple("Xi1530_NT","Xi1530_NT Ntuple","InvMassAB:InvMassABC:InvMassABCD:refMultA:dcaA:dcaB:dcaAB:dcaB2:dcaBB2:VerdistX:VerdistY:VerdistX2:VerdistY2:TriggerIdA:EventVertexZA:Mass2B2Corr:TPCdEdxB2:BetaB2:PolarityB2:Mass2ACorr:Mass2BCorr:BetaA:BetaB:dca_omega:fDCA_Helix_out:RunIdA:n_primaries:n_non_primaries:n_tofmatch_prim");
        Xi1530_NT        ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 22)
    {
        RefMult_QA_NT        = new TNtuple("RefMult_QA_NT","RefMult_QA_NT Ntuple","mult:x:y:z:id:prim:primt:glob:mean_apt:mean_ppt:BG_rate:ZDCx:BBCx:vzVpd:Qx_eta_pos:Qy_eta_pos:Qx_eta_neg:Qy_eta_neg:multB:globB:multC:multD:multE:trig:Qx:Qy:Neta_pos:Neta_neg:N_EP:Qx_QA:Qy_QA:Qx_eta_pos_QA:Qy_eta_pos_QA:Neta_pos_QA:Qx_eta_neg_QA:Qy_eta_neg_QA:Neta_neg_QA:N_EP_QA");
        RefMult_QA_NT        ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 23)
    {
        Outputfile->cd();
        Tree_KaonV0_v2     = NULL;
        alexV0_event_ptr_A = &alexV0_event_A;
        Tree_KaonV0_v2   = new TTree("Kaon_v2_tree" , ALEXV0_EVENT_TREE );
        Tree_KaonV0_v2   ->Branch("Kaon_v2_branch"  , "StAlexV0Event", &alexV0_event_ptr_A );

        //Kaon_NT        = new TNtuple("Kaon_NT","Kaon_NT Ntuple","VerdistX:InvMassABC:MomentumABC:dcaABC:dcaAB:pt2:rap2:dcaA:dcaB:dcaC:MomentumA:MomentumB:MomentumC:nSigmaPionA:nSigmaPionB:nSigmaPionC:dca_Kaon");
        //Kaon_NT        ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 9)
    {
        cout << "Defining Vertex_NT" << endl;
        Vertex_NT       = new TNtuple("Vertex_NT","Vertex_NT Ntuple","x:y:z:ntracks:xo:yo:zo:comb_counter:goodHLT:BBC:VPD:N_primaries_old:N_primaries_new:fraction_comb_counter:eta_sum_old:eta_sum_new:x_mean2:y_mean2:z_mean2:N_primaries_new2:refMult:N_tofmatch_new:N_tofmatch_new2:vector_dist:RunId:N_tofmatch_old:TriggerId:xA:yA:zA:xB:yB:zB:EventId");
        Vertex_NT       ->SetAutoSave( 5000000 );

        cout << "Defining EDisplay_NT" << endl;
        EDisplay_NT       = new TNtuple("EDisplay_NT","EDisplay_NT Ntuple","x:y:z:ntracks:xo:yo:zo:x2:y2:z2:mult:tof:tof2:vecdist:id:npid:px:py:pz:origx:origy:origz:bfield:charge:dca:beta:tofx:tofy:tofz:event:e1:e2:e3:e:e0:btowId:m2:nSigmaPion:nSigmaKaon:nSigmaProton:nSigmaElectron");
        EDisplay_NT       ->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 26)
    {
        Outputfile->cd();
        my_tree = new TTree( PARTICLE_EVENT_TREE, PARTICLE_EVENT_TREE );
        my_event = new StParticleEvent();
        my_tree->Branch( PARTICLE_EVENT_BRANCH, "StParticleEvent", &my_event);
        my_tree->Write();
        my_tree->SetAutoSave(5000000);
    }

    cout << "All Ntuples defined" << endl;


    if(eAnalysisNum == 111)
    {
        heta_EP = new TH1F("heta_EP","heta_EP",200,-6,6);
        hFTPC_BBC_corr = new TH2F("hFTPC_BBC_corr","hFTPC_BBC_corr",20,0,20,200,0,800);
        h_FTPC_phi     = new TH1F("h_FTPC_phi","h_FTPC_phi",100,-TMath::Pi(),TMath::Pi());
        for(Int_t n_ew = 0; n_ew < 2; n_ew++)
        {
            HistName = "hrefMult_BBC_hits_East_West_";
            HistName += n_ew;
            hrefMult_BBC_hits_East_West[n_ew] = new TH2F(HistName.Data(),HistName.Data(),500,0,500,50,0,50);

            for(Int_t n_tiles = 0; n_tiles < 24; n_tiles++)
            {
                HistName = "hBBC_ADC_tiles_Eeast_West_";
                HistName += n_ew;
                HistName += "_";
                HistName += n_tiles;
                hBBC_ADC_tiles_Eeast_West[n_ew][n_tiles] = new TH1F(HistName.Data(),HistName.Data(),2000,0,2000);

                HistName = "pBBC_ADC_tiles_vs_RefMult_Eeast_West_";
                HistName += n_ew;
                HistName += "_";
                HistName += n_tiles;
                pBBC_ADC_tiles_vs_RefMult_Eeast_West[n_ew][n_tiles] = new TProfile(HistName.Data(),HistName.Data(),500,0,500);

                HistName = "hBBC_ADC_tiles_vs_RefMult_Eeast_West_";
                HistName += n_ew;
                HistName += "_";
                HistName += n_tiles;
                hBBC_ADC_tiles_vs_RefMult_Eeast_West[n_ew][n_tiles] = new TH2F(HistName.Data(),HistName.Data(),200,0,500,200,0,1000);
            }
        }
    }

    if ( eAnalysisNum == 124)
    {
        hDecVertex_YZ = new TH2F("hDecVertex_YZ","hDecVertex_YZ",500,-80,80,400,-60,60);
        hDecVertex_YX = new TH2F("hDecVertex_YX","hDecVertex_YX",400,-60,60,400,-60,60);
        for ( Int_t i = 0; i < 3; ++i )
        {
            hInvMassAB[i] = new TH1F(Form("hInvMassAB[%d]",i),Form("hInvMassAB[%d]",i),250,0,0.25);
            hDcaDiff[i] = new TH1F(Form("hDcaDiff[%d]",i),Form("hDcaDiff[%d]",i),60,0,30);
        }
        for ( Int_t i = 0; i < 4; ++i ) hnSigElDist[i] = new TH1F(Form("hnSigElDist[%d]",i),Form("hnSigElDist[%d]",i),460,-2.3,2.3);

        my_tree_single = new TTree("SingleTracksTree","Tree of Single Tracks w/ event structure");
        my_event_single = new MyEventEmbData();
        my_tree_single->Branch("Event","MyEventEmbData",&my_event_single,64000,2);
        my_tree_single->SetAutoSave(5000000);
    }

    if( eAnalysisNum == 31)
    {
        for(Int_t n_mult = 0; n_mult < N_refmult_bins; n_mult++)
        {
            for(Int_t n_w = 0; n_w < 2; n_w++)
            {
                HistName = "h_phi_vs_track_number_m";
                HistName += n_mult;
                HistName += "_c";
                HistName += n_w;
                h_phi_vs_track_number[n_mult][n_w] = new TH2F(HistName.Data(),HistName.Data(),1000,0,1000,100,0,TMath::Pi()*2.0);
                HistName = "h_pT_vs_track_number_m";
                HistName += n_mult;
                HistName += "_c";
                HistName += n_w;
                h_pT_vs_track_number[n_mult][n_w]  = new TH2F(HistName.Data(),HistName.Data(),1000,0,1000,100,0,3.0);
                HistName = "h_mean_pt_vs_phi_";
                HistName += n_mult;
                HistName += "_c";
                HistName += n_w;
                h_mean_pt_vs_phi[n_mult][n_w]  = new TProfile(HistName.Data(),HistName.Data(),100,0,TMath::Pi()*2.0);
            }
        }
    }

    if( eAnalysisNum == 30)
    {
        for(Int_t i = 0; i < 2; i++)
        {
            for(Int_t j = 0; j < 3; j++)
            {
                HistName = "h_phi_pp_corr_";
                HistName += i;
                HistName += "_";
                HistName += j;
                h_phi_pp_corr[i][j] = new TH1D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
                HistName = "h_cos_phi_pp_corr_";
                HistName += i;
                HistName += "_";
                HistName += j;
                h_cos_phi_pp_corr[i][j] = new TH1D(HistName.Data(),HistName.Data(),200,-1,1);
            }
        }


        for(Int_t n_mult = 0; n_mult < N_refmult_bins; n_mult++)
        {
            for(Int_t n_w = 0; n_w < 2; n_w++)
            {
                for(Int_t n_harm = 0; n_harm < n_pp_harm; n_harm++)
                {
                    for(Int_t n_pp = 0; n_pp < 3; n_pp++)
                    {
                        HistName = "p_pp_corr_mult";
                        HistName += n_mult;
                        HistName += "_w";
                        HistName += n_w;
                        HistName += "_harm";
                        HistName += n_harm;
                        HistName += "_pp";
                        HistName += n_pp;
                        p_pp_corr[n_mult][n_w][n_harm][n_pp] = new TProfile2D(HistName.Data(),HistName.Data(),n_pt_bins_pp_corr,0,n_pt_bins_pp_corr,n_pt_bins_pp_corr,0,n_pt_bins_pp_corr);
                    }
                }
            }
        }
    }


    for(Int_t i = 0; i < nMinDCA; i++)
    {
        HistName = "hMinDCA_";
        HistName += i;
        hMinDCA[i] = new TH1F(HistName.Data(),HistName.Data(),400,-200,200);
    }
    for(Int_t i = 0; i < nMinDCA; i++)
    {
        HistName = "hMinDCA_X_";
        HistName += i;
        hMinDCA_X[i] = new TH1F(HistName.Data(),HistName.Data(),400,-200,200);
    }
    for(Int_t i = 0; i < nMinDCA; i++)
    {
        HistName = "hMinDCA_origX_";
        HistName += i;
        hMinDCA_origX[i] = new TH1F(HistName.Data(),HistName.Data(),400,-200,200);
    }
    cMinDCA       = new TCanvas("cMinDCA","cMinDCA",10,10,1200,750);
    cEventDisplay = new TCanvas("cEventDisplay","cEventDisplay",10,10,750,750);
    cVertex_distr = new TCanvas("cVertex_distr","cVertex_distr",10,10,600,750);

    for(Int_t i = 0; i < n_cEventDisplay_array; i++)
    {
        HistName = "cEventDisplay_array_";
        HistName += i;
        cEventDisplay_array[i] = new TCanvas(HistName.Data(),HistName.Data(),10,10,750,750);
        HistName = "cEventDCA_array_";
        HistName += i;
        cEventDCA_array[i] = new TCanvas(HistName.Data(),HistName.Data(),10,10,750,750);
        HistName = "cEvent_MeanX_array_";
        HistName += i;
        cEvent_MeanX_array[i] = new TCanvas(HistName.Data(),HistName.Data(),10,10,750,750);
        HistName = "cEvent_MeanY_array_";
        HistName += i;
        cEvent_MeanY_array[i] = new TCanvas(HistName.Data(),HistName.Data(),10,10,750,750);
        HistName = "cEvent_MeanZ_array_";
        HistName += i;
        cEvent_MeanZ_array[i] = new TCanvas(HistName.Data(),HistName.Data(),10,10,750,750);
    }


    //********************* For 3D graphics **********************************************
    for(Int_t i = 0; i < n_poly_marker_track; i++)
    {
        pTrack[i]   = new TPolyLine3D();
        pTrack[i]   ->SetLineWidth(1);
        pTrack[i]   ->SetLineColor(2);
        pTrack[i]   ->SetLineStyle(1);

        pTrack_extrapolate[i]   = new TPolyLine3D();
        pTrack_extrapolate[i]   ->SetLineWidth(1);
        pTrack_extrapolate[i]   ->SetLineColor(5);
        pTrack_extrapolate[i]   ->SetLineStyle(1);

        pTrack_extrapolate_to_Tof[i]   = new TPolyLine3D();
        pTrack_extrapolate_to_Tof[i]   ->SetLineWidth(1);
        pTrack_extrapolate_to_Tof[i]   ->SetLineColor(5);
        pTrack_extrapolate_to_Tof[i]   ->SetLineStyle(1);

        pTofHit[i] = new TMarker3DBox();
        pTofHit[i] ->SetSize(1.5,1.5,1.5);
        pTofHit[i] ->SetLineColor(4);
        pTofHit[i] ->SetLineStyle(1);
        pTofHit[i] ->SetLineWidth(1);


    }

    XAxis          = new TPolyLine3D(2);
    YAxis          = new TPolyLine3D(2);
    ZAxis          = new TPolyLine3D(2);
    BeamLine       = new TPolyLine3D(2);

    pPrimaryVertex = new TMarker3DBox();
    pPrimaryVertex->SetSize(1.5,1.5,1.5);
    pPrimaryVertex->SetLineColor(6);
    pPrimaryVertex->SetLineStyle(1);
    pPrimaryVertex->SetLineWidth(1);

    pPrimaryVertex2 = new TMarker3DBox();
    pPrimaryVertex2->SetSize(1.5,1.5,1.5);
    pPrimaryVertex2->SetLineColor(7);
    pPrimaryVertex2->SetLineStyle(1);
    pPrimaryVertex2->SetLineWidth(1);

    pPrimaryVertex_old = new TMarker3DBox();
    pPrimaryVertex_old->SetSize(1.5,1.5,1.5);
    pPrimaryVertex_old->SetLineColor(3);
    pPrimaryVertex_old->SetLineStyle(1);
    pPrimaryVertex_old->SetLineWidth(1);

    // Axis points
    XAxis    ->SetPoint(0,0,0,0);
    XAxis    ->SetPoint(1,300,0,0);
    YAxis    ->SetPoint(0,0,0,0);
    YAxis    ->SetPoint(1,0,300,0);
    ZAxis    ->SetPoint(0,0,0,0);
    ZAxis    ->SetPoint(1,0,0,300);
    BeamLine ->SetPoint(0,0,0,-500);
    BeamLine ->SetPoint(1,0,0,500);

    const Int_t n_TPC_points = 50;
    Float_t radius_table[4] = {200,200,3.81,3.81};
    Float_t z_val_table[4]  = {200,-200,200,-200};

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_endcaps[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_endcaps[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        TPC_endcaps[r]->SetLineStyle(0);
        TPC_endcaps[r]->SetLineColor(28);
        TPC_endcaps[r]->SetLineWidth(2);
    }

    BeamLine->SetLineStyle(0);
    BeamLine->SetLineColor(4);
    BeamLine->SetLineWidth(2);

    XAxis->SetLineStyle(0);
    XAxis->SetLineColor(46);
    XAxis->SetLineWidth(2);

    YAxis->SetLineStyle(0);
    YAxis->SetLineColor(3);
    YAxis->SetLineWidth(2);

    ZAxis->SetLineStyle(0);
    ZAxis->SetLineColor(4);
    ZAxis->SetLineWidth(2);
    //*******************************************************************************************


    heta  = new TH1F("heta","heta",200,-1.5,1.5);
    hphi  = new TH1F("hphi","hphi",200,-3.15,3.15);
    hphi_eta_pos  = new TH1F("hphi_eta_pos","hphi_eta_pos",200,-3.15,3.15);
    hphi_eta_neg  = new TH1F("hphi_eta_neg","hphi_eta_neg",200,-3.15,3.15);
    hsign = new TH1F("hsign","hsign",50,-2,2);
    hTPC_dEdx_vs_p = new TH2F("hTPC_dEdx_vs_p","hTPC_dEdx_vs_p",400,-4,4,400,0,0.000015);
    hTOF_vs_p      = new TH2F("hTOF_vs_p","hTOF_vs_p",400,-4,4,400,0,1.5);
    hMass          = new TH1F("hMass","hMass",1000,-2.0,2.0);
    hx_vertex_distr = new TH1F("hx_vertex_distr","hx_vertex_distr",400,-20,20);
    hy_vertex_distr = new TH1F("hy_vertex_distr","hy_vertex_distr",400,-20,20);
    hz_vertex_distr = new TH1F("hz_vertex_distr","hz_vertex_distr",1000,-500,500);
    hEventDCA_new  = new TH1F("hEventDCA_new","hEventDCA_new",200,0,20);
    hEventDCA_old  = new TH1F("hEventDCA_old","hEventDCA_old",200,0,20);
    for(Int_t n = 0; n < n_cEventDisplay_array; n++)
    {
        HistName = "hEvent_fit_points_";
        HistName += n;
        hEvent_fit_points[n] = new TH1F(HistName.Data(),HistName.Data(),100,0,100);
        HistName = "hEvent_dca_";
        HistName += n;
        hEvent_dca[n] = new TH1F(HistName.Data(),HistName.Data(),200,0,20);
        HistName = "hEvent_pt_";
        HistName += n;
        hEvent_pt[n] = new TH1F(HistName.Data(),HistName.Data(),200,0,10);
        HistName = "hEvent_p_";
        HistName += n;
        hEvent_p[n] = new TH1F(HistName.Data(),HistName.Data(),200,0,10);
    }
    hTheta_vs_phi  = new TH2D("hTheta_vs_phi","hTheta_vs_phi",200,-3.2,3.2,200,0.0,3.2);

    hEvent_MeanX = new TH1F("hEvent_MeanX","hEvent_MeanX",600,-100,100);
    hEvent_MeanY = new TH1F("hEvent_MeanY","hEvent_MeanY",600,-100,100);
    hEvent_MeanZ = new TH1F("hEvent_MeanZ","hEvent_MeanZ",800,-400,400);

    hnon_prim_to_prim_ratio     = new TH1D("hnon_prim_to_prim_ratio","hnon_prim_to_prim_ratio",1000,0,5);
    htracks_left_to_right_ratio = new TH1D("htracks_left_to_right_ratio","htracks_left_to_right_ratio",1000,0,5);
    htracks_left_to_right_diff  = new TH1D("htracks_left_to_right_diff","htracks_left_to_right_diff",1000,-5,5);

    //**********************************************************************************************************
    // K* analysis
    hKStar_inv_mass      = new TH1F("hKStar_inv_mass","hKStar_inv_mass",n_inv_mass_bins,start_inv_mass,stop_inv_mass); // 2.5 MeV binning
    hKStar_inv_massB     = new TH1F("hKStar_inv_massB","hKStar_inv_massB",n_inv_mass_bins,start_inv_mass,stop_inv_mass); // 2.5 MeV binning
    hrho_pippim_inv_mass = new TH1F("hrho_pippim_inv_mass","hrho_pippim_inv_mass",n_inv_mass_bins,start_inv_mass,stop_inv_mass); // 2.5 MeV binning

    for(Int_t i = 0; i < nmult_bins; i++)
    {
        for(Int_t j = 0; j < npt_bins; j++)
        {
            for(Int_t k = 0; k < nphi_bins; k++)
            {
                HistName = "h_inv_mass_mult_pt_phi_m";
                HistName += i;
                HistName += "_pt";
                HistName += j;
                HistName += "_phi";
                HistName += k;
                h_inv_mass_mult_pt_phi[i][j][k] = new TH1D(HistName.Data(),HistName.Data(),n_inv_mass_bins,start_inv_mass,stop_inv_mass);
            }
        }
    }

    for(Int_t j = 0; j < npt_bins_D0; j++)
    {
        HistName = "h_inv_mass_D0_pt";
        HistName += "_pt";
        HistName += j;
        h_inv_mass_D0_pt[j] = new TH1D(HistName.Data(),HistName.Data(),1600,-4,4);
    }
    for(Int_t j = 0; j < npt_bins_Sigma1385; j++)
    {
        HistName = "h_inv_mass_Sigma1385_pt";
        HistName += "_pt";
        HistName += j;
        h_inv_mass_Sigma1385_pt[j] = new TH1D(HistName.Data(),HistName.Data(),3200,-4,4);
    }

    for(Int_t k = 0; k < nphi_bins; k++)
    {
        phi_bin_counter[k] = 0;
    }
    hDelta_phi_KStar  = new TH1D("hDelta_phi_KStar","hDelta_phi_KStar",8000,-4,4);
    hEP_Psi_KStar     = new TH1D("hEP_Psi_KStar","hEP_Psi_KStar",400,-4,4);
    hPhi_KStar        = new TH1D("hPhi_KStar","hPhi_KStar",400,-4,4);
    hPhi_vs_Psi       = new TH2D("hPhi_vs_Psi","hPhi_vs_Psi",400,-4,4,400,-4,4);

    for(Int_t j = 0; j < N_pipi_spectra; j++)
    {
        for(Int_t k = 0; k < N_pipi_pt_bins; k++)
        {
            HistName = "h_inv_mass_pipi";
            HistName += "_";
            HistName += j;
            HistName += "_pt";
            HistName += k;
            h_inv_mass_pipi[j][k] = new TH1D(HistName.Data(),HistName.Data(),N_pipi_bins,start_inv_mass_pipi,stop_inv_mass_pipi);
        }
    }
    //**********************************************************************************************************



    //*********************************************************************************************************
    // 7.7 GeV beamtime
    Z_axis_table[0][0]    = -70.0;
    Z_axis_table[0][1]    = -55.0;
    Z_axis_table[0][2]    = -45.0;
    Z_axis_table[0][3]    = -35.0;
    Z_axis_table[0][4]    = -25.0;
    Z_axis_table[0][5]    = -15.0;
    Z_axis_table[0][6]    = -5.0;
    Z_axis_table[0][7]    = 5.0;
    Z_axis_table[0][8]    = 15.0;
    Z_axis_table[0][9]    = 25.0;
    Z_axis_table[0][10]   = 35.0;
    Z_axis_table[0][11]   = 45.0;
    Z_axis_table[0][12]   = 55.0;
    Z_axis_table[0][13]   = 70.0;

    // 11.5 GeV beamtime
    Z_axis_table[1][0]    = -50.0;
    Z_axis_table[1][1]    = -42.0;
    Z_axis_table[1][2]    = -34.0;
    Z_axis_table[1][3]    = -26.0;
    Z_axis_table[1][4]    = -18.0;
    Z_axis_table[1][5]    = -11.0;
    Z_axis_table[1][6]    = -4.0;
    Z_axis_table[1][7]    = 4.0;
    Z_axis_table[1][8]    = 11.0;
    Z_axis_table[1][9]    = 18.0;
    Z_axis_table[1][10]   = 26.0;
    Z_axis_table[1][11]   = 34.0;
    Z_axis_table[1][12]   = 42.0;
    Z_axis_table[1][13]   = 50.0;

    // 39 GeV beamtime
    Z_axis_table[2][0]    = -40.0;
    Z_axis_table[2][1]    = -34.0;
    Z_axis_table[2][2]    = -28.0;
    Z_axis_table[2][3]    = -22.0;
    Z_axis_table[2][4]    = -16.0;
    Z_axis_table[2][5]    = -10.0;
    Z_axis_table[2][6]    = -3.0;
    Z_axis_table[2][7]    = 3.0;
    Z_axis_table[2][8]    = 10.0;
    Z_axis_table[2][9]    = 16.0;
    Z_axis_table[2][10]   = 22.0;
    Z_axis_table[2][11]   = 28.0;
    Z_axis_table[2][12]   = 34.0;
    Z_axis_table[2][13]   = 40.0;

    // 62.4 GeV beamtime
    Z_axis_table[3][0]    = -40.0;
    Z_axis_table[3][1]    = -34.0;
    Z_axis_table[3][2]    = -28.0;
    Z_axis_table[3][3]    = -22.0;
    Z_axis_table[3][4]    = -16.0;
    Z_axis_table[3][5]    = -10.0;
    Z_axis_table[3][6]    = -3.0;
    Z_axis_table[3][7]    = 3.0;
    Z_axis_table[3][8]    = 10.0;
    Z_axis_table[3][9]    = 16.0;
    Z_axis_table[3][10]   = 22.0;
    Z_axis_table[3][11]   = 28.0;
    Z_axis_table[3][12]   = 34.0;
    Z_axis_table[3][13]   = 40.0;

    // 19.6 GeV beamtime
    Z_axis_table[4][0]    = -70.0;
    Z_axis_table[4][1]    = -55.0;
    Z_axis_table[4][2]    = -45.0;
    Z_axis_table[4][3]    = -35.0;
    Z_axis_table[4][4]    = -25.0;
    Z_axis_table[4][5]    = -15.0;
    Z_axis_table[4][6]    = -5.0;
    Z_axis_table[4][7]    = 5.0;
    Z_axis_table[4][8]    = 15.0;
    Z_axis_table[4][9]    = 25.0;
    Z_axis_table[4][10]   = 35.0;
    Z_axis_table[4][11]   = 45.0;
    Z_axis_table[4][12]   = 55.0;
    Z_axis_table[4][13]   = 70.0;

    // 27 GeV beamtime
    Z_axis_table[5][0]    = -70.0;
    Z_axis_table[5][1]    = -55.0;
    Z_axis_table[5][2]    = -45.0;
    Z_axis_table[5][3]    = -35.0;
    Z_axis_table[5][4]    = -25.0;
    Z_axis_table[5][5]    = -15.0;
    Z_axis_table[5][6]    = -5.0;
    Z_axis_table[5][7]    = 5.0;
    Z_axis_table[5][8]    = 15.0;
    Z_axis_table[5][9]    = 25.0;
    Z_axis_table[5][10]   = 35.0;
    Z_axis_table[5][11]   = 45.0;
    Z_axis_table[5][12]   = 55.0;
    Z_axis_table[5][13]   = 70.0;

    // 200 GeV beamtime
    Z_axis_table[6][0]    = -40.0;
    Z_axis_table[6][1]    = -34.0;
    Z_axis_table[6][2]    = -28.0;
    Z_axis_table[6][3]    = -22.0;
    Z_axis_table[6][4]    = -16.0;
    Z_axis_table[6][5]    = -10.0;
    Z_axis_table[6][6]    = -3.0;
    Z_axis_table[6][7]    = 3.0;
    Z_axis_table[6][8]    = 10.0;
    Z_axis_table[6][9]    = 16.0;
    Z_axis_table[6][10]   = 22.0;
    Z_axis_table[6][11]   = 28.0;
    Z_axis_table[6][12]   = 34.0;
    Z_axis_table[6][13]   = 40.0;

    // 14.5 GeV beamtime
    Z_axis_table[7][0]    = -70.0;
    Z_axis_table[7][1]    = -55.0;
    Z_axis_table[7][2]    = -45.0;
    Z_axis_table[7][3]    = -35.0;
    Z_axis_table[7][4]    = -25.0;
    Z_axis_table[7][5]    = -15.0;
    Z_axis_table[7][6]    = -5.0;
    Z_axis_table[7][7]    = 5.0;
    Z_axis_table[7][8]    = 15.0;
    Z_axis_table[7][9]    = 25.0;
    Z_axis_table[7][10]   = 35.0;
    Z_axis_table[7][11]   = 45.0;
    Z_axis_table[7][12]   = 55.0;
    Z_axis_table[7][13]   = 70.0;

    // 200 GeV run14 beamtime
    Z_axis_table[8][0]    = -13.0;
    Z_axis_table[8][1]    = -11.0;
    Z_axis_table[8][2]    = -9.0;
    Z_axis_table[8][3]    = -7.0;
    Z_axis_table[8][4]    = -5.0;
    Z_axis_table[8][5]    = -3.0;
    Z_axis_table[8][6]    = -1.0;
    Z_axis_table[8][7]    = 1.0;
    Z_axis_table[8][8]    = 3.0;
    Z_axis_table[8][9]    = 5.0;
    Z_axis_table[8][10]   = 7.0;
    Z_axis_table[8][11]   = 9.0;
    Z_axis_table[8][12]   = 11.0;
    Z_axis_table[8][13]   = 13.0;


    // Beam time dependend cuts
    fBeamTimeNum     = eBeamTimeNum;
    vertex_z_cut     = vertex_z_array[eBeamTimeNum];
    //*********************************************************************************************************


    if(eBeamTimeNum == 5)
    {
        nsigma_scaling_fac = 1.9; // 27 GeV
    }
    else
    {
        nsigma_scaling_fac = 1.0;
    }



    // K+ TPC dE/dx cut from 0.1 to 0.7 GeV/c
    // To be used only for special purpose
    KaonP_dEdx_cut = new TCutG("KaonP_dEdx_cut",21);
    KaonP_dEdx_cut->SetPoint(0,0.09652918,2.619705e-05);
    KaonP_dEdx_cut->SetPoint(1,0.1342195,1.697274e-05);
    KaonP_dEdx_cut->SetPoint(2,0.1740038,1.127773e-05);
    KaonP_dEdx_cut->SetPoint(3,0.2116942,8.871383e-06);
    KaonP_dEdx_cut->SetPoint(4,0.2745114,6.063984e-06);
    KaonP_dEdx_cut->SetPoint(5,0.3142957,4.860812e-06);
    KaonP_dEdx_cut->SetPoint(6,0.3603617,4.058698e-06);
    KaonP_dEdx_cut->SetPoint(7,0.4001459,3.737852e-06);
    KaonP_dEdx_cut->SetPoint(8,0.4964657,3.497218e-06);
    KaonP_dEdx_cut->SetPoint(9,0.6891053,3.497218e-06);
    KaonP_dEdx_cut->SetPoint(10,0.6765419,3.818064e-06);
    KaonP_dEdx_cut->SetPoint(11,0.6011612,4.219121e-06);
    KaonP_dEdx_cut->SetPoint(12,0.5299683,4.941024e-06);
    KaonP_dEdx_cut->SetPoint(13,0.4462119,5.662926e-06);
    KaonP_dEdx_cut->SetPoint(14,0.375019,7.186943e-06);
    KaonP_dEdx_cut->SetPoint(15,0.328953,8.951595e-06);
    KaonP_dEdx_cut->SetPoint(16,0.2661358,1.256111e-05);
    KaonP_dEdx_cut->SetPoint(17,0.2368211,1.520809e-05);
    KaonP_dEdx_cut->SetPoint(18,0.1991307,2.026141e-05);
    KaonP_dEdx_cut->SetPoint(19,0.1677221,2.611684e-05);
    KaonP_dEdx_cut->SetPoint(20,0.09652918,2.619705e-05);

    fAnalysisNum = eAnalysisNum;

    PicoAlexEventA    = new StPicoAlexEvent();
    PicoAlexEventB    = new StPicoAlexEvent();
    PicoAlexEvent_use = new StPicoAlexEvent();

    for(Int_t i = 0; i < 9; i++)
    {
        dummy_counter_array[i] = 0;
    }

    return 1;
}

Int_t Analysis_Snurf::Make()
{
    cout << "Make started" << endl;
    cout << "******* Start of event loop ********" << endl;
    time (&start_time_ana);


    //----------------------------------------------------------------------------------------------------
    cout << "Define bad run list" << endl;
    Int_t bad_run_numbers[n_bad_run_numbers[eBeamTimeNum]];
    for(Int_t j = 0; j < n_bad_run_numbers[eBeamTimeNum]; j++)
    {
        if(eBeamTimeNum == 0) bad_run_numbers[j] = bad_run_list_7GeV[j];
        if(eBeamTimeNum == 1) bad_run_numbers[j] = bad_run_list_11GeV[j];
        if(eBeamTimeNum == 2) bad_run_numbers[j] = bad_run_list_39GeV[j];
        if(eBeamTimeNum == 3) bad_run_numbers[j] = bad_run_list_62GeV[j];
        if(eBeamTimeNum == 4) bad_run_numbers[j] = bad_run_list_19GeV[j];
        if(eBeamTimeNum == 5) bad_run_numbers[j] = bad_run_list_27GeV[j];
        if(eBeamTimeNum == 6) bad_run_numbers[j] = bad_run_list_200GeV[j];
        if(eBeamTimeNum == 7) bad_run_numbers[j] = bad_run_list_15GeV[j];
        if(eBeamTimeNum == 8) bad_run_numbers[j] = bad_run_list_200GeV_run14[j];
    }
    //----------------------------------------------------------------------------------------------------




    TString HistNameB;

    //Int_t eSE_ME_Flag   = 1; // 0 == same event, 1 == mixed event
    Int_t event_num_sav = 0;
    if(eSE_ME_Flag == 1) cout << "--------------- Mixed event analysis started ---------------" << endl;
    else cout << "--------------- Same event analysis started: " << eSE_ME_Flag << endl;

    Int_t PID_Array_B[N_Ana][N_max_PIDs][N_max_tracks]; // [Ana][PID][number] stores the track numbers for each PID

    Int_t PID_counter_Array_B[N_Ana][N_max_PIDs]; // counter array for PID_Array

    TObjArray* event_ME_array;
    event_ME_array = new TObjArray(N_max_sample);

    //cout << "first_event = " << first_event << ", nEvents = " << nEvents << endl;

    // Event loop
    for(Int_t event_num = first_event; event_num < nEvents; event_num++) // loop over all events
    {
        //cout << "event_num = " << event_num << endl;
        if (!input->GetEntry( event_num )) // take the event -> information is stored in event
            break;  // end of data chunk
        Epico   = picoMaker->picoDst();
        PicoAlexEventA->clearTrackList();
        fillAlexEvent(Epico,PicoAlexEventA);
        PicoAlexEvent_use = PicoAlexEventA;
        //cout << "event_num B = " << event_num << ", event = "  << event << endl;
        if (event_num != 0  &&  event_num % 50 == 0)
            cout << "." << flush;
        if (event_num != 0  &&  event_num % 500 == 0)
        {
            time (&end_time_ana);
            diff_time_ana = difftime (end_time_ana,start_time_ana);
            time (&start_time_ana);
            Float_t event_percent = 100.0*event_num/nEvents;
            cout << " " << event_num << " (" << event_percent << "%) " << "\n" << "==> Processing data, " << flush;
            Long64_t num_events = n_total_entries;
            Double_t events_per_time = 0.0;
            if(diff_time_ana > 0) events_per_time = dummy_counter_loop/diff_time_ana;
            cout << "Total number of events = " << num_events <<
                ", total counter = " << dummy_counter << ", total counter A = " << dummy_counter_A <<
                ", diff time = " << diff_time_ana << ", events/second = " << events_per_time << endl;
            dummy_counter_loop = 0;
        }

        StPicoAlexTrack*   track_A;
        StPicoAlexTrack*   track_B;
        Int_t track_counter   = 0;
        Int_t track_counter_B = 0;
        Int_t   erun_nevents   = PicoAlexEvent_use->numberOfTracks(); // number of tracks in this event
        Float_t eEventVertexX  = PicoAlexEvent_use->primaryVertex().x();
        Float_t eEventVertexY  = PicoAlexEvent_use->primaryVertex().y();
        Float_t eEventVertexZ  = PicoAlexEvent_use->primaryVertex().z();
        //Float_t ebField        = PicoAlexEvent_use->bField();
        Int_t   erefMult       = PicoAlexEvent_use->refMult();
        Int_t   erunId         = PicoAlexEvent_use->runId();
        Float_t eZDCx          = PicoAlexEvent_use->ZDCx(); // ZDC coincidence rate -> luminosity
        Float_t eBBCx          = PicoAlexEvent_use->BBCx(); // BBC coincidence rate -> luminosity

        //----------------------------------------------------------------------------------------------------
        // reference multiplicity correction
        // ******* IMPORTANT ***********
        // Call initEvent(const UShort_t RefMult, const Double_t z) function
        // event-by-event at the beginning before using any other functions
        refmultCorrUtil->init(erunId);
        refmultCorrUtil->initEvent(erefMult, eEventVertexZ, eZDCx);

        // Get centrality bins
        //   see StRefMultCorr.h for the definition of centrality bins
        erefMult_bin16   = refmultCorrUtil->getCentralityBin16();
        erefMult_bin     = refmultCorrUtil->getCentralityBin9();

        // erefMult_bin (9)
        // 0 = 70-80%
        // 1 = 60-70%
        // 2 = 50-60%
        // 3 = 40-50%
        // 4 = 30-40%
        // 5 = 20-30%
        // 6 = 10-20%
        // 7 = 5-10%
        // 8 = 0-5%

        //cout << "erefMult_bin = " << erefMult_bin << endl;

        // Re-weighting corrections for peripheral bins
        //Double_t reweight = refmultCorrUtil->getWeight();

        // Not really necessary for your study but if you want to see the corrected refmult distribution
        // Corrected refmult (z-vertex dependent correction)
        //  NOTE: type should be double or float, not integer
        //Double_t refmultCor = refmultCorrUtil->getRefMultCorr();

        //cout << "erefMult = " << erefMult << ", refmultCor = " << refmultCor << ", reweight = " << reweight
        //    << ", z = " << eEventVertexZ << ", cent16 = " << cent16 << ", cent9 = " << cent9 << ", erunId = " << erunId << endl;
        //----------------------------------------------------------------------------------------------------




        Int_t Last_event = 0;

        if(
           event_num == (nEvents - 1) // last event -> do final mixing with everything which is left over
          )
        {
            Last_event = 1;
        }

        //cout << "event_num = " << event_num << ", eEventVertexX = " << eEventVertexX
        //    << ", eEventVertexY = " << eEventVertexY << ", eEventVertexZ = " << eEventVertexZ << ", N_tracks = " << erun_nevents << endl;

        // Loop only over good events
        if(
           ((eEventVertexX*eEventVertexX + eEventVertexY*eEventVertexY) < Event_radius_cut*Event_radius_cut
            && fabs(eEventVertexZ) < vertex_z_cut)
           || (Last_event  == 1 && eSE_ME_Flag == 1)
           || eAnalysisNum == 22 // refMult QA analysis
           || eAnalysisNum == 9
          )
        {
            //cout << "Accepted" << endl;
            // Calculate the z-axis bin of the event
            Int_t ezaxis_bin = -1;
            if(eEventVertexZ >= Z_axis_table[eBeamTimeNum][0])
            {
                for(Int_t r = 0; r < N_z_axis_bins-1; r++)
                {
                    if(eEventVertexZ >= Z_axis_table[eBeamTimeNum][r] && eEventVertexZ < Z_axis_table[eBeamTimeNum][r+1])
                    {
                        ezaxis_bin = r;
                        break;
                    }
                }
            }

            //cout << "erefMult = " << erefMult << ", erefMult_bin = " << erefMult_bin << ", eEventVertexZ = " << eEventVertexZ << ", z-axis bin = " << ezaxis_bin << endl;
            //cout << "event_num = " << event_num << ", erefMult_bin = " << erefMult_bin << ", ezaxis_bin = "
            //    << ezaxis_bin << ", counter = " << Sample_Counter[erefMult_bin][ezaxis_bin] << endl;


            //***********************************************************************************************************************************************************
            // Event plane calculation for mixed event buffer
            Int_t eep_bin = 0;
            if(eSE_ME_Flag == 1)
            {
                Int_t PID_Array_ME[N_Ana][N_max_PIDs][N_max_tracks]; // [Ana][PID][number] stores the track numbers for each PID
                memset(&PID_Array_ME[0][0][0],0,sizeof(Int_t)*N_Ana*N_max_PIDs*N_max_tracks);
                Int_t PID_counter_Array_ME[N_Ana][N_max_PIDs]; // counter array for PID_Array
                memset(&PID_counter_Array_ME[0][0],0,sizeof(Int_t)*N_Ana*N_max_PIDs);
                //Track loop
                //cout << "********** NEW EVENT " << event_num << ", with " << erun_nevents << " tracks **********" << endl;
                Int_t track_counter_ME = 0;
                for(UShort_t i = 0; i < erun_nevents; ++i) // loop over all tracks of the actual event
                {
                    track_A = PicoAlexEvent_use->track( i ); // take the track
                    fPID(track_A); // Particle identification -> writes multiple PIDs into PIDs_track[Ana][i] Array (global definition) Protons == 0, Kaons == 1, Pions == 2, Electrons == 3
                    // In the following loops the track numbers and for mixed event the event numbers are stored in arrays
                    for(Int_t ipid = 0; ipid < N_max_PIDs_per_track; ipid++) // loop over maximum number of possible PIDs per track (usually 4)
                    {
                        for(Int_t iana = 0; iana < N_Ana; iana++) // Loop over all analyses
                        {
                            if(PIDs_track[iana][ipid] != 0 && PID_counter_Array_ME[iana][PIDs_track[iana][ipid]] < N_max_tracks) // if the PID of this species and analysis is not 0
                            {
                                PID_Array_ME[iana][PIDs_track[iana][ipid]][PID_counter_Array_ME[iana][PIDs_track[iana][ipid]]] = i; // Store the track number for this PID and analysis into the PID_Array
                                PID_counter_Array_ME[iana][PIDs_track[iana][ipid]]++; // Increase the counter for this PID and analysis
                            }
                        }
                    }
                    track_counter_ME++;
                }

                Float_t EP_Qx_ME = 0.0; // x value of event plane vector
                Float_t EP_Qy_ME = 0.0; // y value of event plane vector
                EventPlane_analysis(PID_counter_Array_ME,PID_Array_ME,PicoAlexEvent_use,14,1,track_counter_ME,event_num,EP_Qx_ME,EP_Qy_ME,EP_Qx_eta_pos,EP_Qy_eta_pos,EP_Qx_eta_neg,EP_Qy_eta_neg);
                Float_t phi_event_plane_ME = calc_phi_event_plane_2nd(EP_Qx_ME,EP_Qy_ME);
                eep_bin = (Int_t)((phi_event_plane_ME - phi_EP_start)/phi_EP_step);
                //cout << "Event number = " << event_num << ", event plane angle = " << phi_event_plane_ME << ", eep_bin = " << eep_bin << endl;
            }
            //***********************************************************************************************************************************************************


            Int_t Start_mixing = 0;
            if(eSE_ME_Flag == 1 && erefMult_bin >= 0 && ezaxis_bin >= 0 && eep_bin >= 0 && Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin] < N_max_sample_array[erefMult_bin][ezaxis_bin][eep_bin]) // if mixed event analysis is on and the sample has not reached the maximum number...
            {
                Sample_Event_Num[erefMult_bin][ezaxis_bin][eep_bin][Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin]] = event_num; // stores the event number for each bin [centrality bin][z-axis bin]
                Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin]++; // counts the number of events for each bin [centrality bin][z-axis bin]
            }


            if(eSE_ME_Flag == 1 && erefMult_bin >= 0 && ezaxis_bin >= 0 && eep_bin >= 0 && Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin] >= N_max_sample_array[erefMult_bin][ezaxis_bin][eep_bin]) // sample has reached maximum number -> start mixing of this sample [centrality bin][z-axis bin]
            {
                Start_mixing = 1;
                event_num_sav = event_num;  // event number has to be saved to later restart at the same event
            }

            Int_t N_refmult_bins_use = 1; // usually 1, only for last event it will be changed to maximum number
            Int_t N_z_axis_bins_use  = 1;
            Int_t N_EP_bins_use      = 1;

            if(
               Last_event == 1 // last event -> do final mixing with everything which is left
               && eSE_ME_Flag == 1
              )
            {
                Start_mixing = 1;
                N_refmult_bins_use = N_refmult_bins;
                N_z_axis_bins_use  = N_z_axis_bins;
                N_EP_bins_use      = N_EP_bins;

                //cout << "Start event mixing at last event: " << event_num << endl;
                //cout << "Centrality bin = " << erefMult_bin << ", z-axis bin = " << ezaxis_bin << endl;

                for(Int_t e = 0; e < N_refmult_bins; e++)
                {
                    cout << "------------ ref. mult. bin = " << e << " ----------------" << endl;
                    cout << "[" << flush;
                    for(Int_t f = 0; f < N_z_axis_bins-1; f++)
                    {
                        cout << "------------ z-axis bin = " << f << " ----------------" << endl;
                        cout << "[" << flush;
                        for(Int_t g = 0; g < N_EP_bins; g++)
                        {
                            cout << Sample_Counter[e][f][g] << ", " << flush;
                        }
                        cout << "]" << endl;
                    }
                    cout << "]" << endl;
                    cout << "----------------------------------------------------------" << endl;
                }
                cout << "-----------------------------------------------------------------------------------------" << endl;
            }


            // The following three loops are only activated above 1 for the last event and mixed event!
            for(Int_t lref_mult = 0; lref_mult < N_refmult_bins_use; lref_mult++)  // loop usually does only loop over one ref mult
            {
                for(Int_t lzaxis = 0; lzaxis < N_z_axis_bins_use; lzaxis++)
                {
                    for(Int_t lep_bin = 0; lep_bin < N_EP_bins_use; lep_bin++)
                    {
                        //cout << "lref_mult = " << lref_mult << ", lzaxis = " << lzaxis << endl;
                        if(Last_event == 1 && eSE_ME_Flag == 1)
                        {
                            erefMult_bin = lref_mult;
                            ezaxis_bin   = lzaxis;
                            eep_bin      = lep_bin;
                            cout << "Final mixing: lref_mult = " << lref_mult << ", lzaxis = " << lzaxis << ", lep_bin = " << lep_bin << endl;
                        }

                        if(
                           Start_mixing == 1
                           || eSE_ME_Flag == 0
                          )
                        {

                            Int_t em_number_A = 1; // for same event analysis == 1 -> only one loop, for mixed event -> sample counter
                            Int_t em_number_B = 1; // for same event analysis == 1 -> only one loop, for mixed event -> sample counter

                            if(
                               Start_mixing == 1
                              )
                            {
                                em_number_A = Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin]; // number of entries for these bins
                                em_number_B = Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin]; // number of entries for these bins
                                Sample_Counter[erefMult_bin][ezaxis_bin][eep_bin] = 0; // reset the sample counter

                                // Fill the event_ME_array with the events from this buffer
                                Int_t em_number_buffer = em_number_A;
                                //cout << "Fill buffer" << endl;
                                for(Int_t nc_max_sample = 0; nc_max_sample < em_number_buffer; nc_max_sample++)
                                {
                                    HistNameB = "event_ME_array_";
                                    HistNameB += nc_max_sample;
                                    input->GetEntry( Sample_Event_Num[erefMult_bin][ezaxis_bin][eep_bin][nc_max_sample] );
                                    StPicoDst* Epico_ME = picoMaker->picoDst();
                                    PicoAlexEventA->clearTrackList();
                                    fillAlexEvent(Epico_ME,PicoAlexEventA);
                                    event_ME_array->AddAt((StPicoAlexEvent*)PicoAlexEventA->Clone(HistNameB.Data()),nc_max_sample);
                                    PicoAlexEventA->clearTrackList();
                                }
                                event_ME_array->SetOwner(kTRUE); // TObjArray is now the owner of the content
                            }


                            for(Int_t em_A = 0; em_A < em_number_A; em_A++)
                            {

                                if(
                                   Start_mixing == 1
                                  )
                                {
                                    //input->GetEntry( Sample_Event_Num[erefMult_bin][ezaxis_bin][eep_bin][em_A] );
                                    //event_use = event;
                                    PicoAlexEvent_use = (StPicoAlexEvent*)(event_ME_array->At(em_A));
                                    eEvent_numA    = Sample_Event_Num[erefMult_bin][ezaxis_bin][eep_bin][em_A];
                                    erun_nevents   = PicoAlexEvent_use->numberOfTracks(); // number of tracks in this event
                                    //cout << "(B) Start event mixing at event number: " << event_num << ", refMult bin = " << erefMult_bin << ", z-axis bin = " << ezaxis_bin << ", number of tracks = " << erun_nevents << endl;
                                    //cout << "em_A = " << em_A << ", old tracks = " << erun_nevents << ", new tracks = " << (&event_ME_array[em_A])->numberOfTracks() << endl;
                                }
                                Int_t PID_Array[N_Ana][N_max_PIDs][N_max_tracks]; // [Ana][PID][number] stores the track numbers for each PID
                                memset(&PID_Array[0][0][0],0,sizeof(Int_t)*N_Ana*N_max_PIDs*N_max_tracks);
                                Int_t PID_counter_Array[N_Ana][N_max_PIDs]; // counter array for PID_Array
                                memset(&PID_counter_Array[0][0],0,sizeof(Int_t)*N_Ana*N_max_PIDs);
                                //Track loop
                                //cout << "********** NEW EVENT " << event_num << ", with " << erun_nevents << " tracks **********" << endl;

                                n_primaries     = 0;
                                n_non_primaries = 0;
                                n_tofmatch_prim = 0;
                                Int_t n_tracks_left   = 0;
                                Int_t n_tracks_right  = 0;

                                for(UShort_t i = 0; i < erun_nevents; ++i) // loop over all tracks of the actual event
                                {
                                    track_A = PicoAlexEvent_use->track( i ); // take the track
                                    fPID(track_A); // Particle identification -> writes multiple PIDs into PIDs_track[Ana][i] Array (global definition) Protons == 0, Kaons == 1, Pions == 2, Electrons == 3
                                    // In the following loops the track numbers and for mixed event the event numbers are stored in arrays
                                    for(Int_t ipid = 0; ipid < N_max_PIDs_per_track; ipid++) // loop over maximum number of possible PIDs per track (usually 4)
                                    {
                                        for(Int_t iana = 0; iana < N_Ana; iana++) // Loop over all analyses
                                        {
                                            if(PIDs_track[iana][ipid] != 0 && PID_counter_Array[iana][PIDs_track[iana][ipid]] < N_max_tracks) // if the PID of this species and analysis is not 0
                                            {
                                                PID_Array[iana][PIDs_track[iana][ipid]][PID_counter_Array[iana][PIDs_track[iana][ipid]]] = i; // Store the track number for this PID and analysis into the PID_Array
                                                PID_counter_Array[iana][PIDs_track[iana][ipid]]++; // Increase the counter for this PID and analysis
                                            }
                                        }
                                    }

                                    //*************************************************************************************
                                    // QA

                                    //Float_t Polarity   = (nHitsFit >= 0) ? 1.0 : -1.0; // charge of the particle
                                    Float_t Polarity = static_cast<Float_t>(track_A->charge());
                                    StThreeVectorF p(track_A->pMom());  // primary momentum
                                    Float_t Momentum = p.mag(); // momentum
                                    Float_t Beta     = static_cast<Float_t>(track_A->btofBeta()); // beta


                                    Float_t Mass2 = -100.0;
                                    // calculate mass2
                                    if(track_A->btofMatchFlag() > 0 && track_A->btof() != 0 && Beta != 0)
                                    {
                                        Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
                                    }

                                    //Float_t nSigmaEl   = track_A->nSigmaElectron();
                                    Float_t nSigmaPion = track_A->nSigmaPion();
                                    //Float_t nSigmaKaon = track_A->nSigmaKaon();
                                    //Float_t nSigmaP    = track_A->nSigmaProton();


                                    Float_t dca        = track_A->dca();
                                    Float_t nHitsPoss  = track_A->nHitsMax();
                                    Float_t nHitsFit   = track_A->nHitsFit();

                                    Float_t eta = p.pseudoRapidity();

                                    if(fabs(nSigmaPion) < 2.0 && dca < 2.0 && Momentum > 0.0 && nHitsPoss > 0.0 && eta < 0.9 && (nHitsFit/nHitsPoss) > 0.52)
                                    {
                                        h_nhitsFit ->Fill(p.perp(),nHitsFit);
                                        if(nHitsFit > 39 && p.perp() > 0.2 && p.perp() < 0.6) h_eta_phi[0]->Fill(p.phi(),eta);
                                        if(nHitsFit < 36 && nHitsFit > 30  && p.perp() > 0.2 && p.perp() < 0.6)  h_eta_phi[1]->Fill(p.phi(),eta);
                                        if(nHitsFit > 39 && p.perp() > 0.6 && p.perp() < 1.0) h_eta_phi[2]->Fill(p.phi(),eta);
                                        if(nHitsFit < 36 && nHitsFit > 30  && p.perp() > 0.6 && p.perp() < 1.0)  h_eta_phi[3]->Fill(p.phi(),eta);
                                        if(nHitsFit > 39 && p.perp() > 1.0 && p.perp() < 1.5) h_eta_phi[4]->Fill(p.phi(),eta);
                                        if(nHitsFit < 36 && nHitsFit > 30  && p.perp() > 1.0 && p.perp() < 1.5)  h_eta_phi[5]->Fill(p.phi(),eta);
                                        if(nHitsFit > 39 && p.perp() > 1.5 && p.perp() < 2.0) h_eta_phi[6]->Fill(p.phi(),eta);
                                        if(nHitsFit < 36 && nHitsFit > 30  && p.perp() > 1.5 && p.perp() < 2.0)  h_eta_phi[7]->Fill(p.phi(),eta);
                                        if(nSigmaPion < 0.2)
                                        {
                                            if(nHitsFit > 39 && p.perp() > 1.0 && p.perp() < 1.5) h_eta_phi[8]->Fill(p.phi(),eta);
                                            if(nHitsFit < 36 && nHitsFit > 30  && p.perp() > 1.0 && p.perp() < 1.5)  h_eta_phi[9]->Fill(p.phi(),eta);
                                        }
                                        if(Mass2 < 0.15 && Mass2 > -0.2)
                                        {
                                            h_nhitsFit_TOF->Fill(p.perp(),nHitsFit);
                                        }
                                    }

                                    StPhysicalHelixD helixA = StPhysicalHelixD(track_A->gMom(),track_A->origin(),PicoAlexEvent_use->bField(),track_A->charge());

                                    StThreeVectorD vectoratsA,vectoratsB,vector_diff,vectorprim,primdirA,vectornewA;
                                    vectorprim.set(eEventVertexX,eEventVertexY,eEventVertexZ);
                                    vectoratsA     = helixA.at(0.0);
                                    vectoratsB     = helixA.at(0.0+0.1);
                                    vector_diff = (vectoratsB-vectoratsA)/(vectoratsB-vectoratsA).mag();
                                    Double_t vector_z = vector_diff.z();
                                    primdirA = helixA.cat(helixA.pathLength(vectorprim)); //

                                    Float_t TPCdEdx    = track_A->dEdx();


                                    if(Start_mixing != 1)
                                    {
                                        hTPC_dEdx_vs_p  ->Fill(Polarity*Momentum,TPCdEdx);
                                        hTOF_vs_p       ->Fill(Polarity*Momentum,Beta);
                                        if(Mass2 > 0.0)
                                        {
                                            Float_t Mass = sqrt(Mass2);
                                            hMass  ->Fill(Polarity*Mass);
                                        }
                                    }

                                    //Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm

                                    if(
                                       nHitsPoss               > 0.0
                                       && nHitsFit             >= 15
                                       && (nHitsFit/nHitsPoss) > 0.52
                                       && Momentum             > 0.1
                                       && Momentum             < 10.0
                                       && (eEventVertexX*eEventVertexX + eEventVertexY*eEventVertexY) < Event_radius_cut*Event_radius_cut
                                       //&& fabs(eEventVertexZ) < z_axis_cut
                                      )
                                    {
                                        Double_t radius = sqrt(primdirA.x()*primdirA.x()+primdirA.y()*primdirA.y());
                                        //Double_t Theta  = TMath::ATan2(radius,primdirA.z());
                                        hTheta_vs_phi   ->Fill(primdirA.phi(),TMath::ATan2(radius,primdirA.z()));
                                        //if(dca < 3.0)
                                        {
                                            if(vector_z < 0.0) n_tracks_left++;
                                            if(vector_z > 0.0) n_tracks_right++;
                                        }
                                        if(dca < 3.0)
                                        {
                                            n_primaries++;
                                            if(track_A->btofMatchFlag() > 0 && track_A->btof() != 0 && Beta != 0)
                                            {
                                                n_tofmatch_prim++;
                                            }
                                        }
                                        if(dca >= 3.0)
                                        {
                                            n_non_primaries++;
                                        }
                                    }
                                    //*************************************************************************************


                                    track_counter++;
                                }


                                // Calculate event plane for event A -> event B is only calculated later for mixed event analysis
                                EP_Qx = 0.0; // x value of event plane vector
                                EP_Qy = 0.0; // y value of event plane vector
                                EventPlane_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num,EP_Qx,EP_Qy,EP_Qx_eta_pos,EP_Qy_eta_pos,EP_Qx_eta_neg,EP_Qy_eta_neg);
                                //cout << "EP_Qx = " << EP_Qx << ", EP_Qy = " << EP_Qy << endl;


                                non_prim_to_prim_ratio = 0.0;
                                if(n_primaries > 0 && n_non_primaries > 0)
                                {
                                    non_prim_to_prim_ratio = ((Double_t)n_non_primaries)/((Double_t)n_primaries);
                                    hnon_prim_to_prim_ratio      ->Fill(non_prim_to_prim_ratio);
                                }
                                tracks_left_to_right_ratio = 0.0;
                                tracks_left_to_right_diff  = 0.0;
                                if(n_tracks_right > 0.0 && n_tracks_left > 0.0)
                                {
                                    tracks_left_to_right_ratio = ((Double_t)n_tracks_left)/((Double_t)n_tracks_right);
                                    tracks_left_to_right_diff  = ((Double_t)(n_tracks_left-n_tracks_right))/((Double_t)(n_tracks_left+n_tracks_right));
                                    htracks_left_to_right_ratio  ->Fill(tracks_left_to_right_ratio);
                                    htracks_left_to_right_diff   ->Fill(tracks_left_to_right_diff);
                                }

                                //cout << "n_tracks_left = " << n_tracks_left << ", n_tracks_right = " << n_tracks_right << ", erefMult = " << erefMult << endl;

                                //Int_t use_mixing_rate_number = mixing_rate_number;
                                //if(Start_mixing != 1) // no mixing
                                //{
                                //    use_mixing_rate_number = 1; // only one loop of event B -> second loop will be ignored
                                //}

                                //for(Int_t em_B = em_A+1; em_B < em_A+1+use_mixing_rate_number; em_B++)  // ok!... but inefficient and switching events is not used
                                for(Int_t em_B = 0; em_B < em_number_B; em_B++)  // The loop from 0 takes care that event A and B are switched, e.g. A = 1 && B = 2, later A = 2 && B = 1
                                {
                                    if(Start_mixing == 1 && em_A != em_B)
                                    {
                                        if(em_B < N_max_sample_array[erefMult_bin][ezaxis_bin][eep_bin])
                                        {
                                            //input_B->GetEntry( Sample_Event_Num[erefMult_bin][ezaxis_bin][eep_bin][em_B] );
                                            //event_B_use = event_B;
                                            //cout << "Before em_A = " << em_A << ", em_B = " << em_B << ", erun_nevents = " << erun_nevents
                                            //    << ", event_use = " << event_use << endl;
                                            PicoAlexEventB = (StPicoAlexEvent*)(event_ME_array->At(em_B));
                                            eEvent_numB    = Sample_Event_Num[erefMult_bin][ezaxis_bin][eep_bin][em_B];
                                            Int_t erun_nevents_B   = PicoAlexEventB->numberOfTracks(); // number of tracks in this event
                                            //cout << "After em_A = " << em_A << ", em_B = " << em_B << ", erun_nevents = " << Epico_use->numberOfTracks() << ", erun_nevents_B = " << Epico_B_use->numberOfTracks()
                                            //    << ", event_use = " << event_use << ", event_B_use = " << event_B_use << ", Epico_use = "
                                            //    << Epico_use << ", Epico_B_use = " << Epico_B_use
                                            //    << ", Epico_use(B) = " << event_ME_array->At(em_A) << ", Epico_B_use(B) = " << event_ME_array->At(em_B) << endl;

                                            memset(&PID_Array_B[0][0][0],0,sizeof(Int_t)*N_Ana*N_max_PIDs*N_max_tracks);
                                            memset(&PID_counter_Array_B[0][0],0,sizeof(Int_t)*N_Ana*N_max_PIDs);
                                            //Track loop
                                            //cout << "********** NEW EVENT " << event_num << ", with " << erun_nevents << " tracks **********" << endl;
                                            //cout << "em_A = " << em_A << ", em_B = " << em_B << ", erefMult_bin = " << erefMult_bin
                                            //    << ", ezaxis_bin = " << ezaxis_bin << ", N_sample = " << em_number_A << ", eEvent_numA = " << eEvent_numA << ", eEvent_numB = " << eEvent_numB
                                            //    << ", total_counter = " << total_counter << endl;
                                            for(UShort_t i = 0; i < erun_nevents_B; ++i) // loop over all tracks of the actual event
                                            {
                                                track_B = PicoAlexEventB->track( i ); // take the track
                                                fPID(track_B); // Particle identification -> writes multiple PIDs into PIDs_track[Ana][i] Array (global definition) Protons == 0, Kaons == 1, Pions == 2, Electrons == 3
                                                // In the following loops the track numbers and for mixed event the event numbers are stored in arrays
                                                for(Int_t ipid = 0; ipid < N_max_PIDs_per_track; ipid++) // loop over maximum number of possible PIDs per track (usually 4)
                                                {
                                                    for(Int_t iana = 0; iana < N_Ana; iana++) // Loop over all analyses
                                                    {
                                                        if(PIDs_track[iana][ipid] != 0 && PID_counter_Array_B[iana][PIDs_track[iana][ipid]] < N_max_tracks) // if the PID of this species and analysis is not 0
                                                        {
                                                            PID_Array_B[iana][PIDs_track[iana][ipid]][PID_counter_Array_B[iana][PIDs_track[iana][ipid]]] = i; // Store the track number for this PID and analysis into the PID_Array
                                                            PID_counter_Array_B[iana][PIDs_track[iana][ipid]]++; // Increase the counter for this PID and analysis
                                                        }
                                                    }
                                                }

                                                track_counter_B++;
                                            }
                                        }
                                        else break;
                                    }


                                    // Correlation analysis
                                    //**********************************************************************************************************************************
                                    // For exp data
                                    //cout << "Start of phi analysis, track_counter = " << track_counter
                                    //    << ", N(K+) = " << PID_counter_Array[0][11] << ", N(K-) = " << PID_counter_Array[0][12] << endl;
                                    x_mean = -999.0;
                                    y_mean = -999.0;
                                    z_mean = -999.0;
                                    seed_number = ran_seed_number.Integer(10000);
                                    //Int_t Vertexout = Vertex_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,14,1,track_counter,event_num,x_mean,y_mean,z_mean,1);
                                    //cout << "******************* Start of analysis with " << PID_counter_Array[1][14] << " tracks *******************" << endl;
                                    //Int_t test_func_out = test_func(PID_counter_Array,PID_Array,event_use,14,1,track_counter,event_num);
                                    //Int_t Hadronout   = Hadron_analysis(PID_counter_Array,PID_Array,event_use,14,1,track_counter,event_num);

                                    if(eSE_ME_Flag == 1 && Start_mixing == 1 && em_A != em_B)  // mixed event
                                    {
                                        EP_Qx_B = 0.0; // x value of event plane vector
                                        EP_Qy_B = 0.0; // y value of event plane vector
                                        EventPlane_analysis(PID_counter_Array_B,PID_Array_B,PicoAlexEventB,14,1,track_counter_B,event_num,EP_Qx_B,EP_Qy_B,EP_Qx_eta_pos_B,EP_Qy_eta_pos_B,EP_Qx_eta_neg_B,EP_Qy_eta_neg_B);
                                        Lambda_ME_counter = 0;

                                        if(eAnalysisNum == 1)
                                        {
                                            PhiKPKM_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,12,11,2,track_counter,track_counter_B,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 25)
                                        {
                                            JPsi_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,3,2,2,track_counter,track_counter_B,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 2)
                                        {
                                            Lambda_ppim_analysis(PID_counter_Array,PID_counter_Array_B,PID_counter_Array_B,PID_Array,PID_Array_B,PID_Array_B,
                                                                 PicoAlexEvent_use,PicoAlexEventB,PicoAlexEventB,14,9,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 24)
                                        {
                                            ThetaPlus_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array,PID_Array_B,
                                                                 PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEventB,8,9,14,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 14)
                                        {
                                            Lambda_ppim_analysis(PID_counter_Array,PID_counter_Array_B,PID_counter_Array_B,PID_Array,PID_Array_B,PID_Array_B,
                                                                 PicoAlexEvent_use,PicoAlexEventB,PicoAlexEventB,15,8,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 3)
                                        {
                                            K0S_pippim_analysis(PID_counter_Array,PID_counter_Array_B,PID_counter_Array_B,PID_Array,PID_Array_B,PID_Array_B,
                                                                 PicoAlexEvent_use,PicoAlexEventB,PicoAlexEventB,8,9,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 400)
                                        {
                                            D0_analysis(PID_counter_Array,PID_counter_Array_B,PID_counter_Array_B,PID_Array,PID_Array_B,PID_Array_B,
                                                        PicoAlexEvent_use,PicoAlexEventB,PicoAlexEventB,12,8,12,2,eSE_ME_Flag);
                                            D0_analysis(PID_counter_Array,PID_counter_Array_B,PID_counter_Array_B,PID_Array,PID_Array_B,PID_Array_B,
                                                        PicoAlexEvent_use,PicoAlexEventB,PicoAlexEventB,11,9,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 20)
                                        {
                                            Merged_Omega_Xi_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                                                     PicoAlexEvent_use,PicoAlexEventB,14,9,2,eSE_ME_Flag);
                                            Merged_Omega_Xi_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                                                     PicoAlexEvent_use,PicoAlexEventB,15,8,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 351)
                                        {
                                            //LambdaCPlus_V0_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                            //                        PicoAlexEvent_use,PicoAlexEventB,14,9,8,2,eSE_ME_Flag);
                                            LambdaCPlus_V0_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                                                    PicoAlexEvent_use,PicoAlexEventB,8,9,14,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 200)
                                        {
                                            Omega2250_analysis(PID_counter_Array,PID_counter_Array_B,
                                                               PID_Array,PID_Array_B,
                                                               PicoAlexEvent_use,PicoAlexEventB,
                                                               14,9,9,8,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 201) // mixed event
                                        {
#if 1
                                            // standard analysis
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array_B,
                                                             PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,
                                                             14,9,9,14,9,2,eSE_ME_Flag,0,
                                                             0
                                                            );
#endif

#if 0
                                            // "Omega" instead of Xi
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array_B,
                                                             PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,
                                                             14,9,12,14,9,2,eSE_ME_Flag,0);
#endif

#if 0
                                            // anti-Lambda + Xi-
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array_B,
                                                             PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,
                                                             14,9,9,15,8,2,eSE_ME_Flag,0);
#endif

#if 0
                                            // K0S + K0S
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array_B,
                                                             PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,
                                                             8,9,9,8,9,2,eSE_ME_Flag,0);
#endif

#if 0
                                            d4s_analysis(PID_counter_Array,PID_counter_Array_B,
                                                         PID_Array,PID_Array_B,
                                                         PicoAlexEvent_use,PicoAlexEventB,
                                                         14,9,9,14,9,2,eSE_ME_Flag,0);
                                            //d4s_analysis(PID_counter_Array,PID_counter_Array_B,
                                            //                   PID_Array,PID_Array_B,
                                            //                   PicoAlexEvent_use,PicoAlexEventB,
                                            //                   14,9,9,15,8,2,eSE_ME_Flag,0);
#endif
                                        }
                                        if(eAnalysisNum == 203) // mixed event
                                        {
                                            // standard analysis
                                        }
                                        if(eAnalysisNum == 15)
                                        {
                                            rho_analysis(PID_counter_Array,PID_counter_Array_B,PID_Array,PID_Array_B,
                                                             PicoAlexEvent_use,PicoAlexEventB,8,9,2,track_counter,track_counter_B,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 23)
                                        {
                                            //K_Pi_three_analysis(PID_counter_Array,PID_counter_Array_B,PID_counter_Array_B,PID_Array,PID_Array_B,PID_Array_B,
                                            //                    PicoAlexEvent_use,PicoAlexEventB,PicoAlexEventB,8,8,9,9,2,track_counter,track_counter,event_num,eSE_ME_Flag,eAnalysisNum);
                                        }
                                    }

                                    if(eSE_ME_Flag == 0) // same event
                                    {

                                        //cout << "--------------- em_A = " << em_A << ", em_B = " << em_B << " ----------------" << endl;
                                        if(eAnalysisNum == 1)
                                        {
                                            //cout << "Start Phi analysis, event number = " << event_num << ", lref_mult = " << lref_mult << ", lzaxis = " << lzaxis << ", lep_bin = " << lep_bin
                                            //    << ", eEvent_numA = " << eEvent_numA << ", ntracks = " << erun_nevents
                                            //    << ", EP_Qx = " << EP_Qx << ", EP_Qy = " << EP_Qy << endl;
                                            PhiKPKM_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                             PicoAlexEvent_use,PicoAlexEvent_use,12,11,2,track_counter,track_counter,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 35)
                                        {
                                            LambdaCplus_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                             PicoAlexEvent_use,PicoAlexEvent_use,14,8,12,2,track_counter,track_counter,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 25)
                                        {
                                            JPsi_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                             PicoAlexEvent_use,PicoAlexEvent_use,3,2,2,track_counter,track_counter,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 2)
                                        {
                                            Lambda_ppim_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                 PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,14,9,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 14)
                                        {
                                            Lambda_ppim_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                 PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,15,8,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 102)
                                        {
                                            Lambda_ppim_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                 PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,14,9,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 24)
                                        {
                                            ThetaPlus_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                 PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,8,9,14,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 3)
                                        {
                                            K0S_pippim_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                 PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,8,9,12,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 400)
                                        {
                                            D0_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,12,8,12,2,eSE_ME_Flag);
                                            D0_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                                                PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,11,9,12,2,eSE_ME_Flag);
                                   
                                        }
                                        if(eAnalysisNum == 22)
                                        {
                                            RefMult_QA_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1);
                                        }
                                        if(eAnalysisNum == 360)
                                        {
                                            Ach_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1);
                                        }
                                        if(eAnalysisNum == 8)
                                        {
                                            Hadron_v2_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num);
                                        }

                                        if(eAnalysisNum == 88)
                                        {
                                            Int_t flag_good_run = 1;
                                            for(Int_t bad_run = 0; bad_run < n_bad_run_numbers[eBeamTimeNum]; bad_run++)
                                            {
                                                if(bad_run_numbers[bad_run] == (Int_t)erunId)
                                                {
                                                    flag_good_run = 0;
                                                }
                                            }
                                            if(
                                               flag_good_run == 1
                                              )
                                            {
                                                Hadron_spectra_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num);
                                            }
                                        }

                                        if(eAnalysisNum == 130)
                                        {
                                            Fill_m2_nSigmaP_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num);
                                        }
                                        if(eAnalysisNum == 31)
                                        {
                                            fill_phi_hist_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,14,1,track_counter,event_num);
                                        }
                                        if(eAnalysisNum == 30)
                                        {
                                            pp_correlation_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,14,0,track_counter,event_num);
                                            pp_correlation_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,15,15,0,track_counter,event_num);
                                            pp_correlation_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,15,0,track_counter,event_num);
                                        }
                                        if(eAnalysisNum == 9)
                                        {
                                            comb_counter_global = 0;
                                            Vertex_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,14,1,track_counter,event_num,x_mean,y_mean,z_mean,1);
                                            Event_Display(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num);
                                        }
                                        if(eAnalysisNum == 23)
                                        {
                                            //K_Pi_three_analysis(PID_counter_Array,PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,PID_Array,
                                            //                    PicoAlexEvent_use,PicoAlexEvent_use,PicoAlexEvent_use,8,8,9,9,2,track_counter,track_counter,event_num,eSE_ME_Flag,eAnalysisNum);
                                        }
                                        if(eAnalysisNum == 20)
                                        {
                                            Merged_Omega_Xi_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                                     PicoAlexEvent_use,PicoAlexEvent_use,14,9,2,eSE_ME_Flag);
                                            Merged_Omega_Xi_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                                     PicoAlexEvent_use,PicoAlexEvent_use,15,8,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 351)
                                        {
                                            //LambdaCPlus_V0_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                            //                        PicoAlexEvent_use,PicoAlexEvent_use,14,9,8,2,eSE_ME_Flag);
                                            LambdaCPlus_V0_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                                    PicoAlexEvent_use,PicoAlexEvent_use,8,9,14,2,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 200)
                                        {
                                            //cout << "event_num = " << event_num << endl;
                                            //if(erefMult_bin >= 0 && erefMult_bin < 9 && erefMult_bin == 0)
                                            //{
                                            //    cout << "A dummy_counter = " << dummy_counter_array[erefMult_bin] << ", hist entries = "
                                            //        << h1D_Omega2250_m_cent[erefMult_bin]->GetEntries() << endl;
                                            //}
                                            Omega2250_analysis(PID_counter_Array,PID_counter_Array,
                                                               PID_Array,PID_Array,
                                                               PicoAlexEvent_use,PicoAlexEvent_use,
                                                               14,9,9,8,12,2,eSE_ME_Flag);
                                            //if(erefMult_bin >= 0 && erefMult_bin < 9 && erefMult_bin == 0)
                                            //{
                                            //    cout << "B dummy_counter = " << dummy_counter_array[erefMult_bin] << ", hist entries = "
                                            //        << h1D_Omega2250_m_cent[erefMult_bin]->GetEntries() << endl;
                                            //}
                                        }
                                        if(eAnalysisNum == 201) // same event
                                        {
#if 1
                                            // standard analysis
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array,
                                                         PID_Array,PID_Array,
                                                         PicoAlexEvent_use,PicoAlexEvent_use,
                                                             14,9,9,14,9,2,eSE_ME_Flag,0,0);
#endif

#if 0
                                            // "Omega" instead of Xi
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array,
                                                         PID_Array,PID_Array,
                                                         PicoAlexEvent_use,PicoAlexEvent_use,
                                                             14,9,12,14,9,2,eSE_ME_Flag,0);
#endif

#if 0
                                            // anti-Lambda + Xi-
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array,
                                                         PID_Array,PID_Array,
                                                         PicoAlexEvent_use,PicoAlexEvent_use,
                                                             14,9,9,15,8,2,eSE_ME_Flag,0);
#endif

#if 0
                                            // K0S + K0S
                                            d4s_analysis_new(PID_counter_Array,PID_counter_Array,
                                                         PID_Array,PID_Array,
                                                         PicoAlexEvent_use,PicoAlexEvent_use,
                                                             8,9,9,8,9,2,eSE_ME_Flag,0);
#endif

#if 0
                                            d4s_analysis(PID_counter_Array,PID_counter_Array,
                                                         PID_Array,PID_Array,
                                                         PicoAlexEvent_use,PicoAlexEvent_use,
                                                         14,9,9,15,8,2,eSE_ME_Flag,0);
                                            //d4s_analysis(PID_counter_Array,PID_counter_Array,
                                            //             PID_Array,PID_Array,
                                            //             PicoAlexEvent_use,PicoAlexEvent_use,
                                            //             14,9,9,14,9,2,eSE_ME_Flag,1);
#endif
                                        }
                                        if(eAnalysisNum == 300)
                                        {
                                            //Jet_Analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,0);
                                        }
                                        if(eAnalysisNum == 301)
                                        {
                                            //Jet_Analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,1);
                                        }
                                        if(eAnalysisNum == 302)
                                        {
                                            FillJet_Analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 202) // same event
                                        {
                                            Lambda_mixed_analysis(PID_counter_Array,PID_counter_Array,
                                                         PID_Array,PID_Array,
                                                         PicoAlexEvent_use,PicoAlexEvent_use,
                                                             14,9,9,14,9,2,eSE_ME_Flag,0);
                                        }
                                        if(eAnalysisNum == 15)
                                        {
                                            rho_analysis(PID_counter_Array,PID_counter_Array,PID_Array,PID_Array,
                                                             PicoAlexEvent_use,PicoAlexEvent_use,8,9,2,track_counter,track_counter,event_num,eSE_ME_Flag);
                                        }
                                        if(eAnalysisNum == 112)
                                        {
                                            BBC_analysis(PicoAlexEvent_use);
                                        }
                                        if(eAnalysisNum == 111)
                                        {
                                            EventPlane_analysis_V2(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num);
                                        }
                                        if(eAnalysisNum == 125)
                                        {
                                            Event_Analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,1,track_counter,event_num);
                                        }

                                        if(eAnalysisNum == 26)
                                        {
                                            EventPlane_Patrick_analysis(PID_counter_Array,PID_Array,PicoAlexEvent_use,14,eBeamTimeNum,1,event_num);
                                            MyLeptonTreeAna(PID_counter_Array,PID_Array,PicoAlexEvent_use,3);
                                        }
                                        if(eAnalysisNum == 124 && eSE_ME_Flag == 0 )
                                        {
                                            event_A_ana      = PicoAlexEvent_use;

                                            StThreeVectorF vectorprim_event;
                                            vectorprim_event.set(event_A_ana->primaryVertex().x(),event_A_ana->primaryVertex().y(),event_A_ana->primaryVertex().z());
                                            if( vectorprim_event.perp2() < Event_radius_cut*Event_radius_cut && fabs(vectorprim_event.z()) < vertex_z_cut && event_A_ana   ->isMinBias() )
                                            {
                                                my_event_single->clearTrackLists();
                                                my_event_single->SetHeader((Float_t)event_A_ana->refMult(),(Float_t)event_A_ana->runId());
                                                my_event_single->SetPrim_X((Float_t)event_A_ana->primaryVertex().x());
                                                my_event_single->SetPrim_Y((Float_t)event_A_ana->primaryVertex().y());
                                                my_event_single->SetPrim_Z((Float_t)event_A_ana->primaryVertex().z());

                                                SingleTrackTree(PID_counter_Array,PID_Array,PicoAlexEvent_use,my_event_single,3,2,0); // e- e+
#if 0 // don't produce like sign histos for saving calculation time ( vertex determination ! )
                                                SingleTrackTree(PID_counter_Array,PID_Array,PicoAlexEvent_use,my_event_single,3,3,0); // e- e-
                                                SingleTrackTree(PID_counter_Array,PID_Array,PicoAlexEvent_use,my_event_single,2,2,0); // e+ e+
#endif
                                                PiKP_SingleTrackTree(PID_counter_Array,PID_Array,PicoAlexEvent_use,my_event_single,0); // pi K p

                                                my_tree_single->Fill();
                                                my_event_single->Clear();
                                            }
                                        }
                                    }
                                    //**********************************************************************************************************************************
                                }
                            }

                            //************************************************************************
                            // delete the event_ME_array objects
                            if(
                               Start_mixing == 1
                              )
                            {
                                PicoAlexEvent_use = NULL;
                                PicoAlexEventB    = NULL;
                                event_ME_array->Clear();
                            }
                            //************************************************************************

                        }
                    }
                }
            }
        }
    } // End event loop

    return 1;
}

Int_t Analysis_Snurf::Finish()
{
    cout << "Finish started" << endl;
    cout << "Number of good events = " << good_event << endl;
    cout << "Number of bad events = " << bad_event << endl;
    cout << "Number of total events = " << total_event << endl;
    cout << "Number of events where analysis was started = " << dummy_counter << endl;
    cout << "Number of invariant mass entries = " << dummy_counter_A << endl;


    //****************************** Event display *********************************************************
    cout << "" << endl;
    cout << "----------------------- Event display ---------------------" << endl;
    cout << "vertex = {" << x_mean << ", " << y_mean << ", " << z_mean << "}" << endl;
    cout << "-----------------------------------------------------------" << endl;
    cout << "" << endl;

    pPrimaryVertex->SetPosition(x_mean,y_mean,z_mean);
    cEventDisplay->cd();

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_endcaps[r] ->Draw();
    }

    for(Int_t i = 0; i < n_poly_track_counter; i++)
    {
        pTrack[i]             -> Draw();
        pTrack_extrapolate[i] -> Draw();
    }

    BeamLine       ->Draw();
    XAxis          ->Draw();
    YAxis          ->Draw();
    ZAxis          ->Draw();
    pPrimaryVertex ->Draw();

    cEventDisplay->Update();

    //******************************************************************************************************


    //******************************************************************************************************

    cVertex_distr->Divide(1,3);
    cVertex_distr->cd(1)->SetTicks(1,1);
    cVertex_distr->cd(1);
    hx_vertex_distr->DrawCopy("h");
    cVertex_distr->cd(2);
    hy_vertex_distr->DrawCopy("h");
    cVertex_distr->cd(3);
    hz_vertex_distr->GetXaxis()->SetRangeUser(-200,200);
    hz_vertex_distr->DrawCopy("h");

    //******************************************************************************************************

    cout << "Save data to file" << endl;
    if(eAnalysisNum == 1 || eAnalysisNum == 12 || eAnalysisNum == 25)
    {
        Tree_PhiMeson_v2->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 2 || eAnalysisNum == 14 || eAnalysisNum == 3)
    {
        Lambda_X_NT ->Write();
    }
    if(eAnalysisNum == 400)
    {
        cout << "Write D0 tree" << endl;

        h_dca_diff           ->Write();
        h_decay_length_diff  ->Write();

        Tree_D0_v2->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 24)
    {
        ThetaPlus_NT ->Write();
    }
    if(eAnalysisNum == 20)
    {
        Tree_OmegaPV0_v2   ->Write("",TObject::kOverwrite);
        Tree_OmegaMV0_v2   ->Write("",TObject::kOverwrite);
        Tree_XiPV0_v2      ->Write("",TObject::kOverwrite);
        Tree_XiMV0_v2      ->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 102)
    {
        Tree_Sigma0_v2   ->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 26)
    {
        my_tree->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 8)
    {
        Tree_hadron_v2->Write();
        Tree_hadron_v2->SetAutoSave( 5000000 );
    }
    if(eAnalysisNum == 112)
    {
        BBC_tree->Write();
    }
    if(eAnalysisNum == 11)
    {
        EventPlane_NT   ->AutoSave("SaveSelf");
    }
    if(eAnalysisNum == 125)
    {
        EventAnalysis_NT   ->AutoSave("SaveSelf");
        for(Int_t k = 0; k < 2; k++)
        {
            p_mean_pt_vs_charge_asymm[k] ->Write();
        }
    }
    if(eAnalysisNum == 111)
    {
        heta_EP ->Write();
        hFTPC_BBC_corr->Write();
        h_FTPC_phi->Write();
        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics (5)
        {
            EventPlane_array_NT[n_harm]   ->AutoSave("SaveSelf");
        }
        for(Int_t n_ew = 0; n_ew < 2; n_ew++)
        {
            hrefMult_BBC_hits_East_West[n_ew] ->Write();
            for(Int_t n_tiles = 0; n_tiles < 24; n_tiles++)
            {
                hBBC_ADC_tiles_Eeast_West[n_ew][n_tiles] ->Write();
                pBBC_ADC_tiles_vs_RefMult_Eeast_West[n_ew][n_tiles] ->Write();
                hBBC_ADC_tiles_vs_RefMult_Eeast_West[n_ew][n_tiles] ->Write();
            }
        }
    }
    if(eAnalysisNum == 13)
    {
        PhiKPKM_V0_NT      ->AutoSave("SaveSelf");
    }
    if(eAnalysisNum == 17 || eAnalysisNum == 18)
    {
        Xi1530_NT       ->AutoSave("SaveSelf");
    }
    if(eAnalysisNum == 22)
    {
        RefMult_QA_NT      ->AutoSave("SaveSelf");
    }
    if(eAnalysisNum == 23)
    {
        Tree_KaonV0_v2 ->Write("",TObject::kOverwrite);
        //Kaon_NT     ->AutoSave("SaveSelf");
    }
    if(eAnalysisNum == 9)
    {
        Vertex_NT     ->AutoSave("SaveSelf");
        EDisplay_NT   ->AutoSave("SaveSelf");
        cEventDisplay    ->Write();
        cMinDCA          ->Write();
        for(Int_t i = 0; i < ED_counter; i++)
        {
            cEventDisplay_array[i] ->Write();
            cEventDCA_array[i] ->Write();
        }
    }
    Outputfile       ->cd();
    heta             ->Write();
    hphi_eta_pos     ->Write();
    hphi_eta_neg     ->Write();
    hphi             ->Write();
    hsign            ->Write();
    hTPC_dEdx_vs_p   ->Write();
    hTOF_vs_p        ->Write();
    hMass            ->Write();
    hnon_prim_to_prim_ratio     ->Write();
    htracks_left_to_right_ratio ->Write();
    htracks_left_to_right_diff  ->Write();
    hTheta_vs_phi               ->Write();
    h_nhitsFit->Write();
    h_nhitsFit_TOF->Write();
    for(Int_t xa = 0; xa < 10; xa++)
    {
        h_eta_phi[xa] ->Write();
    }

    if(eAnalysisNum == 22)
    {
        Outputfile      ->mkdir("phicorr");
        Outputfile      ->cd("phicorr");
        for(Int_t i = 0; i < nPhi_corr_days; i++)
        {
            if(hPhi_days_use[i] == 1)
            {
                for(Int_t j = 0; j < nPhi_corr_z_bins; j++)
                {
                    for(Int_t k = 0; k < nPhi_corr_eta_bins; k++)
                    {
                        for(Int_t l = 0; l < 2; l++)
                        {
                            for(Int_t m = 0; m < nPhi_corr_pt_bins; m++)
                            {
                                hPhi_corr[i][j][k][l][m] ->Write();
                            }
                        }
                    }
                }
            }
        }
    }

    if(eAnalysisNum == 6 || eAnalysisNum == 7)  // D0 + KStar or KStarPM analysis
    {
        Outputfile        ->cd("");
        Outputfile        ->mkdir("KStar");
        Outputfile        ->cd("KStar");

        for(Int_t j = 0; j < npt_bins_D0; j++)
        {
            h_inv_mass_D0_pt[j] ->Write();
        }

    }
    if(eAnalysisNum == 130)
    {
        Outputfile->cd();
        for(Int_t nrefmult = 0; nrefmult < N_refmult_bins; nrefmult++)
        {
            for(Int_t npt = 0; npt < N_2D_m2_nSigmaP_pT_bins; npt++)
            {
                for(Int_t ncharge = 0; ncharge < 2; ncharge++)
                {
                    h2D_m2_nSigmaP[nrefmult][npt][ncharge] ->Write();
                }
            }
        }
        Outputfile->Close();
    }

    if(eAnalysisNum == 21) // Sigma1385 analysis
    {
        Outputfile        ->cd("");
        Outputfile        ->mkdir("Sigma1385");
        Outputfile        ->cd("Sigma1385");

        for(Int_t j = 0; j < npt_bins_Sigma1385; j++)
        {
            h_inv_mass_Sigma1385_pt[j] ->Write();
        }

    }

    if ( eAnalysisNum == 222 || eAnalysisNum == 124 )
    {
        Outputfile      ->cd("");
        hDecVertex_YX->Write();
        hDecVertex_YZ->Write();
        for ( Int_t i = 0; i < 3; ++i ) {
            hInvMassAB[i]->Write();
            hDcaDiff[i]->Write();
        }
        for ( Int_t i = 0; i < 4; ++i ) hnSigElDist[i]->Write();
        my_tree_single->Write();
        my_tree_single->SetAutoSave(5000000);
        Outputfile->Close();
        delete my_event_single; my_event_single = 0;
        cout << "Patrick's tree written" << endl;
    }

    if( eAnalysisNum == 31)
    {
        for(Int_t n_mult = 0; n_mult < N_refmult_bins; n_mult++)
        {
            for(Int_t n_w = 0; n_w < 2; n_w++)
            {
                h_phi_vs_track_number[n_mult][n_w] ->Write();
                h_pT_vs_track_number[n_mult][n_w]  ->Write();
                h_mean_pt_vs_phi[n_mult][n_w]      ->Write();
            }
        }
    }

    if(eAnalysisNum == 88)
    {
        Outputfile->cd();
        for(Int_t i = 0; i < 9; i++)
        {
            for(Int_t j = 0; j < npt_bins_pt_spectra; j++)
            {
                for(Int_t k = 0; k < n_hadrons_pt_spectra; k++)
                {
                    h_NewY_vs_NewX_mult_pt[i][j][k] ->Write();

                    for(Int_t l = 0; l < 2; l++)
                    {
                        h_m2_mult_pt[i][j][k][l] ->Write();
                    }
                }
            }
        }
        for(Int_t i = 0; i < 2; i++)
        {
            h_ToF_local_YZ_vs_m2[i] ->Write();
        }
        Outputfile->Close();
    }

    if( eAnalysisNum == 30)
    {
        for(Int_t i = 0; i < 2; i++)
        {
            for(Int_t j = 0; j < 3; j++)
            {
                h_phi_pp_corr[i][j]     ->Write();
                h_cos_phi_pp_corr[i][j] ->Write();
            }
        }

        for(Int_t n_mult = 0; n_mult < N_refmult_bins; n_mult++)
        {
            for(Int_t n_w = 0; n_w < 2; n_w++)
            {
                for(Int_t n_harm = 0; n_harm < n_pp_harm; n_harm++)
                {
                    for(Int_t n_pp = 0; n_pp < 3; n_pp++)
                    {
                        p_pp_corr[n_mult][n_w][n_harm][n_pp] ->Write();
                    }
                }
            }
        }
    }

    if(eAnalysisNum == 15)
    {
        for(Int_t j = 0; j < N_pipi_spectra; j++)
        {
            for(Int_t k = 0; k < N_pipi_pt_bins; k++)
            {
                h_inv_mass_pipi[j][k]->Write();
            }
        }
    }

    if(eAnalysisNum == 200 || eAnalysisNum == 201 || eAnalysisNum == 203)
    {
        //for(Int_t m = 0; m < 9; m++)
        //{
        //    h2D_Omega2250_m_vs_p_cent[m]  ->Write();
        //    h1D_Omega2250_m_cent[m]       ->Write();
        //}
        //h1D_Lambda   ->Write();
        //h1D_Lambda_B ->Write();
        //h1D_Xi       ->Write();
        //d4s_NT       ->Write();

        Tree_d4s->Write("",TObject::kOverwrite);
        for(Int_t i = 0; i < 4; i++)
        {
            h_d4s_InvMassAB[i] ->Write();
        }
    }
    if(eAnalysisNum == 202)
    {
        Tree_Lambda->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 35)
    {
        for(Int_t i = 0; i < N_InvMass_lambdaCplus; i++)
        {
            h_2D_InvMass_LambdaCplus[i] ->Write();
        }
    }
    if(eAnalysisNum == 351)
    {
        for(Int_t i = 0; i < N_InvMass_lambdaCplus_V0_Lambda; i++)
        {
            h_InvMass_AB_LambdaCplus[i]   ->Write();
        }
        for(Int_t i = 0; i < N_InvMass_lambdaCplus_V0; i++)
        {
            h_InvMass_ABC_LambdaCplus[i]  ->Write();
        }

        h_massA_LambdaCplus      ->Write();
        h_massB_LambdaCplus      ->Write();
        h_massC_LambdaCplus      ->Write();
        h_massA_corr_LambdaCplus ->Write();
        h_massB_corr_LambdaCplus ->Write();

    }
    if(eAnalysisNum == 302)
    {
        Tree_JetTrackEvent->Write("",TObject::kOverwrite);
    }
    if(eAnalysisNum == 300 || eAnalysisNum == 301)
    {
        for(Int_t i = 0; i < N_jet_histos; i++)
        {
            h_jet_pt[i]        ->Write();
            h_jet_pt_sub[i]    ->Write();
            h_jet_area[i]      ->Write();
            h_jet_per_event[i] ->Write();
        }
    }

    if(eAnalysisNum == 360)
    {
        for(Int_t i_charge = 0; i_charge < 2; i_charge++)
        {
            p_pt_Ach[i_charge] ->Write();
        }

        h_Ach ->Write();
    }

    cout << "Saved data to file" << endl;

    return 1;
}
