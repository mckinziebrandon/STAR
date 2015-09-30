#include "EP_Det_Sim.h"

static const Double_t Pi = TMath::Pi();
static TF1 *f_RefMult;
static TF1 *f_Spectrum_pt;
static TF1 *f_v2_centrality;
static TF1 *fGauss_v2_smear;
static TF1 *fFlowFit_track;
static TH1D *hPhi_Psi;
static TH1D *hPhi_Psi_pos_eta;
static TH1D *hPhi_Psi_neg_eta;
static TH1D *hPhi;

static const Int_t N_centralities = 10; // 16
//static const Double_t cent_start_values[N_centralities] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
//static const Double_t cent_stop_values[N_centralities]  = {5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};
static const Double_t cent_start_values[N_centralities] = {0,3,6,10,15,20,25,30,35,40};
static const Double_t cent_stop_values[N_centralities]  = {3,6,10,15,20,25,30,35,40,45};
static Double_t cent_mean_values[N_centralities];
static TF1*  f_v2_input[N_centralities];

static Double_t total_RefMult;
static Double_t refmult_start_values[N_centralities];
static Double_t refmult_stop_values[N_centralities];

static TRandom ran;
static TString HistName;

static TFile* Outputfile;

static const char HA_EVENT_TREE[]   = "HA_Events";
static TTree         *Tree_EP_Det_Sim;
static EP_Det_Sim_Event  HA_event;
static EP_Det_Sim_Event *HA_event_ptr;
static EP_Det_Sim_Track *HA_track;
static TRandom3 r3;

static const Int_t N_dNdeta_fit_par = 4; // a, alpha, beta, c
static const Int_t N_dNdetaMod_fit_par = 6; // a, alpha, beta, c, par5, par6
static const Int_t N_dNdeta_mult    = 10; // 0-3%, 3-6%, 6-10%, 10-15%, 15-20%, 20-25%, 25-30%, 30-35%, 35-40%, 40-45%
//static TF1* dNdeta;
static TH1D* hdNdeta_fit_params[N_dNdeta_fit_par];
static TH1D* hdNdetaMod_fit_params[N_dNdetaMod_fit_par];
static TGraphAsymmErrors *Phobos_scan_dNdeta_mult[N_dNdeta_mult];
static TF1* fdNdeta[N_dNdeta_mult];
static TF1* fdNdetaMod[N_dNdeta_mult];
static Double_t NTracks_dNdeta[N_dNdeta_mult];
static Double_t NTracks_dNdeta_full[N_dNdeta_mult];
static Double_t NTracks_refMult[N_dNdeta_mult];
static Double_t NScaleFac_refMult[N_dNdeta_mult];
static Double_t aver_scale_fac_dNdeta = 0.0;


#include "EP_Det_Sim_Func.h"



ClassImp(EP_Det_Sim)
    //_____________________________________________________________________________
    EP_Det_Sim::EP_Det_Sim() {

    }


//_____________________________________________________________________________
EP_Det_Sim::~EP_Det_Sim() {

}



Int_t EP_Det_Sim::Init()
{
    cout << "***************** Initializing objects *******************" << endl;

    cout << "Create output file" << endl;

    r3.SetSeed(0);

    // open appropriate dN/deta root file
    if(eBeamtime == 0)
    {
        dNdeta_file = TFile::Open("./Data/dNdeta_19GeV.root");  // open the file
    }
    if(eBeamtime == 1)
    {
        dNdeta_file = TFile::Open("./Data/dNdeta_200GeV.root");  // open the file
    }

    // obtain hist arrays of fit parameters (per centrality) for dN/deta
    for(Int_t fit_par = 0; fit_par < N_dNdeta_fit_par; fit_par++)
    {
        HistName = "hdNdeta_fit_params";
        HistName += fit_par;
        hdNdeta_fit_params[fit_par] = (TH1D*)dNdeta_file->Get(HistName.Data());
    }
    for(Int_t fit_par = 0; fit_par < N_dNdetaMod_fit_par; fit_par++)
    {
        HistName = "hdNdetaMod_fit_params";
        HistName += fit_par;
        hdNdetaMod_fit_params[fit_par] = (TH1D*)dNdeta_file->Get(HistName.Data());
    }

    // obtain dN/deta fitted function & integrate. Yields NTracks_dNdeta for all nMult values. 
    Int_t nMult_counter = 0;
    for(Int_t nMult = 0; nMult < N_dNdeta_mult; nMult++)
    {
        HistName = "Phobos_scan_dNdeta_mult";
        HistName += nMult;
        Phobos_scan_dNdeta_mult[nMult] = (TGraphAsymmErrors*)dNdeta_file->Get(HistName.Data());

        HistName = "fdNdeta_mult";
        HistName += nMult;
        fdNdeta[nMult] = new TF1(HistName.Data(),dNdetaFunc,eta_low,eta_high,4);
        for(Int_t fit_par = 0; fit_par < N_dNdeta_fit_par; fit_par++)
        {
            fdNdeta[nMult] ->SetParameter(fit_par,hdNdeta_fit_params[fit_par]->GetBinContent(hdNdeta_fit_params[fit_par]->FindBin(nMult)));
        }
        //NTracks_dNdeta[nMult]      = fdNdeta[nMult]->Integral(-0.5,0.5);
        //NTracks_dNdeta_full[nMult] = fdNdeta[nMult]->Integral(-4.5,4.5);

        HistName = "fdNdetaMod_mult";
        HistName += nMult;
        fdNdetaMod[nMult] = new TF1(HistName.Data(),dNdetaModFunc,eta_low,eta_high,6);
        for(Int_t fit_par = 0; fit_par < N_dNdetaMod_fit_par; fit_par++)
        {
            fdNdetaMod[nMult] ->SetParameter(fit_par,hdNdetaMod_fit_params[fit_par]->GetBinContent(hdNdetaMod_fit_params[fit_par]->FindBin(nMult)));
        }
        NTracks_dNdeta[nMult]      = fdNdetaMod[nMult]->Integral(-0.5,0.5);
        NTracks_dNdeta_full[nMult] = fdNdetaMod[nMult]->Integral(eta_low,eta_high);

        if(NTracks_dNdeta[nMult] > 0.0)
        {
            aver_scale_fac_dNdeta += (NTracks_dNdeta_full[nMult]/NTracks_dNdeta[nMult]);
            nMult_counter++;
            cout << "nMult = " << nMult << ", Int[-0.5,0.5]: " << NTracks_dNdeta[nMult] << ", Int[" << eta_low << "," << eta_high << "]: " << NTracks_dNdeta_full[nMult]
                << ", dNdeta full/1 = " << NTracks_dNdeta_full[nMult]/NTracks_dNdeta[nMult] << endl;
        }
    }
    if(nMult_counter > 0) aver_scale_fac_dNdeta /= ((Double_t)nMult_counter);
    cout << "aver_scale_fac_dNdeta = " << aver_scale_fac_dNdeta << endl;

    Outputfile        = new TFile(output_file_name.Data(),"RECREATE");

    //-------------------------------------------------------
    cout << "Create c_RefMult" << endl;
    // Result for a fit to the 19.6 GeV STAR reference multiplicity
    //p0                        =	2201.83      	+/-	22.977
    //p1                        =	-32.9937     	+/-	0.840205
    //p2                        =	0.238481     	+/-	0.0105491
    //p3                        =	-0.000828501 	+/-	5.87267e-05
    //p4                        =	1.29415e-06  	+/-	1.49151e-07
    //p5                        =	-7.12945e-10 	+/-	1.41132e-10

    // Result for a fit to the 200 GeV STAR reference multiplicity
    //p0           1.66357e+04   1.80725e+01   5.01354e-03   1.44206e-06
    //p1          -2.51869e+02   6.98241e-02   6.00504e-05   7.29100e-04
    //p2           1.90088e+00   1.33309e-04   4.53204e-07   4.30237e-01
    //p3          -6.76806e-03   2.42040e-07   1.61363e-09   2.34198e+02
    //p4           8.96086e-06   4.24671e-10   2.13643e-12   1.44647e+05
    //p5           7.34452e-09   7.23539e-13   1.75107e-15   8.15232e+07
    //p6          -3.03284e-11   1.19461e-15   7.23085e-18   4.94759e+10
    //p7           2.12450e-14   1.89589e-18   5.06519e-21   2.64153e+13



    Double_t start_RefMult = 20.0;
    Double_t stop_RefMult  = 360.0;
    f_RefMult               = new TF1("f_RefMult",PolyFitFunc,0.0,10,8);
    for(Int_t f = 0; f < 8; f++)
    {
        f_RefMult->ReleaseParameter(f);
        f_RefMult->SetParError(f,0.0);
        f_RefMult->SetParameter(f,0.0);
    }
    if(eBeamtime == 0) // 19.6 GeV
    {
        f_RefMult->SetParameter(0,2201.83);
        f_RefMult->SetParameter(1,-32.9937);
        f_RefMult->SetParameter(2,0.238481);
        f_RefMult->SetParameter(3,-0.000828501);
        f_RefMult->SetParameter(4,1.29415e-06);
        f_RefMult->SetParameter(5,-7.12945e-10);
        f_RefMult->SetParameter(6,0.0);
        f_RefMult->SetParameter(7,0.0);
        start_RefMult = 20.0;
        stop_RefMult  = 360.0;
    }
    if(eBeamtime == 1) // 200 GeV
    {
        f_RefMult->SetParameter(0,1.66357e+04);
        f_RefMult->SetParameter(1,-2.51869e+02);
        f_RefMult->SetParameter(2,1.90088e+00);
        f_RefMult->SetParameter(3,-6.76806e-03);
        f_RefMult->SetParameter(4,8.96086e-06);
        f_RefMult->SetParameter(5,7.34452e-09);
        f_RefMult->SetParameter(6,-3.03284e-11);
        f_RefMult->SetParameter(7,2.12450e-14);
        start_RefMult = 50.0;
        stop_RefMult  = 600.0;
    }
    f_RefMult->SetRange(start_RefMult,stop_RefMult);

    total_RefMult = f_RefMult->Integral(start_RefMult,stop_RefMult); // 60-65%

    Int_t centrality_bin = 0;
    // Define centralities
    Double_t use_stop_RefMult = stop_RefMult;
    for(Int_t n_refmult = ((Int_t)stop_RefMult); n_refmult >= ((Int_t)start_RefMult); n_refmult--) // loop over reference multiplicity bins
    {
        if(centrality_bin >= N_centralities) break;
        //cent_start_values
        Double_t Int_RefMult  = f_RefMult->Integral(n_refmult,use_stop_RefMult);
        Double_t frac_RefMult = Int_RefMult/(total_RefMult*(1.0/0.62))*100.0; // 100%
        Double_t frac_compare = cent_stop_values[centrality_bin]-cent_start_values[centrality_bin];

        if(frac_RefMult >= frac_compare)
        {
            cout << "frac_RefMult = " << 100.0*frac_RefMult << ", frac_compare  = " << frac_compare << ", n_refmult = " << n_refmult << endl;
            refmult_start_values[centrality_bin] = n_refmult;
            refmult_stop_values[centrality_bin]  = use_stop_RefMult;
            use_stop_RefMult = n_refmult;
            centrality_bin++;
        }
    }

    f_RefMult->SetRange(refmult_start_values[N_centralities-1],stop_RefMult);

    for(Int_t nMult = 0; nMult < N_dNdeta_mult; nMult++)
    {
        Double_t val_ref = 0.0;
        Double_t wei_ref = 0.0;
        for(Int_t nref = refmult_start_values[nMult]; nref < refmult_stop_values[nMult]; nref++)
        {
            val_ref =  f_RefMult->Eval(nref)*nref;
            wei_ref =  f_RefMult->Eval(nref);
        }
        if(wei_ref > 0.0) val_ref /= wei_ref;
        NTracks_refMult[nMult] = val_ref;

        //NTracks_refMult[nMult] = f_RefMult->Integral(refmult_start_values[nMult],refmult_stop_values[nMult]);
        if(NTracks_refMult[nMult] > 0.0)
        {
            NScaleFac_refMult[nMult] = NTracks_dNdeta[nMult]/NTracks_refMult[nMult];
        }
    }

    for(Int_t nMult = 0; nMult < N_dNdeta_mult; nMult++)
    {
        cout << "Centrality bin = " << nMult << ", IntdNdeta(PHOBOS) = " << NTracks_dNdeta[nMult] << ", NrefMult(STAR) = "
            << NTracks_refMult[nMult] << ", scale_fac = " << NScaleFac_refMult[nMult] << endl;
    }
    //-------------------------------------------------------




    //-------------------------------------------------------
    fGauss_v2_smear              = new TF1("fGauss_v2_smear",GaussFitFunc,0.0,10,3);
    fFlowFit_track               = new TF1("fFlowFit_track",FlowFitFunc,-Pi,Pi,5);
    fGauss_v2_smear->SetParameter(0,1.0);     // amplitude
    fGauss_v2_smear->SetParameter(1,1.0);     // mean
    fGauss_v2_smear->SetParameter(2,0.05); // sigma
    //-------------------------------------------------------



    //-------------------------------------------------------
    //dNdeta                = new TF1("dNdeta",dNdetaFunc,-5.0,5.0,4);
    //dNdeta ->SetParameter(0,0.9); // a ~width
    //dNdeta ->SetParameter(1,1.25); // alpha ~defines dip at eta = 0
    //dNdeta ->SetParameter(2,1.0); // beta ~width?
    //dNdeta ->SetParameter(3,230); // c  ~amplitude
    //dNdeta ->SetLineColor(1);
    //dNdeta ->SetLineWidth(1);
    //dNdeta ->SetLineStyle(1);
    //-------------------------------------------------------



    //-------------------------------------------------------
    hPhi_Psi = new TH1D("hPhi_Psi","hPhi_Psi",400,-Pi,Pi);
    hPhi_Psi_pos_eta = new TH1D("hPhi_Psi_pos_eta","hPhi_Psi_pos_eta",400,-Pi,Pi);
    hPhi_Psi_neg_eta = new TH1D("hPhi_Psi_neg_eta","hPhi_Psi_neg_eta",400,-Pi,Pi);
    hPhi     = new TH1D("hPhi","hPhi",400,-3.0*Pi,3.0*Pi);
    //-------------------------------------------------------



    //-------------------------------------------------------
    // Polynomial fit to v2(centrality) 19.6 GeV (central -> peripheral in %)
    //p0                        =	0.00994233   	+/-	0.000987568
    //p1                        =	0.00282346   	+/-	0.000278221
    //p2                        =	-5.81173e-05 	+/-	2.25525e-05
    //p3                        =	7.15709e-07  	+/-	7.42984e-07
    //p4                        =	-7.72361e-09 	+/-	1.0601e-08
    //p5                        =	4.20936e-11  	+/-	5.44115e-11
    f_v2_centrality               = new TF1("f_v2_centrality",PolyFitFunc,0.0,10,7);
    for(Int_t f = 0; f < 7; f++)
    {
        f_v2_centrality->ReleaseParameter(f);
        f_v2_centrality->SetParError(f,0.0);
        f_v2_centrality->SetParameter(f,0.0);
    }
    f_v2_centrality->SetParameter(0,0.00994233);
    f_v2_centrality->SetParameter(1,0.00282346);
    f_v2_centrality->SetParameter(2,-5.81173e-05);
    f_v2_centrality->SetParameter(3,7.15709e-07);
    f_v2_centrality->SetParameter(4,-7.72361e-09);
    f_v2_centrality->SetParameter(5,4.20936e-11);
    f_v2_centrality->SetParameter(6,0.0);
    Double_t start_cent = 0.0;
    Double_t stop_cent  = 45.0;
    f_v2_centrality->SetRange(start_cent,stop_cent);
    //-------------------------------------------------------



    //-------------------------------------------------------
    for(Int_t bin_cent = 0; bin_cent < N_centralities; bin_cent++)
    {
        Double_t bin_cent_val = (cent_start_values[bin_cent]+cent_stop_values[bin_cent])/2.0;
        Double_t v2_cent_val  = f_v2_centrality->Eval(bin_cent_val);
        Double_t v2_scaling_factor = 0.7*v2_cent_val/f_v2_centrality->Eval(5.0);

        HistName = "f_v2_input_";
        HistName += bin_cent;
        //f_v2_input[bin_cent] = (TF1*)v2_pT_Fit->Clone(HistName.Data());
        f_v2_input[bin_cent] = new TF1(HistName.Data(),v2_pT_FitFunc,0.0,3.0,6);
        f_v2_input[bin_cent]->SetParameter(0,3.0); // number-of-constituent quarks
        //f_v2_input[bin_cent]->SetParameter(1,6.78950e-02*0.5*(1.0+((Double_t)bin_cent)*0.1)); // a
        f_v2_input[bin_cent]->SetParameter(1,6.78950e-02*0.5*1.97144); // a
        f_v2_input[bin_cent]->SetParameter(2,3.82098e-01); // b
        f_v2_input[bin_cent]->SetParameter(3,1.59207e-01); // c
        f_v2_input[bin_cent]->SetParameter(4,9.94995e-03); // d
        f_v2_input[bin_cent]->SetParameter(5,v2_scaling_factor*0.8); // f
        cout << "centrality bin = " << bin_cent << ", v2_scaling_factor = " << v2_scaling_factor << endl;
    }
    //-------------------------------------------------------



    //-------------------------------------------------------
    // Polynomial fit parameters for 0-3 GeV/c to 19.6 GeV uncorrected proton dN/dpT spectrum
    //p0                        =	-5.7435e+08  	+/-	196875
    //p1                        =	2.97222e+09  	+/-	1.18041e+06
    //p2                        =	-3.94328e+09 	+/-	2.39286e+06
    //p3                        =	2.33314e+09  	+/-	2.263e+06
    //p4                        =	-6.90021e+08 	+/-	1.09403e+06
    //p5                        =	9.70298e+07  	+/-	262292
    //p6                        =	-4.84697e+06 	+/-	24733.1
    f_Spectrum_pt               = new TF1("f_Spectrum_pt",PolyFitFunc,0.0,10,7);
    for(Int_t f = 0; f < 7; f++)
    {
        f_Spectrum_pt->ReleaseParameter(f);
        f_Spectrum_pt->SetParError(f,0.0);
        f_Spectrum_pt->SetParameter(f,0.0);
    }
    f_Spectrum_pt->SetParameter(0,-5.7435e+08);
    f_Spectrum_pt->SetParameter(1,2.97222e+09);
    f_Spectrum_pt->SetParameter(2,-3.94328e+09);
    f_Spectrum_pt->SetParameter(3,2.33314e+09);
    f_Spectrum_pt->SetParameter(4,-6.90021e+08);
    f_Spectrum_pt->SetParameter(5,9.70298e+07);
    f_Spectrum_pt->SetParameter(6,-4.84697e+06);
    Double_t start_pt = 0.29;
    Double_t stop_pt  = 3.0;
    f_Spectrum_pt->SetRange(start_pt,stop_pt);
    //-------------------------------------------------------



    //-------------------------------------------------------
    Outputfile->cd();
    HA_event_ptr = &HA_event;
    Tree_EP_Det_Sim   = new TTree("Tree_EP_Det_Sim_tree" , HA_EVENT_TREE );
    Tree_EP_Det_Sim   ->Branch("Tree_EP_Det_Sim_branch"  , "EP_Det_Sim_Event", &HA_event_ptr );
    //-------------------------------------------------------



    return 1;
}



Int_t EP_Det_Sim::Make()
{
    cout << "Make started" << endl;
    r3.SetSeed(0);
    gRandom->SetSeed(0);
    cout << "Seed = " << r3.GetSeed() << endl;

    //------------------------------------------------------------------------------------------------------------------------------------
    const Long64_t N_Events = nEvents;

    // Event loop
    for(Long64_t n_evt = 0; n_evt < N_Events; n_evt++) // loop over all events
    {
        if (n_evt != 0  &&  n_evt % 1 == 0)
            cout << "." << flush;
        if (n_evt != 0  &&  n_evt % 50 == 0)
        {
            Double_t event_percent = 100.0*n_evt/N_Events;
            cout << " " << n_evt << " (" << event_percent << "%) " << "\n" << "==> Processing data, " << flush;
        }

        Int_t N_Tracks = (Int_t)(f_RefMult->GetRandom()); // sample random number of tracks out of reference multiplicity distribution
        Int_t N_Cent   = -1; //
        for(Int_t centrality_bin = 0; centrality_bin < N_centralities; centrality_bin++) // determine centrality bin for this event
        {
            if(N_Tracks >= refmult_start_values[centrality_bin] && N_Tracks < refmult_stop_values[centrality_bin])
            {
                N_Cent = centrality_bin;
            }
        }

        if(N_Cent == -1) continue; // skip events which are outside centrality definition range

        Double_t v2_smear       = fGauss_v2_smear   ->GetRandom(); // v2 event-by-event Gaussian smearing
        //v2_smear = 1.0;
        //Double_t v2_event_smear = v2_track + v2_smear; // smeared v2 for this event
        fFlowFit_track ->SetParameter(0,1.0/(2.0*Pi));
        fFlowFit_track ->SetParameter(1,0.0); // v1
        fFlowFit_track ->SetParameter(2,0.0); // v2
        fFlowFit_track ->SetParameter(3,0.0); // v3
        fFlowFit_track ->SetParameter(4,0.0); // v4

        Double_t Qxy[2][2]; // event Q-vector: [Qx,Qy][subA,subB]
        for(Int_t ixy = 0; ixy < 2; ixy++)
        {
            for(Int_t isub = 0; isub < 2; isub++)
            {
                Qxy[ixy][isub] = 0.0;
            }
        }

        const Int_t N_Tracks_event = N_Tracks*aver_scale_fac_dNdeta;
        Double_t iQxy[N_Tracks_event][2]; // track Q-vector [Qx,Qy]

        Double_t Psi_angle = Pi*ran.Rndm()*2.0-Pi; // Event plane angle [-Pi,Pi]

        // Fill event information
        HA_event.clearTrackList();
        HA_event.setN_Tracks(N_Tracks_event);
        HA_event.setN_Cent(N_Cent);
        HA_event.setv2_smear(v2_smear);
        HA_event.setPsi_in(Psi_angle);

        // Track loop
        for(Int_t n_track = 0; n_track < N_Tracks_event; n_track++) // track loop
        {
            //Double_t pT_track       = ran.Rndm();  // random number between 0 and 1
            //pT_track *= 3.0;
            Double_t pT_track       = f_Spectrum_pt->GetRandom();
            Double_t eta_track      = fdNdetaMod[N_Cent]->GetRandom();
            Double_t v1_track       = 0.02*eta_track;
            Double_t v2_track       = f_v2_input[N_Cent]->Eval(pT_track);
            Double_t v1_track_smear = v1_track*v2_smear; // event-by-event smearing: e.g. fluctuations
            Double_t v2_track_smear = v2_track*v2_smear; // event-by-event smearing: e.g. fluctuations
            //fFlowFit_track ->SetParameter(1,v1_track_smear); // set v2 parameter for flow function
            fFlowFit_track ->SetParameter(2,v2_track_smear); // set v2 parameter for flow function
            //cout << "Integral = " << fFlowFit_track->Integral(-Pi,Pi) << endl;
            //fFlowFit_track ->SetParameter(2,0.5); // set v2 parameter for flow function
            Double_t Psi_phi_track  = fFlowFit_track->GetRandom(); // Psi-phi for this track [-Pi,Pi]
            Double_t phi_track      = Psi_phi_track + Psi_angle; // phi angle of the track [-2Pi,2Pi]

            // [-2Pi,2Pi]-> [-Pi,Pi]
            if(phi_track > Pi)  phi_track -= 2.0*Pi;
            if(phi_track < -Pi) phi_track += 2.0*Pi;

            hPhi_Psi->Fill(Psi_phi_track);
            if(eta_track > 4.0) hPhi_Psi_pos_eta->Fill(Psi_phi_track);
            if(eta_track < -4.0) hPhi_Psi_neg_eta->Fill(Psi_phi_track);
            hPhi    ->Fill(phi_track);

            // Calculate pT weight
            Double_t pT_weight = 1.0;
            if(pT_track < 2.0)  pT_weight = pT_track;
            if(pT_track >= 2.0) pT_weight = 2.0;
            //pT_weight = fabs(eta_track);


            // Event plane is assumed to be always at angle = 0, so Psi_phi_track is equal to phi of the particle
            iQxy[n_track][0] = pT_weight*TMath::Cos(2.0*phi_track); // iQx
            iQxy[n_track][1] = pT_weight*TMath::Sin(2.0*phi_track); // iQy

            // Define sub event and fill Qxy
            Int_t sub_event = 0;
            if( (n_track % 2) != 0) sub_event = 0;
            else sub_event = 1;

            if(fabs(eta_track) < 1.0)
            {
                Qxy[0][sub_event] += iQxy[n_track][0];
                Qxy[1][sub_event] += iQxy[n_track][1];
            }

            HA_track = HA_event.createTrack();
            HA_track->setiQx(iQxy[n_track][0]);
            HA_track->setiQy(iQxy[n_track][1]);
            HA_track->setpT(pT_track);
            HA_track->seteta(eta_track);
            HA_track->setv2_in(v2_track_smear);
            HA_track->setPhi_Psi_in(Psi_phi_track);
            HA_track->setPhi_in(phi_track);
            HA_track->setsub_event(sub_event);

            //cout << "n_track = " << n_track << ", Psi_phi_track = " << Psi_phi_track << endl;

        } // end track loop

        Double_t Psi_sub[2] = {0,0}; // [subA,subB]
        Double_t Psi    = calc_phi_event_plane_2nd(Qxy[0][0]+Qxy[0][1],Qxy[1][0]+Qxy[1][1]);
        for(Int_t isub = 0; isub < 2; isub++)
        {
            Psi_sub[isub] = calc_phi_event_plane_2nd(Qxy[0][isub],Qxy[1][isub]);
        }

        HA_event.setPsi_out(Psi);
        HA_event.setPsi_subA(Psi_sub[0]);
        HA_event.setPsi_subB(Psi_sub[1]);
        HA_event.setQx_subA(Qxy[0][0]);
        HA_event.setQy_subA(Qxy[1][0]);
        HA_event.setQx_subB(Qxy[0][1]);
        HA_event.setQy_subB(Qxy[1][1]);

        Tree_EP_Det_Sim->Fill();

    } // end event loop
    //------------------------------------------------------------------------------------------------------------------------------------


    return 1;
}



Int_t EP_Det_Sim::Finish()
{
    //-------------------------------------------------------
    Outputfile->cd();

    Tree_EP_Det_Sim->Write("",TObject::kOverwrite);
    hPhi_Psi ->Write();
    hPhi_Psi_pos_eta ->Write();
    hPhi_Psi_neg_eta ->Write();
    hPhi     ->Write();

    Outputfile->Close();
    cout << "Saved data to file" << endl;
    //-------------------------------------------------------

    return 1;

}
