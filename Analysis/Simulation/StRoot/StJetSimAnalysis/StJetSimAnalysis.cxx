#include "StJetSimAnalysis.h"

//----------------------------------------------------------------------------------------
static TRandom ran;
static TRandom3 r3;
static TString HistName;
static char NoP[50];
static TRandom ran_gen;
static Double_t N_Events_total;
static const Double_t particle_mass = 0.139;
static const Int_t Remove_N_hardest = 0;
static const Double_t min_pt_threshold = 0.0;
static const Double_t max_pt_threshold = 0.5;
static const Double_t jet_delta_eta_cut    = 100.0;
static const Double_t jet_delta_phi_cut    = 45.0*(2.0*Pi/360.0);
static Double_t area_cut;
//----------------------------------------------------------------------------------------


static const char HA_EVENT_TREE[]   = "HA_Events";
static TTree         *Tree_EP_Det_Sim;
static EP_Det_Sim_Event *HA_event;
static EP_Det_Sim_Event *HA_event_ptr;
static EP_Det_Sim_Track *HA_track;

static char* NAME_HA_EVENT_TREE   = "Tree_EP_Det_Sim_tree";
static char* NAME_HA_EVENT_BRANCH = "Tree_EP_Det_Sim_branch";

static TNtuple *NT_ReCoil_Jet;
static Float_t ReCoil_Jet_NTDataArray[14];


#include "StJetSimAnalysis_Func.h"


//------------------------------------------------------------------------------------------------------------------
ClassImp(StJetSimAnalysis)
    StJetSimAnalysis::StJetSimAnalysis()
{

}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
StJetSimAnalysis::~StJetSimAnalysis()
{

}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetSimAnalysis::Init()
{
    cout << "Init started" << endl;
    r3.SetSeed(0);

    area_cut = 0.4*TMath::Pi()*Jet_R*Jet_R;
    h_jet_pt_sub                 = new TH1D("h_jet_pt_sub","h_jet_pt_sub",1400,-40,100.0);
    h_Delta_phi                  = new TH1D("h_Delta_phi","h_Delta_phi",100,-TMath::Pi(),TMath::Pi());
    h_Delta_phi_trigger          = new TH1D("h_Delta_phi_trigger","h_Delta_phi_trigger",100,-TMath::Pi(),TMath::Pi());
    h_jet_pt_vs_Delta_phi        = new TH2D("h_jet_pt_vs_Delta_phi","h_jet_pt_vs_Delta_phi",50,-TMath::Pi(),TMath::Pi(),50,-20,30);
    h_jet_recoil_pt_vs_Delta_phi = new TH2D("h_jet_recoil_pt_vs_Delta_phi","h_jet_recoil_pt_vs_Delta_phi",50,-TMath::Pi(),TMath::Pi(),50,-20,30);
    h_jet_pt_vs_trigger_pt       = new TH2D("h_jet_pt_vs_trigger_pt","h_jet_pt_vs_trigger_pt",50,0,3.5,50,-20,30);
    h_track_pt_vs_delta_phi      = new TH2D("h_track_pt_vs_delta_phi","h_track_pt_vs_delta_phi",50,-TMath::Pi(),TMath::Pi(),50,0,3.5);
    h_trigger_pt_vs_delta_phi    = new TH2D("h_trigger_pt_vs_delta_phi","h_trigger_pt_vs_delta_phi",50,-TMath::Pi(),TMath::Pi(),50,0,3.5);

    //----------------------------------------------------------------------------------------------------
    // Simulation input
    TString Input_list = SEList;
    if (!Input_list.IsNull())   // if input file is ok
    {
        cout << "Open file list " << Input_list << endl;
        ifstream in(Input_list);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_chain  = new TChain( NAME_HA_EVENT_TREE, NAME_HA_EVENT_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = eIndir;
                    addfile += str;
                    Long64_t file_entries;
                    input_chain ->AddFile(addfile.Data(),-1, NAME_HA_EVENT_TREE );
                    file_entries = input_chain->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
        }
    }

    // Set the input tree
    if (!input_chain->GetBranch( NAME_HA_EVENT_BRANCH ))
    {
        cerr << "ERROR: Could not find branch '"
            << NAME_HA_EVENT_BRANCH << "'in tree!" << endl;
    }

    cout << "Set branch address" << endl;
    input_chain ->SetBranchAddress( NAME_HA_EVENT_BRANCH, &HA_event );
    cout << "Branch address was set" << endl;

    N_Events_total = input_chain->GetEntriesFast();
    cout << "Number of events in file(s) = " << N_Events_total << endl;
    if( (eStartEvent+N_Events) > N_Events_total) N_Events = (N_Events_total-eStartEvent);
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    TString outputfile_name = eOutdir;
    outputfile_name += ListName;
    outputfile_name += "_out.root";
    cout << "Output file: " << outputfile_name.Data() << endl;
    Outputfile        = new TFile(outputfile_name.Data(),"RECREATE");
    //----------------------------------------------------------------------------------------------------

    NT_ReCoil_Jet = new TNtuple("NT_ReCoil_Jet", 
                                "NT_ReCoil_Jet Ntuple", 
                                "EventId:JetId:rho:area:Jetphi:Jeteta:Jetpt:TrackId:eta:phi:pt:x:y:z");
    NT_ReCoil_Jet->SetAutoSave( 5000000 );
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetSimAnalysis::Make()
{
    cout << "Make started" << endl;
    r3.SetSeed(0);
    gRandom->SetSeed(0);
    cout << "Seed = " << r3.GetSeed() << endl;

    ran_gen.SetSeed(0);
    ran.SetSeed(0);

    input_chain->GetEntry( 0 );

    //-----------------------------------------------------------------------------
    cout << "Start of event loop" << endl;
    for(Long64_t n_evt = eStartEvent; n_evt < (eStartEvent+N_Events); n_evt++) // loop over all events
    {
        //cout << "n_evt = " << n_evt << endl;
        if (n_evt != 0  &&  n_evt % 50 == 0)
            cout << "." << flush;
        if (n_evt != 0  &&  n_evt % 500 == 0)
        {
            Double_t event_percent = 100.0*(n_evt-eStartEvent)/N_Events;
            cout << " " << n_evt << " (" << event_percent << "%) " << "\n" << "==> Processing data, " << flush;
        }

        if (!input_chain->GetEntry( n_evt )) // take the event -> information is stored in event
            break;  // end of data chunk

        Int_t   N_Tracks   = (Int_t)HA_event->getN_Tracks();
        Int_t   N_Cent     = (Int_t)HA_event->getN_Cent();
        Double_t v2_smear  = (Double_t)HA_event->getv2_smear();
        Double_t Psi_angle = (Double_t)HA_event->getPsi_in();
        Double_t Psi_rec   = (Double_t)HA_event->getPsi_out();

        if(N_Cent >= 2) continue; // only 0-10%

        std::vector<PseudoJet> particles_jet;
        std::vector< std::vector<Double_t> > vec_trigger_tracks; // eta, phi, pT, Delta phi

        for(Int_t n_track = 0; n_track < N_Tracks; n_track++) // track loop
        {
            HA_track = HA_event->getTrack( n_track ); // take the track
            Double_t pT_track             = (Double_t)HA_track->getpT();
            Double_t eta                  = (Double_t)HA_track->geteta();
            Double_t v2_track_smear       = (Double_t)HA_track->getv2_in();
            Double_t Psi_phi_track        = (Double_t)HA_track->getPhi_Psi_in();
            Double_t phi_track            = (Double_t)HA_track->getPhi_in();
            Double_t theta_track          = 2.0*TMath::ATan(TMath::Exp(-eta));

            Double_t Delta_phi_angle = Psi_angle-phi_track;
            if(Delta_phi_angle > TMath::Pi())  Delta_phi_angle -= 2.0*TMath::Pi();
            if(Delta_phi_angle < -TMath::Pi()) Delta_phi_angle += 2.0*TMath::Pi();
            h_Delta_phi ->Fill(Delta_phi_angle);
            h_track_pt_vs_delta_phi ->Fill(Delta_phi_angle,pT_track);

            // Fill trigger track candidates
            if(pT_track > min_pt_threshold && pT_track < max_pt_threshold && fabs(eta) < 1.0)
            {
                std::vector<Double_t> vec_in;
                vec_in.resize(4);
                vec_in[0] = eta;
                vec_in[1] = phi_track;
                vec_in[2] = pT_track;
                vec_in[3] = Delta_phi_angle;
                vec_trigger_tracks.push_back(vec_in);
                //cout << "n_evt: " << n_evt << ", pT_track: " << pT_track << endl;
                
            }

            TLorentzVector TLV_Particle_use;
            TLV_Particle_use.SetPtEtaPhiM(pT_track,eta,phi_track,particle_mass);
            PseudoJet Fill_PseudoJet(TLV_Particle_use.Px(),TLV_Particle_use.Py(),TLV_Particle_use.Pz(),TLV_Particle_use.E());
            particles_jet.push_back(Fill_PseudoJet);

            //cout << "n_track: " << n_track << ", pT_track: " << pT_track << endl;
        } // End of track loop



        //-----------------------------------------------------------------------------
        // Take one random track from the trigger track candidates
#if 0
        for(Int_t i_val = 0; i_val < vec_trigger_tracks.size(); i_val++) cout << "i before: " << i_val << ", eta: " << vec_trigger_tracks[i_val][0] << ", phi: " << vec_trigger_tracks[i_val][1] << ", pT: " << vec_trigger_tracks[i_val][2] << endl;
#endif
        // sort the pT values
        std::sort (vec_trigger_tracks.begin(), vec_trigger_tracks.end(), sortFunc);

#if 0
        for(Int_t i_val = 0; i_val < vec_trigger_tracks.size(); i_val++) cout << "i after: " << i_val << ", eta: " << vec_trigger_tracks[i_val][0] << ", phi: " << vec_trigger_tracks[i_val][1] << ", pT: " << vec_trigger_tracks[i_val][2] << endl;
#endif

        Int_t trigger_index_min   = 0;
        Int_t trigger_index_range = vec_trigger_tracks.size()-trigger_index_min;

        Int_t random_trigger_index = ran_gen.Integer(trigger_index_range); // [0,trigger_index_range-1]
        random_trigger_index += trigger_index_min;

        if(random_trigger_index >= vec_trigger_tracks.size())
        {
            //cout << "ERROR: random_trigger_index out of range!" << endl;
            continue;
        }
        h_trigger_pt_vs_delta_phi ->Fill(vec_trigger_tracks[random_trigger_index][3],vec_trigger_tracks[random_trigger_index][2]);
        h_Delta_phi_trigger       ->Fill(vec_trigger_tracks[random_trigger_index][3]);
        //cout << "n_evt: " << n_evt << ", trigger pt: " << vec_trigger_tracks[random_trigger_index][2] << ", trigger Delta phi: " << vec_trigger_tracks[random_trigger_index][3] << endl;
        //-----------------------------------------------------------------------------



        //-----------------------------------------------------------------------------
        // Start of jet reconstruction

        vector<PseudoJet> jets;
        // choose a jet definition
        JetDefinition jet_def(antikt_algorithm, Jet_R);

        // jet area definition
        Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
        GhostedAreaSpec area_spec(ghost_maxrap);
        AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

        ClusterSequenceArea clust_seq_hard(particles_jet, jet_def, area_def);

        double ptmin = 0.2;
        vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
        Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - Jet_R); // Fiducial cut for jets
        jets = Fiducial_cut_selector(jets_all);

        // background estimation
        JetDefinition jet_def_bkgd(kt_algorithm, Jet_R); // <--
        AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
        Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest)); // <--

        JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
        Subtractor subtractor(&bkgd_estimator);
        bkgd_estimator.set_particles(particles_jet);

        Double_t jet_rho   = bkgd_estimator.rho();
        Double_t jet_sigma = bkgd_estimator.sigma();
        //cout << "jet_rho: " << jet_rho << endl;
        //-----------------------------------------------------------------------------



        //-----------------------------------------------------------------------------
        // Jet loop
        if(vec_trigger_tracks.size() > 0)
        {
            Double_t trigger_pt = vec_trigger_tracks[random_trigger_index][2];
            for(Int_t i_jet = 0; i_jet < jets.size(); i_jet++)
            {
                Float_t jet_pt        = jets[i_jet].perp();
                Float_t jet_area      = jets[i_jet].area();
                Float_t jet_pt_sub    = jets[i_jet].perp() - jet_rho*jet_area;
                Float_t jet_eta       = jets[i_jet].eta();
                Float_t jet_phi_jet   = jets[i_jet].phi(); // 0..2pi

                //cout << "i_jet: " << i_jet << ", jet_eta: " << jet_eta << ", jet_phi_jet: " << jet_phi_jet << endl;

                Float_t jet_delta_eta = fabs(jet_eta + vec_trigger_tracks[random_trigger_index][0]);
                Float_t jet_delta_phi = fabs(jet_phi_jet - vec_trigger_tracks[random_trigger_index][1]);
                if(jet_delta_phi > 2.0*Pi)  jet_delta_phi -= 2.0*Pi;

                Float_t dijet_delta_phi = jet_phi_jet-vec_trigger_tracks[random_trigger_index][1]; // -2*Pi..2*Pi
                if(dijet_delta_phi < 0.0) dijet_delta_phi += 2.0*Pi; // 0..2*Pi
                if(dijet_delta_phi > 1.5*Pi)
                {
                    dijet_delta_phi = -0.5*Pi + (dijet_delta_phi-1.5*Pi); // -0.5*Pi..1.5*Pi
                }

                Double_t jet_EP_Delta_phi = Psi_angle - jet_phi_jet;
                if(jet_EP_Delta_phi > Pi)  jet_EP_Delta_phi -= 2.0*Pi;
                if(jet_EP_Delta_phi < -Pi) jet_EP_Delta_phi += 2.0*Pi;

                if(jet_area > area_cut)
                {
                    h_jet_pt_vs_Delta_phi   ->Fill(jet_EP_Delta_phi,jet_pt_sub);
                    h_jet_pt_vs_trigger_pt  ->Fill(trigger_pt,jet_pt_sub);
                }

                //--------------------------------------------------
                // Recoil jet histograms
                if(
                   jet_delta_eta < jet_delta_eta_cut
                   && fabs(Pi-jet_delta_phi) < jet_delta_phi_cut
                  )
                {
                    if(jet_area > area_cut)
                    {
                        h_jet_pt_sub                 ->Fill(jet_pt_sub);
                        h_jet_recoil_pt_vs_Delta_phi ->Fill(jet_EP_Delta_phi,jet_pt_sub);

                        ReCoil_Jet_NTDataArray[0] = (Float_t)n_evt;         // EventId
                        ReCoil_Jet_NTDataArray[0] = (Float_t)i_jet;         // JetId
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_rho;       // rho
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_area;      // area
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_phi;       // Jetphi
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_eta;       // Jeteta
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_pt_sub;    // Jetpt
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_const_eta; // TrackId
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_const_phi; // eta
                        ReCoil_Jet_NTDataArray[0] = (Float_t)jet_const_pt;  // phi
                         

                        #if 0
                        cout << "i_jet: " << i_jet << ", jet_delta_phi: " << jet_delta_phi*TMath::RadToDeg()
                            << ", trigger phi: " << vec_trigger_tracks[random_trigger_index][1]*TMath::RadToDeg()
                            << ", jet phi: " << jet_phi_jet*TMath::RadToDeg() << endl;
                        #endif
                    }
                }
                //--------------------------------------------------

            } // End of jet loop
        }
        //-----------------------------------------------------------------------------



    } // End of event loop
    //-----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void StJetSimAnalysis::Finish()
{
    cout << "Finish started" << endl;
    Outputfile   ->cd();
    h_jet_pt_sub                 ->Write();
    h_Delta_phi                  ->Write();
    h_Delta_phi_trigger          ->Write();
    h_jet_pt_vs_Delta_phi        ->Write();
    h_jet_pt_vs_trigger_pt       ->Write();
    h_track_pt_vs_delta_phi      ->Write();
    h_trigger_pt_vs_delta_phi    ->Write();
    h_jet_recoil_pt_vs_Delta_phi ->Write();
    Outputfile   ->Close();
}
//------------------------------------------------------------------------------------------------------------------

