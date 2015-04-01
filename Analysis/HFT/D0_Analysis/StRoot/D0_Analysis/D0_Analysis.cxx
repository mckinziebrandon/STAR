//#include <assert.h>
#include <fstream>
#include <string>
#include "StMessMgr.h"
#include "D0_Analysis.h"
#include "TString.h"
#include "TLatex.h"
#include "TPad.h"

ClassImp(D0_Analysis)

using std::ifstream ;
using std::string ;
using std::vector ;

const Int_t N_Analysis = 1;

static const TString V0_TREE[N_Analysis]   = {"D0Events"};
static const TString V0_BRANCH[N_Analysis] = {"Events"};

static char* D0_TREE;
static char* D0_BRANCH;

static Int_t SE_input_flag = 1;
static Int_t ME_input_flag = 0;

static TString HistName;

const Int_t N_h_InvMass    = 6;
const Int_t N_h_InvMass_pt = 6;
static TH1D* h_InvMass[N_h_InvMass][N_h_InvMass][N_h_InvMass][N_h_InvMass][N_h_InvMass][N_h_InvMass_pt];
static TH1D* special_h_InvMass = new TH1D("special_h_InvMass", "special_h_InvMass", 100, 1.4, 2.4);


//____________________________________________________________________________________________________
// Default constructor
D0_Analysis::D0_Analysis()
{
    clear() ;
}

//____________________________________________________________________________________________________
// Default destructor
D0_Analysis::~D0_Analysis()
{
}

//____________________________________________________________________________________________________
void D0_Analysis::clear()
{

}

void D0_Analysis::setInputDir(const TString inputdir)
{
    pinputdir = inputdir.Copy();
    cout << "Input directory was set to: " << pinputdir.Data() << endl;
}
void D0_Analysis::setOutputfile(const TString outputfile)
{
    poutputfile = outputfile.Copy();
    cout << "Output file was set to: " << poutputfile.Data() << endl;
}
void D0_Analysis::setSEList(const TString iSEList)
{
    SEList = iSEList.Copy();
    cout << "Same event list was set to: " << SEList.Data() << endl;
}
void D0_Analysis::setStopEvent_SE(const Long64_t StopEvent_SE)
{
    nStopEvent_SE = StopEvent_SE;
    cout << "nStopEvent_SE = " << nStopEvent_SE << endl;
}
void D0_Analysis::setStartEvent_SE(const Long64_t StartEvent_SE)
{
    nStartEvent_SE = StartEvent_SE;
    cout << "nStartEvent_SE = " << nStartEvent_SE << endl;
}

void D0_Analysis::init()
{

    cout << "Initializing parameters and input/output" << endl;
    Outputfile = new TFile(poutputfile.Data(),"RECREATE");

    D0_TREE       = (char*)V0_TREE[eAnalysis].Data();
    D0_BRANCH     = (char*)V0_BRANCH[eAnalysis].Data();



    //----------------------------------------------------------------------------------------------------
	cout << "Define histograms" << endl;
	for(Int_t A = 0; A < N_h_InvMass; A++){
	for(Int_t B = 0; B < N_h_InvMass; B++){
	for(Int_t X = 0; X < N_h_InvMass; X++){
	for(Int_t Y = 0; Y < N_h_InvMass; Y++){
	for(Int_t AB = 0; AB < N_h_InvMass; AB++){
	for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){
			HistName = "h_InvMass_A";
			HistName += A;
            HistName += "_B";
			HistName += B;
            HistName += "_X";
            HistName += X;
            HistName += "_Y";
            HistName += Y;
            HistName += "_AB";
            HistName += AB;
			HistName += "_pt";
			HistName += i_hist_pt;
			h_InvMass[A][B][X][Y][AB][i_hist_pt] = 
				new TH1D(HistName.Data(),HistName.Data(),100,1.4,2.4);
	}}}}}}


//----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( D0_TREE, D0_TREE );
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
                    input_SE ->AddFile(addfile.Data(),-1, D0_TREE );
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
    if (SE_input_flag == 1 && !input_SE->GetBranch( D0_BRANCH ))
    {
        cerr << "ERROR: Could not find branch '"
            << D0_BRANCH << "'in tree!" << endl;
    }

    D0_event = new StD0Event();

    if(SE_input_flag == 1)
    {
        {
            input_SE  ->SetBranchAddress( D0_BRANCH, &D0_event );

            num_events_SE = input_SE->GetEntriesFast();
            cout << "Number of events in file(s) = " << num_events_SE << endl;
            if(nStartEvent_SE > num_events_SE) nStartEvent_SE = num_events_SE;
            if(nStopEvent_SE  > num_events_SE) nStopEvent_SE  = num_events_SE;

            cout << "New nStartEvent_SE = " << nStartEvent_SE << ", new nStopEvent_SE = " << nStopEvent_SE << endl;
        }
    }
    //----------------------------------------------------------------------------------------------------

}

void D0_Analysis::loop()
{
    //--------------------------------
    // Track properties
    Float_t m2A; // mass2 of particle A
    Float_t m2B;
    Float_t nsA; // nsigma dE/dx of particle A
    Float_t nsB;
    Float_t dcaA; // distance of closest approach of particle A
    Float_t dcaB;
    Float_t iQxA; // Q-vector x-component of particle A
    Float_t iQyA;
    Float_t iQxB;
    Float_t iQyB;
    Float_t etaA;  // pseudo rapidity of particle A
    Float_t etaB;

    Float_t InvAB;  // invariant mass of particles A and B
    Float_t p_t; // transverse momentum of particles AB(C)
    Float_t rap; // rapidity of particles AB(C)
    Float_t phi; // azimuth angle of particles AB(C)
    Float_t theta; // polar angle of particles AB(C)

    Float_t qpA; // momentum times charge
    Float_t qpB;

    Float_t VerdistX;
    Float_t VerdistY;
    Float_t dcaAB;

    // Event based properties
    Float_t EventVertexX;
    Float_t EventVertexY;
    Float_t EventVertexZ;
    Int_t   RunId;
    Float_t refMult;
    Int_t   n_prim;
    Int_t   n_non_prim;
    Int_t   n_tof_prim;
    Float_t EP_Qx_eta_pos_ptw;
    Float_t EP_Qy_eta_pos_ptw;
    Float_t EP_Qx_eta_neg_ptw;
    Float_t EP_Qy_eta_neg_ptw;
    Float_t EP_Qx_ptw;
    Float_t EP_Qy_ptw;
    Int_t   Qtracks_eta_pos;
    Int_t   Qtracks_eta_neg;
    Int_t   Qtracks_full;

    Float_t ZDCx;
    Float_t BBCx;
    Float_t vzVpd;

    UShort_t      fNumTracks;
    //--------------------------------


    //----------------------------------------------------------------------------------------------------
    // Event loop
    for(Int_t SE_ME_loop = 0; SE_ME_loop < 1; SE_ME_loop++) // 0 = same event, 1 = mixed event
    {
        Long64_t start_event_use = 0;
        Long64_t stop_event_use  = 0;
        if(SE_ME_loop == 0 && SE_input_flag == 1)
        {
            input_SE  ->SetBranchAddress( D0_BRANCH, &D0_event );

            cout << "" << endl;
            cout << "------------------- Start looping: Same event -------------------" << endl;
            start_event_use = nStartEvent_SE;
            stop_event_use  = nStopEvent_SE;
            cout << "start_event_use = " << start_event_use << ", stop_event_use = " << stop_event_use << endl;
            input_SE->GetEntry( 0 ); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
        }

        cout << "SE_ME_loop = " << SE_ME_loop << ", start_event_use = " << start_event_use << ", stop_event_use = " << stop_event_use << ", SE_input_flag = " << SE_input_flag << endl;

        for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
        {
            if (counter != 0  &&  counter % 100000 == 0)
                cout << "." << flush;
            if (counter != 0  &&  counter % 1000000 == 0)
            {
                if((stop_event_use-start_event_use) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                    cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (D0 Analysis) " << flush;
                }
            }

            if(SE_ME_loop == 0)
            {
                if (!input_SE->GetEntry( counter )) // take the event -> information is stored in event
                    break;  // end of data chunk
            }

            //---------------------------------------------------------------------------
            fNumTracks           = D0_event->getNumTracks(); // number of tracks in this event
            EventVertexX         = D0_event->getx();
            EventVertexY         = D0_event->gety();
            EventVertexZ         = D0_event->getz();
            RunId                = D0_event->getid();
            refMult              = D0_event->getmult();
            n_prim               = D0_event->getn_prim();
            n_non_prim           = D0_event->getn_non_prim();
            n_tof_prim           = D0_event->getn_tof_prim();
            EP_Qx_eta_pos_ptw    = D0_event->getEP_Qx_eta_pos_ptw();
            EP_Qy_eta_pos_ptw    = D0_event->getEP_Qy_eta_pos_ptw();
            EP_Qx_eta_neg_ptw    = D0_event->getEP_Qx_eta_neg_ptw();
            EP_Qy_eta_neg_ptw    = D0_event->getEP_Qy_eta_neg_ptw();
            EP_Qx_ptw            = D0_event->getEP_Qx_ptw();
            EP_Qy_ptw            = D0_event->getEP_Qy_ptw();
            Qtracks_eta_pos      = D0_event->getQtracks_eta_pos();
            Qtracks_eta_neg      = D0_event->getQtracks_eta_neg();
            Qtracks_full         = D0_event->getQtracks_full();
            ZDCx                 = D0_event->getZDCx();
            BBCx                 = D0_event->getBBCx();
            vzVpd                = D0_event->getvzVpd();

            //----------------------------------------------------------------------------------------------------
            if(
               fabs(EventVertexZ) < 20.0
              )
            {
                for(UShort_t i_track = 0; i_track < fNumTracks; ++i_track) // loop over all tracks of the actual event
                {
                    D0_track      = D0_event->getTrack( i_track ); // take the track
                    m2A                     = D0_track->getm2A();  // squared mass of Kaon candidate
                    m2B                     = D0_track->getm2B();  // squared mass of Pion candidate
                    nsA                     = D0_track->getnsA();  // n*sigma (dE/dx) of Kaon candidate
                    nsB                     = D0_track->getnsB();  // n*sigma (dE/dx) of Pion candidate
                    dcaA                    = D0_track->getdcaA(); // distance of closest (dca) approach of Kaon candidate to the primary vertex
                    dcaB                    = D0_track->getdcaB(); // distance of closest (dca) approach of Pion candidate to the primary vertex
                    iQxA                    = D0_track->getiQxA();
                    iQyA                    = D0_track->getiQyA();
                    iQxB                    = D0_track->getiQxB();
                    iQyB                    = D0_track->getiQyB();
                    etaA                    = D0_track->getetaA(); // pseudo rapdity of the Kaon candidate
                    etaB                    = D0_track->getetaB(); // pseudo rapdity of the Pion candidate
                    InvAB                   = D0_track->getInvAB(); // invariant mass of Kaon and Pion candidates
                    p_t                     = D0_track->getpt(); // invariant momentum of Kaon + Pion
                    rap                     = D0_track->getrap(); // true Rapidity, not pseudo-Rapidity of Kaon + Pion!
                    phi                     = D0_track->getphi(); // azimuthal angle of Kaon + Pion
                    theta                   = D0_track->gettheta(); // polar angle of Kaon + Pion
                    qpA                     = D0_track->getqpA(); // momentum*charge of Kaon candidate
                    qpB                     = D0_track->getqpB(); // momentum*charge of Pion candidate
                    VerdistX                = D0_track->getVerdistX(); // distance between primary and decay vertex
                    VerdistY                = D0_track->getVerdistY(); // distance of closest approach of mother particle to primary vertex
                    dcaAB                   = D0_track->getdcaAB(); // distance of closest approach between Kaon and Pion

		std::vector<Double_t> dcaA_cut, dcaB_cut, dcaAB_cut, VerdistX_cut, VerdistY_cut, pt_cut(N_h_InvMass_pt);
        pt_cut = {0, 0.5, 1, 1.5, 2.5, 5, 10};

		// initialize vectors with cut ranges
        for(Int_t i = 0; i < N_h_InvMass; i++){
            dcaA_cut.push_back(0.002 + i * 0.0033);		// 20 - 185 microns
            dcaB_cut.push_back(0.002 + i * 0.0033);		// 20 - 185 microns
            VerdistX_cut.push_back(0.005 + i * 0.0083);	// 50 - 465 microns
            VerdistY_cut.push_back(0.0185 - i * 0.0033);// 185 - 20 microns
            dcaAB_cut.push_back(0.0185 - i * 0.0033);	// 185 - 20 microns 
        }

        /* checking vector elements
        for(Int_t i = 0; i < N_h_InvMass; i++){
            std::cout << "dcaA_cut[" << i << "] = " << dcaA_cut[i];
            std::cout << ";; dcaB_cut[" << i << "] = " << dcaB_cut[i];
            std::cout << ";; dcaAB_cut[" << i << "] = " << dcaA_cut[i];
            std::cout << ";; VerdistX_cut[" << i << "] = " << VerdistX_cut[i];
            std::cout << ";; VerdistY_cut[" << i << "] = " << VerdistY_cut[i] << endl;
        }
        */

        // ensure all values are within cut range
        Bool_t in_cut_range = true;
        if(
            fabs(dcaA) < dcaA_cut.front()   ||
            fabs(dcaB) < dcaB_cut.front()   ||
            dcaAB      > dcaAB_cut.front()  || 
            VerdistX   < VerdistX_cut.front()   || 
            VerdistY   > VerdistY_cut.front()   || 
            p_t  > pt_cut.back()
        )   in_cut_range = false;

        // loop counter values are indices of cut vectors
        if (in_cut_range){
            for(Int_t A = 0; A < N_h_InvMass; A++){
                if (fabs(dcaA) < dcaA_cut[A]) break;
            for(Int_t B = 0; B < N_h_InvMass; B++){
                if (fabs(dcaB) < dcaB_cut[B]) break;
            for(Int_t X = 0; X < N_h_InvMass; X++){
                if (VerdistX   < VerdistX_cut[X]) break;
            for(Int_t Y = 0; Y < N_h_InvMass; Y++){
                if (VerdistY   > VerdistY_cut[Y]) break;
            for(Int_t AB = 0; AB < N_h_InvMass; AB++){
                if (dcaAB      > dcaAB_cut[AB]) break;
            for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){
                if(p_t  < pt_cut[i_hist_pt+1] &&
                    p_t >= pt_cut[i_hist_pt] ){
                        h_InvMass[A][B][X][Y][AB][i_hist_pt] ->Fill(InvAB);
                        break;
                }

            }}}}}}
        }

        // fill special hist
        if(
            fabs(dcaA) > 0.004  &&
            fabs(dcaB) > 0.004  &&
            dcaAB      < 0.01  &&
            VerdistX   > 0.02   &&
            VerdistY   < 0.01  &&
            fabs(qpA)  > 0.8    &&
            fabs(qpB)  > 0.8 
            ){
                special_h_InvMass->Fill(InvAB);
            }

                } // end of track loop
            }
            //----------------------------------------------------------------------------------------------------

        } // end of event loop
    } // end of same event <-> mixed event loop

    cout << "" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------------" << endl;


}

void D0_Analysis::finalize()
{

	cout << "" << endl;
	cout << "Write output" << endl;
	Outputfile      ->cd();
    special_h_InvMass->Write();
	Outputfile->mkdir("h_InvMass");
	Outputfile->cd("h_InvMass");

    Int_t n_histograms = 0;
	// it seems pretty inefficient to do all this looping again . . .
	for(Int_t A = 0; A < N_h_InvMass; A++){
        std::cout << "A = " << A << endl;
	for(Int_t B = 0; B < N_h_InvMass; B++){
        std::cout << "B = " << B << endl;
	for(Int_t X = 0; X < N_h_InvMass; X++){
	for(Int_t Y = 0; Y < N_h_InvMass; Y++){
	for(Int_t AB = 0; AB < N_h_InvMass; AB++){
	for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++)
	{
		    h_InvMass[A][B][X][Y][AB][i_hist_pt]->Write();
            n_histograms += 1;
    }
		
	}}}}}


    cout << "\n TOTAL NUMBER OF HISTOGRAMS = " << n_histograms << "\n\n";
	Outputfile      ->cd();
	cout << "Close output file" << endl;
	cout << "" << endl;
	Outputfile      ->Close();

}

