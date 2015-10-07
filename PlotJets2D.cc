
#include "functions.h"

void PlotJets2D(Long64_t event_plot = 0, Int_t flag_high_pT = 0, Int_t Radius = 0)
{
    cout << "PlotJets started" << endl;

    // 0,4,9

    //-------------------------------------------------
    cout << "Open inputfile" << endl;
    TString addfile;
    addfile = "test.root";

    TChain* input_SE  = new TChain( "NT_ReCoil_Jet" , "NT_ReCoil_Jet" );
    input_SE ->AddFile(addfile.Data(),-1, "NT_ReCoil_Jet" );

    cout << "Set branch addresses" << endl;
    Float_t EventId,JetId,rho,area,Jetphi,Jeteta,Jetpt,TrackId,eta,phi,pt,x,y,z;
    input_SE->SetBranchAddress("EventId",&EventId);
    input_SE->SetBranchAddress("JetId",&JetId);
    input_SE->SetBranchAddress("rho",&rho);
    input_SE->SetBranchAddress("area",&area);
    input_SE->SetBranchAddress("Jetphi",&Jetphi);
    input_SE->SetBranchAddress("Jeteta",&Jeteta);
    input_SE->SetBranchAddress("Jetpt",&Jetpt); // jet pt - rho*A
    input_SE->SetBranchAddress("TrackId",&TrackId);
    input_SE->SetBranchAddress("eta",&eta);
    input_SE->SetBranchAddress("phi",&phi);
    input_SE->SetBranchAddress("pt",&pt);
    input_SE->SetBranchAddress("x",&x);
    input_SE->SetBranchAddress("y",&y);
    input_SE->SetBranchAddress("z",&z);

    Long64_t N_entries = input_SE->GetEntries();
    cout << "Entries in tree = " << N_entries << endl;
    //-------------------------------------------------



    //-------------------------------------------------
    TH2D* h_2D_jet_eta_vs_phi         = new TH2D("h_2D_jet_eta_vs_phi","h_2D_jet_eta_vs_phi",64,-1.0,2.0*Pi+1.0,64,-1.5,1.5);
    TH2D* h_2D_jet_eta_vs_phi_no_fill = new TH2D("h_2D_jet_eta_vs_phi_no_fill","h_2D_jet_eta_vs_phi_no_fill",400,-1.0,2.0*Pi+1.0,200,-1.5,1.5);
    //-------------------------------------------------



    //-------------------------------------------------
    const Double_t radius_2D_jet = 4.0;
    //-------------------------------------------------



    //-------------------------------------------------
    cout << "Start event loop" << endl;

    const Int_t N_max_tracks = 10000;
    TPolyMarker* PM_jet_real_tracks[N_max_tracks];
    TPolyMarker* PM_jet_area_tracks[N_max_tracks];

    const Int_t N_max_jets = 120;
    TPolyMarker* PM_jet_center[N_max_jets];     // [jet index]
    Double_t Array_Jet_info[N_max_jets][3];     // [jet][eta,phi,pt]
    Double_t Array_track_info[N_max_tracks][3]; // [jet][eta,phi,pt]
    TString TS_jet_label[N_max_jets];           // [jet index]
    Double_t Jet_center_array[N_max_jets][2];   // [jet index][phi,eta]

    TLorentzVector   lorentz_vector;
    Double_t bField = 4.98948e-14; //magnetic field ?
    Double_t jet_delta_phi_cut = 45.0*(2.0*Pi/360.0);;

    Long64_t start_event_use = 0;
    Long64_t stop_event_use  = N_entries;
    Long64_t event_counter   = -1;//?
    Long64_t EventId_old     = -1;//?
    Long64_t track_counter   = 0;
    Long64_t track_real_use_counter = 0;
    Long64_t track_area_use_counter = 0;
    Color_t  track_color     = kAzure-2;

    Float_t leading_pt     = 0.0;
    Float_t leading_pt_eta = 0.0;
    Float_t leading_pt_phi = 0.0;
    Int_t   Jet_index      = -1;
    Int_t   JetId_old      = -1;

    Color_t jet_color[N_max_jets] = {2,kAzure-2,kGreen+1,kMagenta+2,kCyan+2,kOrange+7,kTeal-8,kYellow-2,kBlue-9,kPink+7,kCyan-3
    ,kRed+2,kViolet-2,kSpring-1,kYellow-6,kRed-2,kGreen-10,kBlue-3,kOrange-5,kCyan+4,kRed,kAzure-3,kGreen+2,kMagenta+3,kCyan+3,kOrange+8,kTeal-9,kYellow-3,
    kBlue-8,kPink+8,kCyan-4,kRed+3,kViolet-3,kSpring-2,kYellow-7,kRed-3,kGreen-9,kBlue-4,kOrange-5,kCyan+5,
    2,kAzure-2,kGreen+1,kMagenta+2,kCyan+2,kOrange+7,kTeal-8,kYellow-2,kBlue-9,kPink+7,kCyan-3
    ,kRed+2,kViolet-2,kSpring-1,kYellow-6,kRed-2,kGreen-10,kBlue-3,kOrange-5,kCyan+4,kRed,kAzure-3,kGreen+2,kMagenta+3,kCyan+3,kOrange+8,kTeal-9,kYellow-3,
    kBlue-8,kPink+8,kCyan-4,kRed+3,kViolet-3,kSpring-2,kYellow-7,kRed-3,kGreen-9,kBlue-4,kOrange-5,kCyan+5,
    2,kAzure-2,kGreen+1,kMagenta+2,kCyan+2,kOrange+7,kTeal-8,kYellow-2,kBlue-9,kPink+7,kCyan-3
    ,kRed+2,kViolet-2,kSpring-1,kYellow-6,kRed-2,kGreen-10,kBlue-3,kOrange-5,kCyan+4,kRed,kAzure-3,kGreen+2,kMagenta+3,kCyan+3,kOrange+8,kTeal-9,kYellow-3,
    kBlue-8,kPink+8,kCyan-4,kRed+3,kViolet-3,kSpring-2,kYellow-7,kRed-3,kGreen-9,kBlue-4,kOrange-5,kCyan+5};

    TPolyMarker* PM_jet_real_tracks_trig = new TPolyMarker();



    //-------------------------------------------------------------------------------
    Int_t event_counter_high_pT = 0;
    for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
    {
        if (counter != 0  &&  counter % 1000 == 0)
            cout << "." << flush;
        if (counter != 0  &&  counter % 10000 == 0)
        {
            if((stop_event_use-start_event_use) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (jet) " << flush;
            }
        }

        if (!input_SE->GetEntry( counter )) // take the event -> information is stored in event
            break;  // end of data chunk


        if((Long64_t)EventId != EventId_old)
        {
            event_counter++;
            EventId_old = EventId;
        }
        if(pt >= 9.0)
        {
            cout << "Event with a high pt track found: " << EventId << endl;
            event_counter_high_pT++;
            if(event_counter_high_pT > event_plot) break;
        }
    }
    if(flag_high_pT) event_plot = event_counter;

    event_counter   = -1;
    EventId_old     = -1;
    //-------------------------------------------------------------------------------



    //-------------------------------------------------------------------------------
    for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
    {
        if (counter != 0  &&  counter % 1000 == 0)
            cout << "." << flush;
        if (counter != 0  &&  counter % 10000 == 0)
        {
            if((stop_event_use-start_event_use) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (jet) " << flush;
            }
        }

        if (!input_SE->GetEntry( counter )) // take the event -> information is stored in event
            break;  // end of data chunk


        if(event_counter > event_plot) break;

        if((Long64_t)EventId != EventId_old)
        {
            event_counter++;
            EventId_old = EventId;
            track_counter = 0;
        }

        
        // if we wanted to see nth jet event and this is the nth jet event
        if(event_counter == event_plot && track_counter == 0)
        {
            leading_pt     = 0.0;
            leading_pt_phi = 0.0;
            Jet_index      = -1;

            cout << "Get highest pt particle in event for event = " << event_counter << ", counter = " << counter << ", Event_Id = " << EventId << endl;
            Int_t   leading_pt_counter = 0;
            for(Long64_t counterB = counter; counterB < stop_event_use; counterB++)
            {
                input_SE->GetEntry( counterB );
                //cout << "counterB = " << counterB << ", EventId = " << EventId << ", EventId_old = " << EventId_old << endl;
                if((Long64_t)EventId != EventId_old) break;
                if(fabs(pt) > leading_pt)
                {
                    leading_pt = fabs(pt);
                    leading_pt_counter = counterB;
                    cout << "leading_pt_counter = " << leading_pt_counter << ", leading_pt = " << leading_pt << endl;
                }
            }
            input_SE->GetEntry( leading_pt_counter );
            leading_pt_phi = phi;
            leading_pt_eta = eta;
            cout << "leading_pt_phi = " << leading_pt_phi << endl;

            PM_jet_real_tracks_trig->SetNextPoint(leading_pt_phi,leading_pt_eta);
            PM_jet_real_tracks_trig->SetMarkerStyle(20);
            PM_jet_real_tracks_trig->SetMarkerColor(2);
            PM_jet_real_tracks_trig->SetMarkerSize(0.5*fabs(pt)+0.6);

            input_SE->GetEntry( counter );
        }

        Float_t jet_delta_phi = fabs(Jetphi - leading_pt_phi);

        if(jet_delta_phi > 2.0*Pi)  jet_delta_phi -= 2.0*Pi;
        //if(jet_delta_phi > Pi) jet_delta_phi -= Pi;

        Int_t flag_plot_jet = 0;
        if(
           event_counter == event_plot
           //&& (Int_t)JetId == 0
           //&& fabs(Pi-jet_delta_phi) < jet_delta_phi_cut
           //&& fabs(pt) > 0.5
           //&& Jet_index < N_max_jets
          )
        {
            if((Int_t)JetId != JetId_old)
            {
                Jet_index++;
                JetId_old = (Int_t)JetId;
                flag_plot_jet = 1;

                PM_jet_center[Jet_index] = new TPolyMarker();
                PM_jet_center[Jet_index]->SetNextPoint(Jetphi,Jeteta);
                PM_jet_center[Jet_index]->SetMarkerStyle(5);
                //PM_jet_center[Jet_index]->SetMarkerColor(jet_color[Jet_index]);
                PM_jet_center[Jet_index]->SetMarkerColor(1);
                PM_jet_center[Jet_index]->SetMarkerSize(1.0);

                Jet_center_array[Jet_index][0] = Jetphi;
                Jet_center_array[Jet_index][1] = Jeteta;

                Array_Jet_info[Jet_index][0] = Jeteta;
                Array_Jet_info[Jet_index][1] = Jetphi;
                Array_Jet_info[Jet_index][2] = Jetpt;

                sprintf(NoP,"%3.1f",Jetpt);
                HistName = NoP;
                HistName += "GeV/c";
                TS_jet_label[Jet_index] = HistName;
                //cout << "Jet_index: " << Jet_index << ", Jetpt: " << Jetpt << endl;
            }
            else
            {
                flag_plot_jet = 0;
            }

            if(fabs(pt) > 0.15)
            {
#if 0
                cout << "track_counter: " << track_counter << ", pt: " << pt << ", Jetpt_sub: " << Jetpt
                    << ", jet_pt: " << Jetpt+(rho*area) << ", Jetphi: " << Jetphi
                    << ", EventId: " << EventId << ", rho: " << rho << endl;
#endif

                h_2D_jet_eta_vs_phi->Fill(phi,eta,fabs(pt));
                if(track_real_use_counter < N_max_tracks)
                {
                    PM_jet_real_tracks[track_real_use_counter] = new TPolyMarker();
                    PM_jet_real_tracks[track_real_use_counter]->SetNextPoint(phi,eta);
                    PM_jet_real_tracks[track_real_use_counter]->SetMarkerStyle(20);
                    PM_jet_real_tracks[track_real_use_counter]->SetMarkerColor(jet_color[Jet_index]);
                    PM_jet_real_tracks[track_real_use_counter]->SetMarkerSize(0.6*TMath::Power(fabs(pt),0.5)+0.8);

                    Array_track_info[track_real_use_counter][0] = eta;
                    Array_track_info[track_real_use_counter][1] = phi;
                    Array_track_info[track_real_use_counter][2] = fabs(pt);
                }
                track_real_use_counter++;
            }
            else
            {
                if(track_area_use_counter < N_max_tracks)
                {
                    PM_jet_area_tracks[track_area_use_counter] = new TPolyMarker();
                    PM_jet_area_tracks[track_area_use_counter]->SetNextPoint(phi,eta);
                    PM_jet_area_tracks[track_area_use_counter]->SetMarkerStyle(24);
                    PM_jet_area_tracks[track_area_use_counter]->SetMarkerColor(jet_color[Jet_index]);
                    PM_jet_area_tracks[track_area_use_counter]->SetMarkerSize(0.6);
                }
                track_area_use_counter++;
            }
        }
        track_counter++;
    }
    //-------------------------------------------------------------------------------



    //-------------------------------------------------
    TCanvas* c_2D_jet_eta_vs_phi = new TCanvas("c_2D_jet_eta_vs_phi","c_2D_jet_eta_vs_phi",10, 10, 850*1.5, 350*1.5); // ww, wh
    c_2D_jet_eta_vs_phi->cd();
    c_2D_jet_eta_vs_phi->cd()->SetTicks(1,1);
    c_2D_jet_eta_vs_phi->cd()->SetGrid(0,0);
    c_2D_jet_eta_vs_phi->cd()->SetRightMargin(0.1);
    c_2D_jet_eta_vs_phi->cd()->SetLeftMargin(0.22);
    c_2D_jet_eta_vs_phi->cd()->SetBottomMargin(0.2);
    c_2D_jet_eta_vs_phi->cd()->SetLogy(0);

    h_2D_jet_eta_vs_phi->SetStats(0);
    h_2D_jet_eta_vs_phi->SetTitle("");
    h_2D_jet_eta_vs_phi->GetXaxis()->SetTitleOffset(1.2);
    h_2D_jet_eta_vs_phi->GetYaxis()->SetTitleOffset(1.2);
    h_2D_jet_eta_vs_phi->GetZaxis()->SetTitleOffset(0.5);
    h_2D_jet_eta_vs_phi->GetXaxis()->SetLabelSize(0.065);
    h_2D_jet_eta_vs_phi->GetYaxis()->SetLabelSize(0.065);
    h_2D_jet_eta_vs_phi->GetZaxis()->SetLabelSize(0.065);
    h_2D_jet_eta_vs_phi->GetXaxis()->SetTitleSize(0.065);
    h_2D_jet_eta_vs_phi->GetYaxis()->SetTitleSize(0.065);
    h_2D_jet_eta_vs_phi->GetZaxis()->SetTitleSize(0.065);
    h_2D_jet_eta_vs_phi->GetXaxis()->SetNdivisions(505,'N');
    h_2D_jet_eta_vs_phi->GetYaxis()->SetNdivisions(505,'N');
    h_2D_jet_eta_vs_phi->GetZaxis()->SetNdivisions(505,'N');
    h_2D_jet_eta_vs_phi->GetXaxis()->CenterTitle();
    h_2D_jet_eta_vs_phi->GetYaxis()->CenterTitle();
    h_2D_jet_eta_vs_phi->GetZaxis()->CenterTitle();
    h_2D_jet_eta_vs_phi->GetXaxis()->SetTitle("#phi (rad)");
    h_2D_jet_eta_vs_phi->GetYaxis()->SetTitle("#eta");
    h_2D_jet_eta_vs_phi->GetZaxis()->SetTitle("p_{T} (GeV/c)");
    h_2D_jet_eta_vs_phi->GetXaxis()->SetRangeUser(0.0-0.5,2.0*Pi+0.5);
    h_2D_jet_eta_vs_phi->GetYaxis()->SetRangeUser(-1.1,1.1);
    //h_2D_jet_eta_vs_phi->RebinX(2);
    //h_2D_jet_eta_vs_phi->RebinY(2);
    h_2D_jet_eta_vs_phi->DrawCopy("Lego2");
    //-------------------------------------------------



    //-------------------------------------------------
    TCanvas* c_Jupiter = new TCanvas("Jupiter","Jupiter",10,10, 760*2, 400*2);
    c_Jupiter->Divide(2, 2);

    h_2D_jet_eta_vs_phi_no_fill->SetStats(0);
    h_2D_jet_eta_vs_phi_no_fill->SetTitle("");
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetTitleOffset(1.0);
    h_2D_jet_eta_vs_phi_no_fill->GetXaxis()->SetLabelSize(0.065);
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetLabelSize(0.065);
    h_2D_jet_eta_vs_phi_no_fill->GetXaxis()->SetTitleSize(0.065);
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetTitleSize(0.065);
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetTitleOffset(0.5);
    h_2D_jet_eta_vs_phi_no_fill->GetXaxis()->SetNdivisions(505,'N');
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetNdivisions(505,'N');
    h_2D_jet_eta_vs_phi_no_fill->GetXaxis()->CenterTitle();
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->CenterTitle();
    h_2D_jet_eta_vs_phi_no_fill->GetXaxis()->SetTitle("#phi (rad)");
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetTitle("#eta");
    h_2D_jet_eta_vs_phi_no_fill->GetXaxis()->SetRangeUser(0.0-0.5,2.0*Pi+0.5);
    h_2D_jet_eta_vs_phi_no_fill->GetYaxis()->SetRangeUser(-1.1,1.1);

    // define different phi & eta acceptance cuts
    const Int_t lo(0), hi(1);
    Float_t phi_acc[4][2] = {{0.0,    1.0*Pi},  
                             {1.0*Pi, 2.0*Pi},  
                             {0.5*Pi, 1.5*Pi},  
                             {0.0,    2.0*Pi}}; 
    Float_t eta_acc[4][2] = {{-1.0,   0.0},
                             {0.0,    1.0},
                             {-0.5,   0.5}, 
                             {-1.0,   1.0}};

    // plot for different acceptances of eta & phi
    for (Int_t i_canvas = 1; i_canvas <= 4; i_canvas++)
    {
        c_Jupiter->cd(i_canvas);
        c_Jupiter->cd(i_canvas)->SetTicks(1,1);
        c_Jupiter->cd(i_canvas)->SetGrid(0,0);
        c_Jupiter->cd(i_canvas)->SetRightMargin(0.1);
        c_Jupiter->cd(i_canvas)->SetLeftMargin(0.22);
        c_Jupiter->cd(i_canvas)->SetBottomMargin(0.2);
        c_Jupiter->cd(i_canvas)->SetLogy(0);

        h_2D_jet_eta_vs_phi_no_fill->DrawCopy("colz");

        // Draw acceptance boundaries for tracks
        // c=color, w=width, s=style
        //       x1,     x2,     y1,    y2, c, w, s 
        PlotLine(   phi_acc[i_canvas][lo], phi_acc[i_canvas][hi], 
                    eta_acc[i_canvas][lo], eta_acc[i_canvas][lo], 1, 1, 2); 
        PlotLine(   phi_acc[i_canvas][lo], phi_acc[i_canvas][hi],  
                    eta_acc[i_canvas][hi], eta_acc[i_canvas][hi], 1, 1, 2); 
        PlotLine(   phi_acc[i_canvas][lo], phi_acc[i_canvas][lo],    
                    eta_acc[i_canvas][lo], eta_acc[i_canvas][hi], 1, 1, 2); 
        PlotLine(   phi_acc[i_canvas][hi], phi_acc[i_canvas][hi],    
                    eta_acc[i_canvas][lo], eta_acc[i_canvas][hi], 1, 1, 2); 

        // Draw acceptance boundaries for jet centroids
        // c=color, w=width, s=style
        //       x1,     x2,      y1,       y2,      c, w, s 
        Double_t R = jet_radius;
        PlotLine(   phi_acc[i_canvas][lo],      phi_acc[i_canvas][hi], 
                    eta_acc[i_canvas][lo]+R,    eta_acc[i_canvas][lo]+R, 2, 1, 2); 
        PlotLine(   phi_acc[i_canvas][lo],      phi_acc[i_canvas][hi],  
                    eta_acc[i_canvas][hi]-R,    eta_acc[i_canvas][hi]-R, 2, 1, 2); 
        PlotLine(   phi_acc[i_canvas][lo],      phi_acc[i_canvas][lo],  
                    eta_acc[i_canvas][lo]+R,    eta_acc[i_canvas][hi]-R, 2, 1, 2); 
        PlotLine(   phi_acc[i_canvas][hi],      phi_acc[i_canvas][hi], 
                    eta_acc[i_canvas][lo]+R,    eta_acc[i_canvas][hi]-R, 2, 1, 2); 

        TBox* Box_fiducial_cut_upper = new TBox();
        Box_fiducial_cut_upper->SetFillColor(kRed-8);
        Box_fiducial_cut_upper->SetFillStyle(3001);
        Box_fiducial_cut_upper->SetX1(phi_acc[i_canvas][lo]);
        Box_fiducial_cut_upper->SetX2(phi_acc[i_canvas][hi]);
        Box_fiducial_cut_upper->SetY1(eta_acc[i_canvas][hi]-R);
        Box_fiducial_cut_upper->SetY2(eta_acc[i_canvas][hi]);
        Box_fiducial_cut_upper->Draw("same");

        TBox* Box_fiducial_cut_lower = new TBox();
        Box_fiducial_cut_lower->SetFillColor(kRed-8);
        Box_fiducial_cut_lower->SetFillStyle(3001);
        Box_fiducial_cut_lower->SetX1(phi_acc[i_canvas][lo]);
        Box_fiducial_cut_lower->SetX2(phi_acc[i_canvas][hi]);
        Box_fiducial_cut_lower->SetY1(eta_acc[i_canvas][lo]+R);
        Box_fiducial_cut_lower->SetY2(eta_acc[i_canvas][lo]);
        Box_fiducial_cut_lower->Draw("same");


        //PM_jet_real_tracks_trig->Draw();

        for(Int_t i_track = 0; i_track < track_real_use_counter; i_track++)
        {
            PM_jet_real_tracks[i_track]->Draw();
        }
        for(Int_t i_track = 0; i_track < track_area_use_counter; i_track++)
        {
            PM_jet_area_tracks[i_track]->Draw();
        }
        for(Int_t i_jet = 0; i_jet < Jet_index+1; i_jet++)
        {
            PM_jet_center[i_jet] ->Draw();

            //cout << "i_jet: " << i_jet << ", jet pt: " << TS_jet_label[i_jet].Data() << endl;
            plotTopLegend((char*)TS_jet_label[i_jet].Data(),Jet_center_array[i_jet][0]-0.2,Jet_center_array[i_jet][1]-0.15,0.035,1,0.0,42,0,0);
        }
    }
    //-------------------------------------------------



    //-------------------------------------------------
    TCanvas* c_2D_circle_jet = new TCanvas("c_2D_circle_jet","c_2D_circle_jet",10,10,800,800);
    c_2D_circle_jet->cd();
    c_2D_circle_jet->cd()->SetTicks(1,1);
    c_2D_circle_jet->cd()->SetGrid(0,0);
    c_2D_circle_jet->cd()->SetRightMargin(0.1);
    c_2D_circle_jet->cd()->SetLeftMargin(0.22);
    c_2D_circle_jet->cd()->SetBottomMargin(0.2);
    c_2D_circle_jet->cd()->SetLogy(0);

    TH1F* h_frame_2D_circle_jet = c_2D_circle_jet->cd(1)->DrawFrame(-15,-15,15,15,"h_frame_2D_circle_jet");
    h_frame_2D_circle_jet->GetXaxis()->CenterTitle();
    h_frame_2D_circle_jet->GetYaxis()->CenterTitle();
    h_frame_2D_circle_jet->GetXaxis()->SetTitle("");
    h_frame_2D_circle_jet->GetYaxis()->SetTitle("");
    h_frame_2D_circle_jet->GetXaxis()->SetLabelSize(0.0);
    h_frame_2D_circle_jet->GetYaxis()->SetLabelSize(0.0);
    h_frame_2D_circle_jet->GetXaxis()->SetTitleSize(0.0);
    h_frame_2D_circle_jet->GetYaxis()->SetTitleSize(0.0);
    h_frame_2D_circle_jet->SetStats(0);
    h_frame_2D_circle_jet->SetTitle("");
    h_frame_2D_circle_jet->GetXaxis()->SetAxisColor(10);
    h_frame_2D_circle_jet->GetYaxis()->SetAxisColor(10);


    TEllipse TE_sigma_high(0.0,0.0,radius_2D_jet+3.0-0.01,radius_2D_jet+3.0-0.01,0,360); // x1, y1, r1, r2, phimin(0), phimax(360)
    TE_sigma_high.SetFillColor(kGray);
    TE_sigma_high.SetFillStyle(3001);
    TE_sigma_high.DrawClone();

    TEllipse TE_sigma_low(0.0,0.0,radius_2D_jet-2.0,radius_2D_jet-2.0,0,360); // x1, y1, r1, r2, phimin(0), phimax(360)
    TE_sigma_low.SetFillColor(10);
    TE_sigma_low.SetFillStyle(3001);
    TE_sigma_low.DrawClone();



    //--------------------------------------------
    // Find leading and sub-leading jets
    Int_t max_jet_pt_index[2] = {0,0};
    Double_t max_jet_pt[2] = {-999.0,-999.0};
    for(Int_t i_jet = 0; i_jet < Jet_index+1; i_jet++)
    {
        Double_t jet_pt = Array_Jet_info[i_jet][2]; // has now lenght of jet pt
        if(jet_pt > max_jet_pt[0])
        {
            max_jet_pt[1] = max_jet_pt[0];
            max_jet_pt_index[1] = max_jet_pt_index[0];
            max_jet_pt[0] = jet_pt;
            max_jet_pt_index[0] = i_jet;
            
        }
        else
        {
            if(jet_pt > max_jet_pt[1])
            {
                max_jet_pt[1] = jet_pt;
                max_jet_pt_index[1] = i_jet;
            }
        }
    }
    //--------------------------------------------

    Double_t jet_radius     = 0.3;
    Double_t jet_radius_phi = jet_radius*TMath::RadToDeg();

    //--------------------------------------------
    for(Int_t i_jet_leading = 0; i_jet_leading < 2; i_jet_leading++)
    {
        Double_t jet_phi = Array_Jet_info[max_jet_pt_index[i_jet_leading]][1]*TMath::RadToDeg();
        Double_t jet_pt  = Array_Jet_info[max_jet_pt_index[i_jet_leading]][2];

        //TEllipse TE_leading_jet_radius_phi(0.0,0.0,radius_2D_jet+jet_pt,radius_2D_jet+jet_pt,jet_phi-jet_radius_phi,jet_phi+jet_radius_phi); // x1, y1, r1, r2, phimin(0), phimax(360)
        TEllipse TE_leading_jet_radius_phi(0.0,0.0,radius_2D_jet+3.6,radius_2D_jet+3.6,jet_phi-jet_radius_phi,jet_phi+jet_radius_phi); // x1, y1, r1, r2, phimin(0), phimax(360)
        TE_leading_jet_radius_phi.SetFillColor(kAzure+2);
        TE_leading_jet_radius_phi.SetFillStyle(3001);
        TE_leading_jet_radius_phi.DrawClone();
    }
    //--------------------------------------------



    Draw_Circle_Detector_2D(radius_2D_jet,radius_2D_jet,40,1,kGray+1,1,1,0,0); // r_in, r_out, n_radii, n_delpha_phi, color, style, width, x, y
    Draw_Circle_Detector_2D(radius_2D_jet+3.0,radius_2D_jet+3.0,40,2,kGray,1,2,0,0); // r_in, r_out, n_radii, n_delpha_phi, color, style, width, x, y
    Draw_Circle_Detector_2D(radius_2D_jet-2.0,radius_2D_jet-2.0,40,2,kGray,1,2,0,0); // r_in, r_out, n_radii, n_delpha_phi, color, style, width, x, y


    for(Int_t i_jet = 0; i_jet < Jet_index+1; i_jet++)
    {
        Double_t x_jet_start = radius_2D_jet*TMath::Cos(Array_Jet_info[i_jet][1]);
        Double_t y_jet_start = radius_2D_jet*TMath::Sin(Array_Jet_info[i_jet][1]);

        TVector2 xy_jet_vec(x_jet_start,y_jet_start);
        TVector2 xy_jet_vec_dir = xy_jet_vec;
        xy_jet_vec_dir /= radius_2D_jet; // normalized direction vector
        xy_jet_vec_dir *= Array_Jet_info[i_jet][2]; // has now lenght of jet pt

        TVector2 xy_jet_vec_border[4];
        xy_jet_vec_border[0] = xy_jet_vec;
        xy_jet_vec_border[1] = xy_jet_vec;
        xy_jet_vec_border[2] = xy_jet_vec_dir;
        xy_jet_vec_border[3] = xy_jet_vec_dir;

        xy_jet_vec_border[0] = xy_jet_vec_border[0].Rotate(-jet_radius);
        xy_jet_vec_border[1] = xy_jet_vec_border[1].Rotate(jet_radius);
        xy_jet_vec_border[2] = xy_jet_vec_border[2].Rotate(-jet_radius);
        xy_jet_vec_border[3] = xy_jet_vec_border[3].Rotate(jet_radius);


#if 0
        cout << "i_jet: " << i_jet << ", jetpt: " << Array_Jet_info[i_jet][2] << endl;
#endif

        Double_t x_jet_stop = x_jet_start + xy_jet_vec_dir.Px();
        Double_t y_jet_stop = y_jet_start + xy_jet_vec_dir.Py();

        PlotLine(x_jet_start,x_jet_stop,y_jet_start,y_jet_stop,2,4,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle

#if 0
        if(i_jet == max_jet_pt_index[0] || i_jet == max_jet_pt_index[1])
        {
            PlotLine(xy_jet_vec_border[0].Px(),xy_jet_vec_border[0].Px()+xy_jet_vec_border[2].Px(),xy_jet_vec_border[0].Py(),xy_jet_vec_border[0].Py()+xy_jet_vec_border[2].Py(),kAzure+2,1,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
            PlotLine(xy_jet_vec_border[1].Px(),xy_jet_vec_border[1].Px()+xy_jet_vec_border[3].Px(),xy_jet_vec_border[1].Py(),xy_jet_vec_border[1].Py()+xy_jet_vec_border[3].Py(),kAzure+2,1,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
        }
#endif

    }
    for(Int_t i_track = 0; i_track < track_real_use_counter; i_track++)
    {
        Double_t x_jet_start = radius_2D_jet*TMath::Cos(Array_track_info[i_track][1]);
        Double_t y_jet_start = radius_2D_jet*TMath::Sin(Array_track_info[i_track][1]);

        TVector2 xy_jet_vec(x_jet_start,y_jet_start);
        TVector2 xy_jet_vec_dir = xy_jet_vec;
        xy_jet_vec_dir /= radius_2D_jet; // normalized direction vector
        xy_jet_vec_dir *= Array_track_info[i_track][2]; // has now lenght of jet pt


#if 0
        cout << "i_track: " << i_track << ", trackpt: " << Array_track_info[i_track][2] << endl;
#endif

        Double_t x_jet_stop = x_jet_start + xy_jet_vec_dir.Px();
        Double_t y_jet_stop = y_jet_start + xy_jet_vec_dir.Py();

        PlotLine(x_jet_start,x_jet_stop,y_jet_start,y_jet_stop,kGray+2,2,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
    }
    //-------------------------------------------------




}
