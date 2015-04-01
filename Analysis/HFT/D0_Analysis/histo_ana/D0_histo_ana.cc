
#include "functions.h"

void D0_histo_ana()
{
    cout << "D0_hist_ana started" << endl;

    SetRootGraphicStyle();

    TString HistName;

    TFile* infiles[2];
    infiles[0] = TFile::Open("./Data/all_same_event.root");
    infiles[1] = TFile::Open("./Data/all_mixed_event.root");

    cout << "Input files opened" << endl;

    // NOTE: ensure these values match D0_Analysis.cxx
    const Int_t N_h_InvMass    = 6;
    const Int_t N_h_InvMass_pt = 6;

    // create InvMass histograms:
    // --> first five indices are topological cut indices
    // --> then pT cut index
    // --> final index is 3 (same/mixed/both)
    TH1D* h_InvMass[N_h_InvMass][N_h_InvMass][N_h_InvMass][N_h_InvMass][N_h_InvMass][N_h_InvMass_pt][3];
    TH1D* h_max_sig[N_h_InvMass_pt];
    TH1D* h_max_mixed[N_h_InvMass_pt];
    TH1D* h_max_same[N_h_InvMass_pt];
    Double_t max_signif[N_h_InvMass_pt];
    for (Int_t i = 0; i < N_h_InvMass_pt; i++){ max_signif[i] = 0; }

    // fill all same/mixed hists with data and name them
	for(Int_t A = 0; A < N_h_InvMass; A++){
        for(Int_t B = 0; B < N_h_InvMass; B++){
        for(Int_t X = 0; X < N_h_InvMass; X++){
        for(Int_t Y = 0; Y < N_h_InvMass; Y++){
        for(Int_t AB = 0; AB < N_h_InvMass; AB++){
        for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){
        for(Int_t i_SE_ME = 0; i_SE_ME < 2; i_SE_ME++){
            HistName = "h_InvMass/h_InvMass_A";
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
            // assign each histogram to corresponding data from input files
            if (((TH1D*)infiles[i_SE_ME]->Get(HistName.Data()))->GetEntries())
			    h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME] = (TH1D*)infiles[i_SE_ME]->Get(HistName.Data());
            else {
                h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME] = NULL;
                continue;}
            
            HistName += "_";
            HistName += i_SE_ME;
            h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->SetName(HistName.Data());
            h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->Sumw2();
    }}}}}}}


    Double_t    Int_SE_ME[2]; 
    Bool_t      hist_exists = true;
    // normalize histograms and get subtracted result
	for(Int_t A = 0; A < N_h_InvMass; A++){
        for(Int_t B = 0; B < N_h_InvMass; B++){
        for(Int_t X = 0; X < N_h_InvMass; X++){
        for(Int_t Y = 0; Y < N_h_InvMass; Y++){
        for(Int_t AB = 0; AB < N_h_InvMass; AB++){
        for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){
        //std::cout << "; i_hist_pt = " << i_hist_pt << endl;
            /* ignoring this for right now
            HistName = "c_InvMass_cut";
            HistName += i_hist;
            HistName += "_pt";
            HistName += i_hist_pt;
            c_InvMass[i_hist][i_hist_pt] = new TCanvas(HistName.Data(),HistName.Data(),10,10,600,600);
            c_InvMass[i_hist][i_hist_pt]->cd()->SetFillColor(10);
            c_InvMass[i_hist][i_hist_pt]->cd()->SetTopMargin(0.1);
            c_InvMass[i_hist][i_hist_pt]->cd()->SetBottomMargin(0.2);
            c_InvMass[i_hist][i_hist_pt]->cd()->SetRightMargin(0.05);
            c_InvMass[i_hist][i_hist_pt]->cd()->SetLeftMargin(0.2);
            c_InvMass[i_hist][i_hist_pt]->cd()->SetLogy(0);
            */

            hist_exists = true;
            Int_SE_ME   = {0};
            // -------- normalization loop --------
            for(Int_t i_SE_ME = 0; i_SE_ME < 2; i_SE_ME++){
                
                // leave loop immediately if this hist has no data
                if (h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]==NULL){
                    hist_exists = false;
                    break;
                }
            
                // normalize regions outside D0 mass range
                Int_t start_int     = h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->FindBin(1.65);
                Int_t stop_int      = h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->FindBin(1.75);
                Int_SE_ME[i_SE_ME]  = h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->Integral(start_int,stop_int);

                start_int           = h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->FindBin(1.95);
                stop_int            = h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->FindBin(2.05);
                Int_SE_ME[i_SE_ME] += h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->Integral(start_int,stop_int);

                // mixed event
                if (i_SE_ME == 1)
                {

                    // normalize mixed event to same event
                    if(Int_SE_ME[i_SE_ME] > 0.0) h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->Scale(Int_SE_ME[0]/Int_SE_ME[i_SE_ME]);
                    HistName = h_InvMass[A][B][X][Y][AB][i_hist_pt][0] ->GetName();
                    HistName += "_sub";
                    h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME+1] = (TH1D*)h_InvMass[A][B][X][Y][AB][i_hist_pt][0]->Clone(HistName.Data());
                    h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME+1]->Add(h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME],-1.0);
                }
            } // -------- end normalization loop --------

            // restart loop if this hist has no data

            if (!hist_exists)
                continue;


            // calculate significance
            Double_t signal, background, this_signif;
            Int_t start = h_InvMass[A][B][X][Y][AB][i_hist_pt][1]->FindBin(1.85);
            Int_t stop  = h_InvMass[A][B][X][Y][AB][i_hist_pt][1]->FindBin(1.89);
            signal      = h_InvMass[A][B][X][Y][AB][i_hist_pt][2]->Integral(start, stop);
            background  = h_InvMass[A][B][X][Y][AB][i_hist_pt][1]->Integral(start, stop);
            if (signal + background > 0){
                this_signif = signal / TMath::Sqrt(signal + background);
                if (this_signif > max_signif[i_hist_pt]){
                    max_signif[i_hist_pt]   = this_signif;
                    h_max_sig[i_hist_pt]    = (TH1D*)h_InvMass[A][B][X][Y][AB][i_hist_pt][2]->Clone();
                    h_max_mixed[i_hist_pt]  = (TH1D*)h_InvMass[A][B][X][Y][AB][i_hist_pt][1]->Clone();
                    h_max_same[i_hist_pt]   = (TH1D*)h_InvMass[A][B][X][Y][AB][i_hist_pt][0]->Clone();
                    std::cout << "New max significance = " << max_signif[i_hist_pt];
                    std::cout << ", at i_hist_pt = " << i_hist_pt << endl;
                }
            }


    }}}}}}


    TCanvas * c_InvMass_sub[N_h_InvMass_pt];

    for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){
            HistName = "c_InvMass_sub_pt";
            HistName += i_hist_pt;
            c_InvMass_sub[i_hist_pt] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1000,1000);
            c_InvMass_sub[i_hist_pt]->Divide(2, 1);

            c_InvMass_sub[i_hist_pt]->cd(2)->SetFillColor(10);
            c_InvMass_sub[i_hist_pt]->cd(2)->SetTopMargin(0.1);
            c_InvMass_sub[i_hist_pt]->cd(2)->SetBottomMargin(0.2);
            c_InvMass_sub[i_hist_pt]->cd(2)->SetRightMargin(0.05);
            c_InvMass_sub[i_hist_pt]->cd(2)->SetLeftMargin(0.2);
            c_InvMass_sub[i_hist_pt]->cd(2)->SetGrid(0, 0);

            c_InvMass_sub[i_hist_pt]->cd(1)->SetFillColor(10);
            c_InvMass_sub[i_hist_pt]->cd(1)->SetTopMargin(0.1);
            c_InvMass_sub[i_hist_pt]->cd(1)->SetBottomMargin(0.2);
            c_InvMass_sub[i_hist_pt]->cd(1)->SetRightMargin(0.05);
            c_InvMass_sub[i_hist_pt]->cd(1)->SetLeftMargin(0.2);
            c_InvMass_sub[i_hist_pt]->cd(1)->SetLogy(0);
            c_InvMass_sub[i_hist_pt]->cd(1)->SetGrid(0, 0);

            h_max_sig[i_hist_pt]->SetStats(0);
            h_max_sig[i_hist_pt]->SetTitle("");
            h_max_sig[i_hist_pt]->GetYaxis()->SetTitleOffset(1.5);
            h_max_sig[i_hist_pt]->GetXaxis()->SetLabelSize(0.06);
            h_max_sig[i_hist_pt]->GetYaxis()->SetLabelSize(0.06);
            h_max_sig[i_hist_pt]->GetXaxis()->SetTitleSize(0.06);
            h_max_sig[i_hist_pt]->GetYaxis()->SetTitleSize(0.06);
            h_max_sig[i_hist_pt]->GetXaxis()->SetNdivisions(505,'N');
            h_max_sig[i_hist_pt]->GetYaxis()->SetNdivisions(505,'N');
            h_max_sig[i_hist_pt]->GetXaxis()->CenterTitle();
            h_max_sig[i_hist_pt]->GetYaxis()->CenterTitle();
            h_max_sig[i_hist_pt]->GetXaxis()->SetTitle("M(K,#pi) (GeV/c^{2})");
            h_max_sig[i_hist_pt]->GetYaxis()->SetTitle("counts");
            h_max_sig[i_hist_pt]->Sumw2();
            h_max_sig[i_hist_pt]->Rebin(1);
            h_max_sig[i_hist_pt]->GetXaxis()->SetRangeUser(1.64,2.06);
            h_max_sig[i_hist_pt]->SetLineColor(1);

            h_max_mixed[i_hist_pt]->SetStats(0);
            h_max_mixed[i_hist_pt]->SetTitle("");
            h_max_mixed[i_hist_pt]->GetYaxis()->SetTitleOffset(1.5);
            h_max_mixed[i_hist_pt]->GetXaxis()->SetLabelSize(0.06);
            h_max_mixed[i_hist_pt]->GetYaxis()->SetLabelSize(0.06);
            h_max_mixed[i_hist_pt]->GetXaxis()->SetTitleSize(0.06);
            h_max_mixed[i_hist_pt]->GetYaxis()->SetTitleSize(0.06);
            h_max_mixed[i_hist_pt]->GetXaxis()->SetNdivisions(505,'N');
            h_max_mixed[i_hist_pt]->GetYaxis()->SetNdivisions(505,'N');
            h_max_mixed[i_hist_pt]->GetXaxis()->CenterTitle();
            h_max_mixed[i_hist_pt]->GetYaxis()->CenterTitle();
            h_max_mixed[i_hist_pt]->GetXaxis()->SetTitle("M(K,#pi) (GeV/c^{2})");
            h_max_mixed[i_hist_pt]->GetYaxis()->SetTitle("counts");
            h_max_mixed[i_hist_pt]->Sumw2();
            h_max_mixed[i_hist_pt]->Rebin(1);
            h_max_mixed[i_hist_pt]->GetXaxis()->SetRangeUser(1.64,2.06);
            h_max_mixed[i_hist_pt]->SetLineColor(kBlue+2);
            h_max_mixed[i_hist_pt]->SetFillStyle(3001);
            h_max_mixed[i_hist_pt]->SetFillColor(kBlue+2);

            h_max_same[i_hist_pt]->SetStats(0);
            h_max_same[i_hist_pt]->SetTitle("");
            h_max_same[i_hist_pt]->GetYaxis()->SetTitleOffset(1.5);
            h_max_same[i_hist_pt]->GetXaxis()->SetLabelSize(0.06);
            h_max_same[i_hist_pt]->GetYaxis()->SetLabelSize(0.06);
            h_max_same[i_hist_pt]->GetXaxis()->SetTitleSize(0.06);
            h_max_same[i_hist_pt]->GetYaxis()->SetTitleSize(0.06);
            h_max_same[i_hist_pt]->GetXaxis()->SetNdivisions(505,'N');
            h_max_same[i_hist_pt]->GetYaxis()->SetNdivisions(505,'N');
            h_max_same[i_hist_pt]->GetXaxis()->CenterTitle();
            h_max_same[i_hist_pt]->GetYaxis()->CenterTitle();
            h_max_same[i_hist_pt]->GetXaxis()->SetTitle("M(K,#pi) (GeV/c^{2})");
            h_max_same[i_hist_pt]->GetYaxis()->SetTitle("counts");
            h_max_same[i_hist_pt]->Sumw2();
            h_max_same[i_hist_pt]->Rebin(1);
            h_max_same[i_hist_pt]->GetXaxis()->SetRangeUser(1.64,2.06);
            h_max_same[i_hist_pt]->SetLineColor(kRed+2);

           /* 
            h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->SetLineColor(kAzure+3);
                    h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->SetFillColor(kAzure+2);
                    h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->SetFillStyle(3001);
                    h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME]->DrawCopy("h same");
                    */

            TLegend * leg1 = new TLegend(0.75, 0.75, 0.95, 0.95);
            leg1->AddEntry(h_max_same[i_hist_pt], "Same Event", "l");
            leg1->AddEntry(h_max_mixed[i_hist_pt], "Mixed Event", "lf");

            TLegend * leg2 = new TLegend(0.75, 0.75, 0.95, 0.95);
            leg2->AddEntry(h_max_sig[i_hist_pt], "Signal (Same - Mixed)", "l");

            c_InvMass_sub[i_hist_pt]->cd(1);
            h_max_same[i_hist_pt]->DrawCopy("h");
            h_max_mixed[i_hist_pt]->DrawCopy("h same");
            leg1->Draw("same");

            c_InvMass_sub[i_hist_pt]->cd(2);
            h_max_sig[i_hist_pt]->DrawCopy("h");
            leg2->Draw("same");
            PlotLine(1.64,2.06,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
    }
}
