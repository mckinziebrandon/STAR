#include "../functions.h"

void plot_600(){

const Int_t N_h_InvMass_pt = 6;

    TFile * inFile = TFile::Open("signal_hists_600.root");

    TH1D * h_same   = (TH1D*)inFile->Get("h_max_same");
    TH1D * h_mixed  = (TH1D*)inFile->Get("h_max_mixed");
    TH1D * h_sig    = (TH1D*)inFile->Get("h_max_sig");

    TCanvas * c_InvMass_sub[N_h_InvMass_pt];
    for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){

        HistName = "c_InvMass_sub_pt";
        HistName += i_hist_pt;
        c_InvMass_sub[i_hist_pt] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1000,800);
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

        TLegend * leg1 = new TLegend(0.65, 0.8, 0.95, 0.95);
        leg1->SetFillColor(0);
        leg1->AddEntry(h_max_same[i_hist_pt], "Same Event", "l");
        leg1->AddEntry(h_max_mixed[i_hist_pt], "Mixed Event", "lf");

        TLegend * leg2 = new TLegend(0.65, 0.8, 0.95, 0.95);
        leg2->SetFillColor(0);
        leg2->AddEntry(h_max_sig[i_hist_pt], "Signal (Same - Mixed)", "l");

        c_InvMass_sub[i_hist_pt]->cd(1);
            h_max_same[i_hist_pt]->DrawCopy("h");
            h_max_mixed[i_hist_pt]->DrawCopy("h same");
            leg1->Draw("same");

        c_InvMass_sub[i_hist_pt]->cd(2);
            h_max_sig[i_hist_pt]->DrawCopy("h");
            leg2->Draw("same");
            sign_val->Draw("same");
            PlotLine(1.64,2.06,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle


    }





//h_InvMass[A][B][X][Y][AB][i_hist_pt][i_SE_ME] = (TH1D*)infiles[i_SE_ME]->Get(HistName.Data());
}
