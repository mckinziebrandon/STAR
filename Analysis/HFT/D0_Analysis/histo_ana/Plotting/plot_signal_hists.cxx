/*----------------------------------------------------------------------------

File Name:      plot_signal_hists.cxx
Author:         Brandon McKinzie
Description:    Plotting customization for optimized D0 signal histograms.
Input:          Fully optimized D0 signals histograms for each pT bin.
Output:         Pretty plots.

----------------------------------------------------------------------------*/


#include "../functions.h"
#include "TPaveLabel.h"
#include "TPaveText.h"

void plot_signal_hists(){

    const Int_t N_h_InvMass_pt = 6;
    const Int_t N_cuts = 5;

    TFile * inFile = TFile::Open("./Data/signal_hists_narrowpt3.root");
    TFile * outFile = new TFile("./Plots/final_hists_narrowpt3.root", "RECREATE");
    outFile->cd();

    TH1D * h_same[N_h_InvMass_pt];      
    TH1D * h_mixed[N_h_InvMass_pt];
    TH1D * h_sig[N_h_InvMass_pt];
    TH1D * h_sig_vals = (TH1D*)inFile->Get("h_sig_vals;1"); // each binContent = a significance value
    TH1D * h_cuts[N_h_InvMass_pt];                          // each binContent = a topological cut value
    TString HistName;

    TCanvas * c_InvMass_sub[N_h_InvMass_pt];
    //for(Int_t i_hist_pt = 0; i_hist_pt < N_h_InvMass_pt; i_hist_pt++){
    for(Int_t i_hist_pt = 3; i_hist_pt < 4; i_hist_pt++){

        // retrieve histograms from input file
        HistName = "sig_pt_";
            HistName += i_hist_pt;
            HistName += ";1";
            h_sig[i_hist_pt] = (TH1D*)inFile->Get(HistName.Data());
        HistName = "mixed_pt_";
            HistName += i_hist_pt;
            HistName += ";1";
            h_mixed[i_hist_pt] = (TH1D*)inFile->Get(HistName.Data());
        HistName = "same_pt_";
            HistName += i_hist_pt;
            HistName += ";1";
            h_same[i_hist_pt] = (TH1D*)inFile->Get(HistName.Data());
        HistName = "h_cuts_pt_";
            HistName += i_hist_pt;
            HistName += ";1";
            h_cuts[i_hist_pt] = (TH1D*)inFile->Get(HistName.Data());

            

        // canvas customization
        HistName = "c_InvMass_sub_pt";
            HistName += i_hist_pt;
            c_InvMass_sub[i_hist_pt] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1200,600);
            c_InvMass_sub[i_hist_pt]->Divide(3, 1);
        c_InvMass_sub[i_hist_pt]->cd(3)->SetFillColor(10);
            c_InvMass_sub[i_hist_pt]->cd(3)->SetTopMargin(0.1);
            c_InvMass_sub[i_hist_pt]->cd(3)->SetBottomMargin(0.2);
            c_InvMass_sub[i_hist_pt]->cd(3)->SetRightMargin(0.05);
            c_InvMass_sub[i_hist_pt]->cd(3)->SetLeftMargin(0.2);
            c_InvMass_sub[i_hist_pt]->cd(3)->SetGrid(0, 0);
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

        // obtain significance value
        Double_t sig = h_sig_vals->GetBinContent(i_hist_pt+1);
        TString  significance = "Significance: ";
                 significance += sig;
            
        // obtain topological cut values
        TString     cut_value;
        TString     cut_names[N_cuts] = {"dcaA", "dcaB", "VerdistX", "cos(#theta) #times 10^{4}", "dcaAB"};
        TPaveText*  cut_text = new TPaveText(0.05, 0.1, 0.95, 0.8);
        cut_text->AddText("Optimized Topological Cuts:");
        for (Int_t i_cut = 0; i_cut < N_cuts; i_cut++){
            cut_value  = cut_names[i_cut];
            cut_value += ": ";
            sprintf(NoP, "%3.0f", h_cuts[i_hist_pt]->GetBinContent(i_cut + 1));
            cut_value += NoP;
            cut_text->AddText(cut_value.Data());
        }
        
        // create plot legends
        TLegend * leg1 = new TLegend(0.62, 0.85, 0.98, 0.95);
            leg1->SetFillColor(0);
            leg1->AddEntry(h_same[i_hist_pt], "Same Event", "l");
            leg1->AddEntry(h_mixed[i_hist_pt], "Mixed Event", "lf");
        TLegend * leg2 = new TLegend(0.62, 0.85, 0.98, 0.95);
            leg2->SetFillColor(0);
            leg2->AddEntry(h_sig[i_hist_pt], "Signal (Same - Mixed)", "l");
            leg2->AddEntry((TObject*)0, significance.Data(), "");

        // draw on canvas
        c_InvMass_sub[i_hist_pt]->cd(1);
            h_same[i_hist_pt]->DrawCopy("h");
            h_mixed[i_hist_pt]->DrawCopy("hsame");
            leg1->Draw("same");
        c_InvMass_sub[i_hist_pt]->cd(2);
            h_sig[i_hist_pt]->DrawCopy("h");
            leg2->Draw("same");
            PlotLine(1.64,2.06,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
        c_InvMass_sub[i_hist_pt]->cd(3);
            cut_text->Draw();

        // save each canvas to output
        c_InvMass_sub[i_hist_pt]->Write();

    }

outFile->Close();
}
