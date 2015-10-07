
using namespace std;

#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
//#include <iostream.h>
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"
#include "TChain.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/GSLSimAnMinimizer.h"
//#include "Math/Functor.h"
//#include "TMinuitMinimizer.h"
//#include <GSLMultiMinimizer.h>
//#include "Math/GSLMultiMinimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "TRotation.h"
#include "TSVDUnfold.h"
#include "TSystemDirectory.h"
#include "TMatrixT.h"
#include "TVector2.h"
#include "TEllipse.h"
#include "TPaveLabel.h"
#include "TBox.h"
#include "TSystemDirectory.h"




//------------------------------------------------------------------------------------------------------------
static const Float_t Pi = TMath::Pi();
static TRandom ran;
static TString HistName, HistNameB;
static char NoP[50];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
static const Int_t N_v2_vs_pt_BW = 14;
static Double_t mt_m0_cut;
static Double_t pt_cut_BW_global;
static Int_t flag_v2_BW_use[N_v2_vs_pt_BW];
static TGraphAsymmErrors* tgae_v2_stat_BW[N_v2_vs_pt_BW];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
TVector3 Calc_Circle_Center(TVector3 v_pointA, TVector3 v_pointB, TVector3 v_pointC)
{
    TVector3 v_center_of_circle;

    // Calculate center of circle
    Double_t xy_circle[2] = {0.0,0.0};
    Double_t x12  = v_pointA.X() - v_pointB.X();
    Double_t x13  = v_pointA.X() - v_pointC.X();
    Double_t y12  = v_pointA.Y() - v_pointB.Y();
    Double_t y13  = v_pointA.Y() - v_pointC.Y();
    Double_t x12p = v_pointA.X() + v_pointB.X();
    Double_t x13p = v_pointA.X() + v_pointC.X();
    Double_t y12p = v_pointA.Y() + v_pointB.Y();
    Double_t y13p = v_pointA.Y() + v_pointC.Y();

    xy_circle[0] = (y12*(x13*x13p + y13*y13p) - y13*(x12*x12p + y12*y12p))/(2.0*(x13*y12 - x12*y13));
    xy_circle[1] = (x13*(x12*x12p + y12*y12p) - x12*(x13*x13p + y13*y13p))/(2.0*(x13*y12 - x12*y13));

    v_center_of_circle.SetXYZ(xy_circle[0],xy_circle[1],0.0);
    return v_center_of_circle;
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void Convert_NDC_to_User(TPad* iPad, Double_t x_NDC, Double_t y_NDC,
                         Double_t &x_User, Double_t &y_User,
                         Int_t flag_log_x, Int_t flag_log_y)
{
    // Convert NDC coordinates to User coordinates, if x/y axes are in log scale set flag_log to 1
    Double_t xr = iPad->GetX2()-iPad->GetX1();
    Double_t yr = iPad->GetY2()-iPad->GetY1();
    x_User = x_NDC*xr + iPad->GetX1();
    y_User = y_NDC*yr + iPad->GetY1();

    if(flag_log_x == 1) x_User = TMath::Power(10,x_User);
    if(flag_log_y == 1) y_User = TMath::Power(10,y_User);
    cout << "X1: " << iPad->GetX1() << ", X2: " << iPad->GetX2() << ", Y1: " <<
        iPad->GetY1() << ", Y2: " << iPad->GetY2() << ", xr: " << xr << ", yr: " << yr << endl;
}

void Convert_User_to_NDC(TPad* iPad, Double_t x_User, Double_t y_User, Double_t &x_NDC, Double_t &y_NDC)
{
    // Convert User coordinates to NDC coordinates
    Double_t xr = iPad->GetX2()-iPad->GetX1();
    Double_t yr = iPad->GetY2()-iPad->GetY1();
    x_NDC = (x_User-iPad->GetX1())/ xr;
    y_NDC = (y_User-iPad->GetY1())/ xr;
}
//------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color
    // align: 1 left aligned, 32, right aligned

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistLine2(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Float_t x_start, Float_t x_stop, Float_t y_min, Float_t y_max)
{
    TPolyLine* Hist_line = new TPolyLine();
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y >= y_min && y < y_max && y != 0 && x >= x_start && x <= x_stop)
        {
            //cout << "x = " << x << ", y = " << y << endl;
            Hist_line->SetNextPoint(x,y);
        }
    }
    Hist_line -> Draw();
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistErrorBand(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    // Modified March 14th, which year?
    Int_t N_points = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y != 0 && x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }
    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line       = new TGraph(N_total_points);
    TGraph* Hist_line_upper = new TGraph(N_points);
    TGraph* Hist_line_lower = new TGraph(N_points);
    Hist_line ->SetLineWidth(LineWidth);
    Hist_line ->SetLineStyle(LineStyle);
    Hist_line ->SetLineColor(Line_Col);
    Hist_line ->SetFillStyle(FillStyle);
    Hist_line ->SetFillColor(FillColor);

    Hist_line_upper ->SetLineWidth(LineWidth);
    Hist_line_upper ->SetLineStyle(LineStyle);
    Hist_line_upper ->SetLineColor(Line_Col);
    Hist_line_lower ->SetLineWidth(LineWidth);
    Hist_line_lower ->SetLineStyle(LineStyle);
    Hist_line_lower ->SetLineColor(Line_Col);

    Int_t N_point = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y     = Histo->GetBinContent(i);
        Double_t y_err = Histo->GetBinError(i);
        if(y != 0)
        {
            Double_t x = Histo->GetBinCenter(i);
            if(x >= x_start && x <= x_stop)
            {
                //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
                Hist_line->SetPoint(N_point,x,y-y_err);
                Hist_line_upper->SetPoint(N_point,x,y+y_err);
                Hist_line_lower->SetPoint(N_point,x,y-y_err);
                Hist_line->SetPoint(N_total_points-2-N_point,x,y+y_err);
                if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-y_err);
                N_point++;
            }
        }
    }
    Hist_line       -> Draw("f");
    Hist_line_upper -> Draw("l");
    Hist_line_lower -> Draw("l");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotGraphErrorBand(TGraphAsymmErrors* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    Int_t N_points = 0;
    for(Int_t i_point = 0; i_point < Histo->GetN(); i_point++)
    {
        Double_t x,y;
        Histo->GetPoint(i_point,x,y);

        if(x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }

    Int_t N_total_points = N_points*2+1;
    cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line       = new TGraph(N_total_points);
    TGraph* Hist_line_upper = new TGraph(N_points);
    TGraph* Hist_line_lower = new TGraph(N_points);
    Hist_line ->SetLineWidth(LineWidth);
    Hist_line ->SetLineStyle(LineStyle);
    Hist_line ->SetLineColor(Line_Col);
    Hist_line ->SetFillStyle(FillStyle);
    Hist_line ->SetFillColor(FillColor);

    Hist_line_upper ->SetLineWidth(LineWidth);
    Hist_line_upper ->SetLineStyle(LineStyle);
    Hist_line_upper ->SetLineColor(Line_Col);
    Hist_line_lower ->SetLineWidth(LineWidth);
    Hist_line_lower ->SetLineStyle(LineStyle);
    Hist_line_lower ->SetLineColor(Line_Col);

    Int_t N_point = 0;
    for(Int_t i_point = 0; i_point < Histo->GetN(); i_point++)
    {
        Double_t x,y;
        Histo->GetPoint(i_point,x,y);
        Double_t Xhigh = Histo->GetErrorXhigh(i_point);
        Double_t Xlow  = Histo->GetErrorXlow(i_point);
        Double_t Yhigh = Histo->GetErrorYhigh(i_point);
        Double_t Ylow  = Histo->GetErrorYlow(i_point);

        if(x >= x_start && x <= x_stop)
        {
            //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
            Hist_line->SetPoint(N_point,x,y-Ylow);
            Hist_line_upper->SetPoint(N_point,x,y+Yhigh);
            Hist_line_lower->SetPoint(N_point,x,y-Ylow);
            Hist_line->SetPoint(N_total_points-2-N_point,x,y+Yhigh);
            if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-Ylow);
            N_point++;
        }
    }
    Hist_line       -> Draw("f");
    Hist_line_upper -> Draw("l");
    Hist_line_lower -> Draw("l");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t One_over_x_FitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    y = 0.0;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    //if(x != 0.0 && !((par0/x) < 0.0 && par2 < 1.0)) y = pow(par0/x,par2) + par1;
    if(x != 0.0) y = par0*TMath::Power(1.0/x,par2) + par1;
    //if(x != 0.0) y =par0/x + par1 + 0.0001*par2;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t TwoGaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + par3*TMath::Gaus(x,par4,par5,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussPolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, pol0, pol1, pol2, pol3, pol4, pol5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    pol0  = par[3];
    pol1  = par[4];
    pol2  = par[5];
    pol3  = par[6];
    pol4  = par[7];
    pol5  = par[8];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussAreaPolyFitFunc(Double_t* x_val, Double_t* par)
{
    // From http://www.originlab.com/doc/Origin-Help/Gaussian-Function-FitFunc
    // Is using the area instead of amplitude
    // Area: area
    // xc: center
    // w: FWHM
    Double_t x, y, Area, xc, w, pol0, pol1, pol2, pol3, pol4, pol5;
    Area  = par[0];
    xc    = par[1];
    w     = par[2];
    pol0  = par[3];
    pol1  = par[4];
    pol2  = par[5];
    pol3  = par[6];
    pol4  = par[7];
    pol5  = par[8];
    x = x_val[0];
    y = Area*TMath::Exp(-4.0*TMath::Log(2)*TMath::Power(x-xc,2)/TMath::Power(w,2))/(w*TMath::Sqrt(TMath::Pi()/(4.0*TMath::Log(2)))) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t FlowFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t phi, y, v0, v1, v2, v3, v4;
    v0  = par[0];
    v1  = par[1];
    v2  = par[2];
    v3  = par[3];
    v4  = par[4];
    phi = x_val[0];
    y = v0 * (1.0 + 2.0*v1*TMath::Cos(phi) + 2.0*v2*TMath::Cos(2.0*phi)
              + 2.0*v3*TMath::Cos(3.0*phi) + 2.0*v4*TMath::Cos(4.0*phi));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod_x(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Double_t LevyFitFunc_pT(Double_t* x_val, Double_t* par)
{
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2 vs. pT
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT, a, b, c, d, n;
    pT = x_val[0];
    n  = par[0]; // number-of-constituent quarks
    a  = par[1];
    b  = par[2];
    c  = par[3];
    d  = par[4];

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. pT/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT_ncq, a, b, c, d, n;
    pT_ncq = x_val[0];
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT_ncq - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc_Poly(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0, par6, par7;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass
    par6   = par[6];
    par7   = par[7];

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = (a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d) + par6*mT_ncq + par7*mT_ncq*mT_ncq;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_m0, a, b, c, d, n, m0;
    mT_m0  = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_m0 + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncE(Double_t* x_val, Double_t* par)
{
    // Original
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        //for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            //r = r_start + j*delta_r;
            r = 0.01;

            //r = R; //

            rho = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_radial(Double_t* x_val, Double_t* par)
{
    // Blast wave function with radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            r = r_start + j*delta_r;

            //r = R; //

            rho  = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_no_mass_array(Double_t* x_val, Double_t* par)
{
    // Removed mass array
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Double_t Mass = par[5];
    Double_t mt = TMath::Sqrt(pt*pt + Mass*Mass);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 +=
            delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0
                                                                                         + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 +=
            delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 +
                                                                    2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] =
    {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 +=
            delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0
                                                                                         + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 +=
            delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 +
                                                                    2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



/*
//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14]  = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt        = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi    = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;

    T          = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------
*/


//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncC(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 100;
    Int_t nbins_r      = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}
    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncD(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t NCQ[14]  = {2.0,3.0,3.0,3.0,3.0,3.0,3.0,2.0,3.0,3.0,2.0,2.0,2.0,2.0};
    pt /= NCQ[PID];
    Mass[PID] /= NCQ[PID];
    //pt *= NCQ[PID];
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 25;
    Int_t nbins_r      = 25;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}

    v2 *= NCQ[PID];
    //v2 /= NCQ[PID];
    //v2 *= (1.0/(1.0 + TMath::Exp((pt-4.0/NCQ[PID])/0.6)));

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void BlastWaveSimultaneous(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
    Int_t npfits      = 0;
    Double_t chi2     = 0.0;
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};

    for(Int_t i = 0; i < 14; i++) // loop over PIDs
    {
        p[5] = i; // PID
        //Double_t pt_cut = TMath::Sqrt((Mass[i]+mt_m0_cut)*(Mass[i]+mt_m0_cut) - Mass[i]*Mass[i]); // mt-m0 cut
        Double_t pt_cut = pt_cut_BW_global;
        if(flag_v2_BW_use[i] == 1)
        {
            for(Int_t ix = 0; ix < tgae_v2_stat_BW[i]->GetN(); ix++)
            {
                Double_t x[] = {tgae_v2_stat_BW[i]->GetX()[ix]};
                if(x[0] < pt_cut)
                {
                    Double_t y      = tgae_v2_stat_BW[i]->GetY()[ix];
                    Double_t ye     = tgae_v2_stat_BW[i]->GetErrorYhigh(ix);
                    //ye = 0.01;
                    //cout << "ix = " << ix << ", x = " << x[0] << ", y = " << y << ", ye = " << ye << endl;
                    Double_t bw_val = BlastWaveFitFunc(x,p);
                    //Double_t bw_val = 0.1;
                    //ye += 0.003;
                    Double_t diff   = (y - bw_val)/ye;
                    chi2 += diff*diff;
                    npfits++;
                }
                else break;
            }
        }
    }

    fval = chi2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Hist_interpolate_and_error(TH1F* hist, Double_t x, Double_t &Int_val, Double_t &Int_err)
{
    // Linear interpolation of a histogram
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[2]   = {0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetEntries() == 0) // no data poins to interpolate
    {
        return_val = -1;
        Int_val = 0;
        Int_err = 0;
        //cout << "No entries in histogram" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t binx = 1; binx < hist->GetNbinsX(); binx++)
        {
            Double_t bin_error_val   = hist->GetBinError(binx);
            Double_t bin_x_pos       = hist->GetBinCenter(binx);
            if(bin_error_val != 0)
            {
                err_counter++;
                bin_entries[1] = hist->GetBinContent(binx);
                bin_error[1]   = hist->GetBinError(binx);
                bin_x_val[1]   = hist->GetBinCenter(binx);
                if(bin_x_pos >= x)
                {
                    binx_high = binx;
                    flag_max  = 1;
                    break;
                }
                else flag_max = 0;
            }
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val = bin_entries[1];
            Int_err = bin_error[1];
            return return_val;
        }
        for(Int_t binx_low = binx_high; binx_low > 0; binx_low--)
        {
            bin_entries[0] = hist->GetBinContent(binx_low);
            bin_error[0]   = hist->GetBinError(binx_low);
            bin_x_val[0]   = hist->GetBinCenter(binx_low);
            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        if(bin_error[0] != 0 && bin_error[1] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;
                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[1];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;
                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val = bin_entries[0];
                Int_err = bin_error[0];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t TGraphAsymmErrors_interpolate_and_error(TGraphAsymmErrors* hist, Double_t x, Double_t &Int_val, Double_t &Int_err_low,Double_t &Int_err_high)
{
    // V2: 03.12.2012 -> bug fixed to calculate error bars
    // Linear interpolation of a TGraphAsymmErrors
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[4]   = {0,0,0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetN() == 0) // no data poins to interpolate
    {
        return_val   = -1;
        Int_val      = 0;
        Int_err_low  = 0;
        Int_err_high = 0;
        //cout << "No entries in TGraphAsymmErrors" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t epoint = 0; epoint < hist->GetN(); epoint++)
        {
            hist->GetPoint(epoint,bin_x_val[1],bin_entries[1]);
            bin_error[2] = hist->GetErrorYlow(epoint);
            bin_error[3] = hist->GetErrorYhigh(epoint);

            err_counter++;
            if(bin_x_val[1] >= x)
            {
                binx_high = epoint;
                flag_max  = 1;
                break;
            }
            else flag_max = 0;
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val      = bin_entries[1];
            Int_err_low  = bin_error[2];
            Int_err_high = bin_error[3];
            return return_val;
        }
        for(Int_t epoint = binx_high; epoint >= 0; epoint--)
        {
            hist->GetPoint(epoint,bin_x_val[0],bin_entries[0]);
            bin_error[0] = hist->GetErrorYlow(epoint);
            bin_error[1] = hist->GetErrorYhigh(epoint);

            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        //cout << "bin0 = " << bin_error[0] << ", bin2 = " << bin_error[2] << endl;

        if(bin_error[0] != 0 && bin_error[2] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;

                //cout << "x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 << ", y2 = " << y2 << endl;

                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[2];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_low = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);

                y1_err = bin_error[1];
                y2_err = bin_error[3];

                termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                termC  = -(x-x2)/(x2-x1);
                termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_high = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val      = bin_entries[0];
                Int_err_low  = bin_error[0];
                Int_err_high = bin_error[1];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Get_T_muB_from_SHM(Double_t sqrt_sNN, Double_t &T, Double_t &muB)
{
    // Calculates baryon chemical potential muB and chemical freeze-out temperature T
    // based on fits to statistical hardronization model (SHM) calculations
    // from here: http://arxiv.org/pdf/hep-ph/0511094.pdf, Phys.Rev. C73 (2006) 034905
    const Double_t a = 0.166; // +/- 0.002
    const Double_t b = 0.139; // +/- 0.016
    const Double_t c = 0.053; // +/- 0.021
    const Double_t d = 1.308; // +/- 0.028
    const Double_t e = 0.273; // +/- 0.008

    if(sqrt_sNN > 0.0)
    {
        muB = d/(1.0 + e*sqrt_sNN);
        T   = a - b*muB*muB - c*muB*muB*muB*muB;
    }
    else
    {
        muB = -999.0;
        T   = -999.0;
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Get_muB(Double_t sqrt_sNN, Int_t flag_input, Double_t &muB_err)
{
    // Returns baryon chemical potential as a function of the center-of-mass energy sqrt(sNN)
    // from http://arxiv.org/pdf/1111.2406v1.pdf
    // flag_input: 0 = parametrization, 1 = STAR fits, 2 = corrected parametrization for 0-80% (see below)

    // corrected parametrization: Centrality data taken from: http://arxiv.org/pdf/0808.2041.pdf
    Double_t centrality_vals[9] = {75,65,55,45,35,25,15,7.5,2.5};
    Double_t muB_Energy[2][2][9] = // [62.4,200][values,error][70-80%,...,0-5%]
    {
        {   // Au+Au @ 62.4 GeV
            {37.7,42.5,47.0,51.3,54.2,54.5,59.4,61.0,62.7}, // values
            {6.5,5.8,5.1,5.2,5.2,5.2,5.4,5.7,6.0}  // errors
        },
        {   // Au+Au @ 200 GeV
            {14.1,15.3,17.7,18.9,18.6,21.3,21.0,22.8,21.9}, // values
            {4.2,4.2,4.2,4.2,4.2,4.2,4.2,4.5,4.5}  // errors
        }
    };

    Double_t muB     = 0.0;
    muB_err          = 0.0;
    Double_t a       = 1.482;  // +/- 0.0037 GeV
    Double_t b       = 0.3517; // +/- 0.009 GeV^-1

    if(!(flag_input == 0 || flag_input == 1 || flag_input == 2))
    {
        cout << "WARNING: Get_muB function, flag_input is wrong. Either 0 or 1. Parametrization is used." << endl;
        flag_input = 0;
    }

    if(sqrt_sNN >= 0.0 && (flag_input == 0 || flag_input == 2))
    {
        muB = a/(1.0+b*sqrt_sNN);
    }
    if(flag_input == 2) // correct the muB values for 0-80%
    {
        Double_t scale_fac[2][2];  // [62.4,200][val,error]
        Double_t muB_mean[2][2]; // [62.4,200][val,error]
        for(Int_t ebeam = 0; ebeam < 2; ebeam++)
        {
            for(Int_t val_err = 0; val_err < 2; val_err++) // calculate mean muB + error
            {
                muB_mean[ebeam][val_err] = 0.0; // approximation for 0-80%
                Double_t total_weight = 0.0;
                for(Int_t cent = 0; cent < 9; cent++)
                {
                    Double_t weight = 1.0;
                    if(cent >= 7) weight = 0.5;
                    if(val_err == 0) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent];
                    if(val_err == 1) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent]*weight*muB_Energy[ebeam][val_err][cent];
                    total_weight += weight;
                }
                if(total_weight > 0.0)
                {
                    if(val_err == 0) muB_mean[ebeam][val_err] /= total_weight;
                    if(val_err == 1)
                    {
                        muB_mean[ebeam][val_err] = TMath::Sqrt(muB_mean[ebeam][val_err]);
                        muB_mean[ebeam][val_err] /= total_weight;
                    }
                }
            }
            scale_fac[ebeam][0] = muB_mean[ebeam][0]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0); // <muB>/0-10% = 0-80%/0-10%
            Double_t term1 = muB_mean[ebeam][1]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0);
            Double_t term2 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][8]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            Double_t term3 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][7]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            scale_fac[ebeam][1] = TMath::Sqrt(term1*term1+term2*term2+term3*term3); // error on scaling factor
        }
        Double_t mean_scale_fac[2] = {(scale_fac[0][0]+scale_fac[1][0])/2.0,TMath::Sqrt(scale_fac[0][1]*scale_fac[0][1]+scale_fac[1][1]*scale_fac[1][1])/2.0}; // [value,error]
        Double_t original_muB = muB;
        muB     = muB * mean_scale_fac[0]; // 0-80% estimation
        muB_err = muB * mean_scale_fac[1]; // 0-80% error estimation
        //cout << "muB = " << original_muB << ", muB(0-80%) = " << muB << ", muB_err(0-80%) = " << muB_err << ", <scale> = " << mean_scale_fac[0] << ", <scale_err> = " << mean_scale_fac[1] << endl;
        //cout << "scale_fac(62 GeV) = " << scale_fac[0][0] << ", scale_fac(200 GeV) = " << scale_fac[1][0] << endl;
        //cout << "scale_fac_err(62 GeV) = " << scale_fac[0][1] << ", scale_fac_err(200 GeV) = " << scale_fac[1][1] << endl;
    }

    if(flag_input == 1)
    {
        // Preliminary muB parameters from Lokesh with strange particles for 0-80% Au+Au collisions
        // 7.7 :  3.61377e-01   1.51360e-02
        // 11.5:  2.59978e-01   1.06312e-02
        // 39  :  8.81846e-02   4.73062e-03
        // 200 :  22 -> Thats a guess

        // Fittet in TmuB_fit.cc macro with One_over_x_FitFunc function
        Double_t par0 = 2.13495e+00;
        Double_t par1 = 4.35916e-04;
        Double_t par2 = 8.67894e-01;
        muB = par0*TMath::Power(1.0/sqrt_sNN,par2) + par1;
    }

    return muB;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    //gStyle->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_new_Symbol(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_hist_new_Symbol(TH1D* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TH1D* ge_clone_A = (TH1D*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TH1D* ge_clone_B = (TH1D*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TH1D* ge_clone_C = (TH1D*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_Point_new_Symbol(Double_t x_val, Double_t y_val, Double_t x_min_err, Double_t x_max_err,
                                Double_t y_min_err, Double_t y_max_err,
                                Int_t style, Int_t color, Float_t size)
{
    TGraphAsymmErrors* tgae = new TGraphAsymmErrors();
    tgae->SetPoint(0,x_val,y_val);
    tgae->SetPointError(0,x_min_err,x_max_err,y_min_err,y_max_err);

    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_header(FILE* HTML_file)
{
    fprintf(HTML_file,"%s \n","<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">");
    fprintf(HTML_file,"%s \n","<html>");
    fprintf(HTML_file,"%s \n","<head>");


    fprintf(HTML_file,"%s \n","  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">");
    fprintf(HTML_file,"%s \n","  <title>STAR Physics Database</title>");
    fprintf(HTML_file,"%s \n","</head>");
    fprintf(HTML_file,"%s \n"," <body bgcolor=\"white\">");

    fprintf(HTML_file,"%s \n","<h2><span class=\"normaltext\"><font color=\"#ff8c00\" face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\">");
    fprintf(HTML_file,"%s \n","Observation of an energy-dependent difference in elliptic flow between particles and anti-particles in relativistic heavy ion collisions");
    fprintf(HTML_file,"%s \n","</font>");

    fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\"><br> <br> </font></span></h2>");

    //fprintf(HTML_file,"%s \n","<hr><h3><font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">Figure 1");
    //fprintf(HTML_file,"%s \n","</font></h3>");
    //fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">");
    //fprintf(HTML_file,"%s \n","<A HREF=\"fig1.png\"> png </A> | <A HREF=\"fig1.eps\"> eps </A> </br>");
    fprintf(HTML_file,"%s \n","</font>");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_title(TString Label, FILE* HTML_file)
{
    TString label_html = "<hr><h3><font face=\\\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\\\" color=\\\"#0000cd\\\">";
    label_html += Label.Data();
    label_html += " </font> </h3> </hr>";
    //label_html += "\"";
    fprintf(HTML_file,"%s \n",label_html.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_tables(TString labelX, TString labelY, TGraphAsymmErrors* data_stat, TGraphAsymmErrors* syst_errors, TGraphAsymmErrors* glob_syst_errors, TString Label, FILE* HTML_file)
{
    fprintf(HTML_file,"%s %s %s \n","<tr><td>",Label.Data(),"</td> </tr>");
    TString label_html = "<table border=1> <tr><th> ";
    label_html += labelX.Data();
    label_html += " </th> <th> ";
    label_html += labelY.Data();
    label_html += " </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>";

    //fprintf(HTML_file,"%s \n","<table border=1> <tr><th> p<sub>T</sub> (GeV/c) </th> <th> v<sub>2</sub> </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>");

    fprintf(HTML_file,"%s \n",label_html.Data());
    Int_t flag_glob_syst = 1;
    if(glob_syst_errors == syst_errors) flag_glob_syst = -1;

    for(Int_t epoint = 0; epoint < data_stat->GetN(); epoint++)
    {
        Double_t x_val, y_val, x_syst, stat_error, low_syst, high_syst, glob_syst;
        data_stat                     ->GetPoint(epoint,x_val,y_val);
        stat_error = data_stat        ->GetErrorYhigh(epoint);
        high_syst  = syst_errors      ->GetErrorYhigh(epoint);
        low_syst   = syst_errors      ->GetErrorYlow(epoint);
        glob_syst  = glob_syst_errors ->GetErrorYhigh(0);
        if(flag_glob_syst == -1) glob_syst = 0.0;

        fprintf(HTML_file,"%s %f %s %f %s %f %s %f %s %f %s %f %s \n","<tr> <td> ",x_val," </td> <td> ",y_val," </td> <td> ",stat_error," </td> <td> ",
                low_syst," </td> <td> ",high_syst," </td> <td> ",glob_syst,"</td></tr>");

    }
    fprintf(HTML_file,"%s \n","</table><hr><table border=1>");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
									 TVector3 &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.Mag()>0))
    {
      return -1000000.;
    }
  
  TVector3 diff = base-point;

  TVector3 cross = dir.Cross(diff);
  
  return cross.Mag()/dir.Mag();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector3 calculateDCA_vec_StraightToPoint(TVector3 &base, TVector3 &dir, TVector3 &point)
{
  // calculates the minimum distance vector of a point to a straight given as parametric straight x = base + n * dir

    TVector3 diff = base-point;
    TVector3 dir_norm = dir;
    dir_norm *= (1.0/dir.Mag());
    Double_t proj_val = diff.Dot(dir_norm);
    TVector3 proj_dir = dir_norm;
    proj_dir *= proj_val;

    TVector3 dist_vec = proj_dir - diff;

    return dist_vec;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t get_HFT_det_index(Int_t sensor_id)
{
    // Determines detector index (see definition below) based on the sonsor id
    // 0 = inner pixel left
    // 1 = outer pixel left
    // 2 = IST left
    // 3 = inner pixel right
    // 4 = outer pixel right
    // 5 = IST right

    Int_t HFT_det_index = -1;
    if(
       (sensor_id >= 1   && sensor_id <= 10)  ||
       (sensor_id >= 41  && sensor_id <= 50)  ||
       (sensor_id >= 81  && sensor_id <= 90)  ||
       (sensor_id >= 121 && sensor_id <= 130) ||
       (sensor_id >= 161 && sensor_id <= 170)
      )
    {
        HFT_det_index = 0;
    }

    if(
       (sensor_id >= 11   && sensor_id <= 40)  ||
       (sensor_id >= 51   && sensor_id <= 80)  ||
       (sensor_id >= 91   && sensor_id <= 120) ||
       (sensor_id >= 131  && sensor_id <= 160) ||
       (sensor_id >= 171  && sensor_id <= 200)
      )
    {
        HFT_det_index = 1;
    }
    if(
       sensor_id >= 1001 && sensor_id <= 1072
      )
    {
        HFT_det_index = 2;
    }
    if(
       (sensor_id >= 1+200   && sensor_id <= 10+200)  ||
       (sensor_id >= 41+200  && sensor_id <= 50+200)  ||
       (sensor_id >= 81+200  && sensor_id <= 90+200)  ||
       (sensor_id >= 121+200 && sensor_id <= 130+200) ||
       (sensor_id >= 161+200 && sensor_id <= 170+200)
      )
    {
        HFT_det_index = 3;
    }
    if(
       (sensor_id >= 11+200   && sensor_id <= 40+200)  ||
       (sensor_id >= 51+200   && sensor_id <= 80+200)  ||
       (sensor_id >= 91+200   && sensor_id <= 120+200) ||
       (sensor_id >= 131+200  && sensor_id <= 160+200) ||
       (sensor_id >= 171+200  && sensor_id <= 200+200)
      )
    {
        HFT_det_index = 4;
    }
    if(
       sensor_id >= 1073 && sensor_id <= 1145
      )
    {
        HFT_det_index = 5;
    }
    return HFT_det_index;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_histogram(TH1D* hist, TString plot_option, TString label_x, TString label_y, Double_t x_start, Double_t x_stop,
                    Double_t y_start, Double_t y_stop, Double_t title_offset_x, Double_t title_offset_y, Double_t label_offset_x,
                    Double_t label_offset_y, Double_t label_size, Double_t title_size, Int_t N_div)
{
    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(title_offset_x);
    hist->GetYaxis()->SetTitleOffset(title_offset_y);
    hist->GetXaxis()->SetLabelOffset(label_offset_x);
    hist->GetYaxis()->SetLabelOffset(label_offset_y);
    hist->GetXaxis()->SetLabelSize(label_size);
    hist->GetYaxis()->SetLabelSize(label_size);
    hist->GetXaxis()->SetTitleSize(title_size);
    hist->GetYaxis()->SetTitleSize(title_size);
    hist->GetXaxis()->SetNdivisions(N_div,'N');
    hist->GetYaxis()->SetNdivisions(N_div,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitle(label_x.Data());
    hist->GetYaxis()->SetTitle(label_y.Data());
    if(!(x_start == 0 && x_stop == 0)) hist->GetXaxis()->SetRangeUser(x_start,x_stop);
    if(!(y_start == 0 && y_stop == 0)) hist->GetYaxis()->SetRangeUser(y_start,y_stop);
    hist->DrawCopy(plot_option.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t get_HFT_pixel_sector(Int_t sensor_id, Int_t &IST_ladder)
{
    // Determines pixel sector based on the sonsor id

    // PXL
    Int_t HFT_pxl_sector = -1;
    if(sensor_id >= 1 && sensor_id <= 400)
    {
        HFT_pxl_sector = (Int_t)((sensor_id-1)/40);
    }

    // IST -> 24 ladders, 6 sensors per ladder -> 144 sensors in total, starting from sensor_id = 1001
    if(sensor_id >= 1001 && sensor_id <= 1144)
    {
        IST_ladder = (Int_t)((sensor_id-1001)/6);
    }
    else IST_ladder = -1;

    return HFT_pxl_sector;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_STAR_3D()
{

    TGeoManager *geom = new TGeoManager("geom","My 3D Project");
    //------------------Creat materials------------------------------
    TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
    TGeoMaterial *Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
    Fe->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_outer_tube = new TGeoMaterial("M_outer_tube",55.84,26.7,7.87);
    M_outer_tube->SetTransparency(93); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_IDS = new TGeoMaterial("M_IDS",55.84,26.7,7.87);
    M_IDS       ->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_beampipe = new TGeoMaterial("M_beampipe",55.84,26.7,7.87);
    M_beampipe       ->SetTransparency(70); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_Pixel_support = new TGeoMaterial("M_Pixel_support",55.84,26.7,7.87);
    M_Pixel_support    ->SetTransparency(70); // higher value means more transparent, 100 is maximum


    //------------------Create media---------------------------------
    TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);
    TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);
    TGeoMedium *Me_outer_tube = new TGeoMedium("Me_outer_tube",1,M_outer_tube);
    TGeoMedium *Me_IDS        = new TGeoMedium("Me_IDS",1,M_IDS);
    TGeoMedium *Me_beampipe   = new TGeoMedium("Me_beampipe",1,M_beampipe);
    TGeoMedium *Me_Pixel_support   = new TGeoMedium("Me_Pixel_support",1,M_Pixel_support);

    //------------------Create TOP volume----------------------------
    TGeoVolume *top = geom->MakeBox("top",Air,500,500,500);
    geom->SetTopVolume(top);
    geom->SetTopVisible(0);
    // If you want to see the boundary, please input the number, 1 instead of 0.
    // Like this, geom->SetTopVisible(1);


    TGeoVolume *inner_field_tube       = geom->MakeTube("inner_field_tube",Iron,49.5,50.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *outer_field_tube       = geom->MakeTube("outer_field_tube",Me_outer_tube,199.5,200.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_central_part       = geom->MakeTube("IDS_central_part",Me_IDS,42.8/2.0,43.0/2.0,56.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_side_parts         = geom->MakeTube("IDS_side_parts",Me_IDS,79.3/2.0,79.5/2.0,(222.7-64.0)/2.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_connection_parts_R = geom->MakeCone("IDS_connection_parts_R",Me_IDS,(64.0-56.0)/2.0,42.8/2.0,43.0/2.0,79.3/2.0,79.5/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *IDS_connection_parts_L = geom->MakeCone("IDS_connection_parts_L",Me_IDS,(64.0-56.0)/2.0,79.3/2.0,79.5/2.0,42.8/2.0,43.0/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *beampipe_central_part       = geom->MakeTube("beampipe_central_part",Me_beampipe,4.05/2.0,4.15/2.0,141.5);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_side_parts         = geom->MakeTube("beampipe_side_parts",Me_beampipe,9.52/2.0,9.62/2.0,100.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_connection_parts_R = geom->MakeCone("beampipe_connection_parts_R",Me_beampipe,(191.5-141.5)/2.0,4.05/2.0,4.15/2.0,9.52/2.0,9.62/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *beampipe_connection_parts_L = geom->MakeCone("beampipe_connection_parts_L",Me_beampipe,(191.5-141.5)/2.0,9.52/2.0,9.62/2.0,4.05/2.0,4.15/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *Pixel_support       = geom->MakeTube("Pixel_support",Me_Pixel_support,21.8/2.0,22.0/2.0,56.0);  // r_min, r_max, dz (half of total length)

    inner_field_tube       ->SetLineColor(4);
    outer_field_tube       ->SetLineColor(kRed-8);
    IDS_central_part       ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_side_parts         ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_R ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_L ->SetLineColor(2);  // Inner Detector Support (IDS)

    beampipe_central_part       ->SetLineColor(3);  // (beampipe)
    beampipe_side_parts         ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_R ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_L ->SetLineColor(3);  // (beampipe)

    Pixel_support ->SetLineColor(kYellow-3);  // (pixel support)

    top->AddNodeOverlap(inner_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,64.0+(222.7-64.0)/2.0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,-(64.0+(222.7-64.0)/2.0)));
    top->AddNodeOverlap(IDS_connection_parts_R,1,new TGeoTranslation(0,0,56.0+(64.0-56.0)/2.0));
    top->AddNodeOverlap(IDS_connection_parts_L,1,new TGeoTranslation(0,0,-(56.0+(64.0-56.0)/2.0)));

    top->AddNodeOverlap(beampipe_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,191.5+100.0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,-(191.5+100.0)));
    top->AddNodeOverlap(beampipe_connection_parts_R,1,new TGeoTranslation(0,0,141.4+(191.5-141.5)/2.0));
    top->AddNodeOverlap(beampipe_connection_parts_L,1,new TGeoTranslation(0,0,-(141.4+(191.5-141.5)/2.0)));

    top->AddNodeOverlap(Pixel_support,1,new TGeoTranslation(0,0,0));

    top->DrawClone("ogl");



    const Int_t n_TPC_points = 50;
    TPolyLine3D   *TPC_endcaps[4];
    TPolyLine3D   *TPC_tube[4];
    TPolyLine3D   *TPC_tube_lines[n_TPC_points+1];

    Float_t radius_table[4] = {200,200,3.81,3.81};
    Float_t z_val_table[4]  = {200,-200,200,-200};

    Float_t radius_table_tube[4] = {50,50,50,50};
    Float_t z_val_table_tube[4]  = {200,-200,100,-100};

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
        TPC_endcaps[r]->SetLineColor(28); // 28
        TPC_endcaps[r]->SetLineWidth(2);
        TPC_endcaps[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_tube[r] = new TPolyLine3D();
        Float_t radius   = radius_table_tube[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table_tube[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_tube[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
            if(r == 0 && (t%4 == 0))
            {
                TPC_tube_lines[t] = new TPolyLine3D();
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_val_table_tube[r+1]);
                TPC_tube_lines[t]->SetLineStyle(0);
                TPC_tube_lines[t]->SetLineColor(28); // 28
                TPC_tube_lines[t]->SetLineWidth(1);
                //TPC_tube_lines[t]->DrawClone("ogl");
            }
        }
        TPC_tube[r]->SetLineStyle(0);
        TPC_tube[r]->SetLineColor(28); // 28
        TPC_tube[r]->SetLineWidth(2);
        TPC_tube[r]->DrawClone("ogl");
    }

    TPolyLine3D   *BeamLine;
    BeamLine       = new TPolyLine3D(2);
    BeamLine   ->SetPoint(0,0,0,-550);
    BeamLine   ->SetPoint(1,0,0,550);
    BeamLine   ->SetLineStyle(0);
    BeamLine   ->SetLineColor(4);
    BeamLine   ->SetLineWidth(2);
    BeamLine   ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TH1F* Modify_hist(TH1F* hist_input, Double_t mean_shift, Double_t width_fac, Double_t width_shift, Long64_t sample_val)
{
    // Modifies an input histogram by scaling it with width_fac (with scaling center width_shift), shifting it with mean_shift
    // sample val is a threshold at which a pure sampling of bin contents is changed to a sampling with max of sample_val and weights
    // this is important if the histogram spans over many magnitudes -> pure sampling would take too much time

    TString histname = hist_input->GetName();
    histname += "_out";
    TH1F* hist_output = (TH1F*)hist_input->Clone(histname.Data());
    hist_output->Reset();

    TRandom3 r3b_hist;
    r3b_hist.SetSeed(0); // seed for random number generator changes every second

    // Get Integral
    Double_t integral_input = hist_input ->Integral(1,hist_input->GetNbinsX());

    // Get lowest bin content
    Double_t bin_cont_min = hist_input ->GetBinContent(hist_input->GetMaximumBin());
    for(Int_t ibin = 1; ibin <= hist_input->GetNbinsX(); ibin++)
    {
        Double_t bin_cont  = hist_input ->GetBinContent(ibin);
        if(bin_cont > 0.0 && bin_cont < bin_cont_min) bin_cont_min = bin_cont;
    }

    for(Int_t ibin = 1; ibin <= hist_input->GetNbinsX(); ibin++)
    {
        Double_t bin_cont  = hist_input ->GetBinContent(ibin);
        Double_t bin_width = hist_input ->GetBinWidth(ibin);
        Double_t bin_cent  = hist_input ->GetBinCenter(ibin);

        Double_t weight_sample = 1.0;
        if(bin_cont_min > 0.0) bin_cont /= bin_cont_min;
        if(bin_cont > sample_val && sample_val > 0)
        {
            weight_sample = ((Double_t)bin_cont)/((Double_t)sample_val);
            bin_cont      = sample_val;
        }

        for(Long64_t iloop = 0; iloop < (Long64_t)bin_cont; iloop++)
        {
            Double_t random = r3b_hist.Rndm();
            random *= bin_width;
            random += bin_cent - (bin_width/2.0) - width_shift;
            random *= width_fac;
            random += mean_shift + width_shift;
            hist_output ->Fill(random,weight_sample);
        }
    }

    if(bin_cont_min > 0.0) hist_output->Scale(bin_cont_min);

    return hist_output;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_Detector_2D(Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1,
                             const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1,
                             Float_t x_offset = 0.0, Float_t y_offset = 0.0
                            )
{
    Float_t z = 0.0;
    const Int_t n_points = 50;
    TPolyLine   *tp_Circles[n_radii];
    TPolyLine   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Circles[r]->SetLineStyle(line_style);
        tp_Circles[r]->SetLineColor(color); // 28
        tp_Circles[r]->SetLineWidth(line_width);
        tp_Circles[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Radial[r]->SetLineStyle(0);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(1);
        tp_Radial[r]->DrawClone("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Decay_mom(Double_t M_mother, Double_t M_daughterA, Double_t M_daughterB)
{
    Double_t mom = 0.0;

    Double_t E_term = M_mother*M_mother - M_daughterB*M_daughterB + M_daughterA*M_daughterA;
    Double_t E_daughterA = 0.0;
    if(E_term > 0.0)
    {
        E_daughterA = E_term/(2.0*M_mother);
    }
    else return -1.0;

    mom = TMath::Sqrt(E_daughterA*E_daughterA - M_daughterA*M_daughterA);

    return mom;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TGraphAsymmErrors* Calc_feed_down_v2(TGraphAsymmErrors* tgae_input_v2, Int_t PID, Int_t Energy, Double_t T_BW, Double_t rho0_BW,
                        Double_t rhoa_BW, Double_t s2_BW)
{
    // PID:
    // 0 = pi
    // 1 = charged K
    // 2 = p
    // 3 = anti-p
    // 4 = Lambda
    // 5 = anti-Lambda
    // 6 = K0S

    // Energy:
    // 0 = 7.7 GeV
    // 1 = 11.5 GeV
    // 2 = 19.6 GeV
    // 3 = 27 GeV
    // 4 = 39 GeV
    // 5 = 62.4 GeV
    // 6 = 200 GeV
    // 7 = 2.76 TeV

    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);

    const Int_t N_Energies           = 8;
    const Int_t N_Particles_v2_input = 7;  // 0     1   2   3      4            5      6     7     8    9      10        11          12            13     14    15   16   17       18            19          20    21   22    23      24
    const Int_t N_Particles          = 25; // pi+, phi, K-, K+, anti-Lambda, anti-Xi, K0S, Lambda, Xi-, rho, eta(958), omega(782), Delta(1232)++, Sigma0, Xi0, Omega, p, anti-p, anti-Sigma0, anti-Delta++, gamma, phi, N*  eta(547)  N*
    const Int_t N_feed_down_channels = 7;  // max value for all particles
    const Int_t N_decay_products     = 3;  // max number of decay daughters

    Int_t PID_in_table[N_Particles_v2_input] = {0,2,16,17,7,4,6};
    Int_t PID_in = PID_in_table[PID];

    const Int_t feed_down_PID_table[N_Particles_v2_input][N_feed_down_channels] =
    {
        //{9,12,22,23,11,10,24}, // pi
        {24,10,23,12,22,11,-1}, // pi   // removed 22 (maybe not so bad), 11 is the bad guy
        {1,-1,-1,-1,-1,-1,-1}, // K
        {7,12,22,-1,-1,-1,-1}, // p
        {4,19,-1,-1,-1,-1,-1}, // anti-p
        {8,13,14,15,-1,-1,-1}, // Lambda
        {5,18,-1,-1,-1,-1,-1}, // anti-Lambda
        {21,-1,-1,-1,-1,-1,-1}, // K0S
    };

    const Int_t decay_products[N_Particles][N_decay_products] =
    {
        {-1,-1,-1}, // pi
        {2,3,-1},   // phi -> K+ + K-
        {-1,-1,-1}, // K-
        {-1,-1,-1}, // K+
        {0,17,-1}, // anti-Lambda -> pi + anti-p
        {0,4,-1}, // anti-Xi -> pi + anti-Lambda
        {0,0,-1}, // K0S -> pi + pi
        {0,16,-1}, // Lambda -> pi + proton
        {0,7,-1}, // Xi- -> pi + Lambda
        {0,0,-1}, // rho -> pi + pi
        {0,0,23}, // eta(958) -> pi + pi + eta
        {0,0,0}, // omega(782) -> pi+ + pi- + pi0
        {0,16,-1}, // Delta(1232)++ -> pi + proton
        {20,7,-1}, // Sigma0 -> gamma + Lambda
        {0,7,-1}, // Xi0 -> pi + Lambda
        {2,7,-1}, // Omega -> K + Lambda
        {-1,-1,-1}, // p
        {-1,-1,-1}, // anti-p
        {20,4,-1}, // anti-Sigma0 -> gamma + anti-Lambda
        {0,17,-1}, // anti-Delta++ -> pi + anti-proton
        {0,0,-1}, // gamma
        {6,6,-1}, // phi -> K0S + K0L
        {0,16,-1}, // N* -> pi + p
        {0,0,0}, // eta(547) -> pi+ + pi- + pi0
        {0,0,16}, // N* -> pi + pi + p
    };
    //                                                   0    1    2   3    4     5     6     7     8     9    10    11    12    13    14     15  16  17   18     19  20   21    22     23     24
    const Double_t decay_momentum_table[N_Particles] = {0.0,0.127,0.0,0.0,0.101,0.139,0.206,0.101,0.139,0.363,0.232,0.327,0.229,0.074,0.135,0.211,0.0,0.0,0.074,0.229,0.0,0.110,0.405, 0.174, 0.380};
    //                                          pi+     phi      K-      K+   anti-Lambda anti-Xi  K0S     Lambda    Xi-    rho    eta(958) omega Delta  Sigma0    Xi0    Omega     p      anti-p   aSigma0 aDelta gamma phi   N*   eta(547) N*
    const Double_t mass_table[N_Particles] = {0.13957,1.01946,0.493677,0.493677,1.115683,1.32131,0.497648,1.115683,1.32131,0.7755,0.95778,0.78265,1.232,1.192642,1.31483,1.67245,0.938272,0.938272,1.192642,1.232,0.0,1.01946,1.480,0.54751,1.480};

    const Long64_t N_sample = 200000;

    // Particle yields from SHM predictions -> THERMUS (Macro_THERMUS_predict.cc)
    Double_t Res_scale        = 1.0;  // rescale for N*
    Double_t Res_scale_NStar  = 0.4;  // rescale for N* 0.4
    Double_t Res_scale_Lambda = 0.2;  // rescale for Lambda   -> estimated dca efficiency
    Double_t Res_scale_Omega  = 0.4;  // rescale for Omega    -> estimated dca efficiency
    Double_t Res_scale_Xi     = 0.6;  // rescale for Xi       -> estimated dca efficiency
    Double_t Res_rho          = 0.5;  // rescale for rho      -> in medium decays, less flow
    Double_t Res_eta          = 0.23; // rescale for eta(547) -> branching ratio
    Double_t Res_omega        = 0.3; // rescale for eta(547) -> branching ratio 0.89 but some decay in medium due to short live time
    Double_t Res_eta_958      = 0.45; // rescale for eta(958) -> branching ratio
    Double_t particle_yields[N_Energies][N_Particles] =
    {
        // 7.7
        {18.842, 0.0809822, 0.984182, 2.03403, Res_scale_Lambda*0.012904, Res_scale_Xi*0.00180397, 1.49986, Res_scale_Lambda*1.97031, Res_scale_Xi*0.110732, Res_rho*1.08777, Res_eta_958*0.0766547,Res_omega*0.883634, 1.58338, 0.395304, 0.110232, Res_scale_Omega*0.00556682, 11.2899, 0.0307144,0.00256662,0.00469493,0.0,0.0809822,Res_scale*1.58338,Res_eta*0.529146,Res_scale_NStar*1.58338},
        // 11.5
        {24.7732,0.166343,1.92581,3.05456,Res_scale_Lambda*0.0665792,Res_scale_Xi*0.00848224,2.45721,Res_scale_Lambda*2.20594,Res_scale_Xi*0.154628,Res_rho*1.89612,Res_eta_958*0.152232,Res_omega*1.59568,1.47253,0.416078,0.155051,Res_scale_Omega*0.0105926,10.252,0.173699,0.0124554,0.0263733,0.0,0.166343,Res_scale*1.47253,Res_eta*0.843459,Res_scale_NStar*1.58338},
        // 19.6
        {31.4526,0.260379,3.02932,3.90747,Res_scale_Lambda*0.221877,Res_scale_Xi*0.0257517,3.41128,Res_scale_Lambda*1.98909,Res_scale_Xi*0.164424,Res_rho*2.69009,Res_eta_958*0.233442,Res_omega*2.30929,1.14568,0.358709,0.165846,Res_scale_Omega*0.0140665,7.91359,0.635199,0.0397831,0.0947728,0.0,0.260379,Res_scale*1.14568,Res_eta*1.13144,Res_scale_NStar*1.58338},
        // 27
        {33.7056,0.296082,3.42981,4.23668,Res_scale_Lambda*0.348153,Res_scale_Xi*0.0395389,3.76868,Res_scale_Lambda*1.7408,Res_scale_Xi*0.149879,Res_rho*2.97559,Res_eta_958*0.26391,Res_omega*2.56822,0.963581,0.309601,0.151293,Res_scale_Omega*0.0135663,6.66259,1.019,0.0616097,0.151298,0.0,0.296082,Res_scale*0.963581,Res_eta*1.23154,Res_scale_NStar*1.58338},
        // 39
        {35.2158,0.3212,3.70381,4.47285,Res_scale_Lambda*0.497978,Res_scale_Xi*0.0558188,4.01892,Res_scale_Lambda*1.48852,Res_scale_Xi*0.131371,Res_rho*3.17241,Res_eta_958*0.285248,Res_omega*2.7473,0.802959,0.262355,0.132621,Res_scale_Omega*0.0123063,5.56283,1.47752,0.0873729,0.218829,0.0,0.3212,Res_scale*0.802959,Res_eta*1.29967,Res_scale_NStar*1.58338},
        // 62.4
        {36.2148,0.337546,3.91726,4.58277,Res_scale_Lambda*0.659449,Res_scale_Xi*0.0724835,4.17698,Res_scale_Lambda*1.27079,Res_scale_Xi*0.115194,Res_rho*3.29889,Res_eta_958*0.299095,Res_omega*2.86261,0.667969,0.222697,0.116366,Res_scale_Omega*0.0111575,4.62849,1.99556,0.115131,0.294752,0.0,0.337546,Res_scale*0.667969,Res_eta*1.34311,Res_scale_NStar*1.58338},
        // 200
        {36.9634,0.348295,4.15469,4.54569,Res_scale_Lambda*0.877946,Res_scale_Xi*0.0925872,4.27385,Res_scale_Lambda*1.03121,Res_scale_Xi*0.0978851,Res_rho*3.38142,Res_eta_958*0.308185,Res_omega*2.93795,0.520315,0.179972,0.0990989,Res_scale_Omega*0.00999385,3.59353,2.76675,0.152901,0.406583,0.0,0.348295,Res_scale*0.520315,Res_eta*1.37131,Res_scale_NStar*1.58338},
        // 2760
        {36.9634,0.348295,4.15469,4.54569,Res_scale_Lambda*0.877946,Res_scale_Xi*0.0925872,4.27385,Res_scale_Lambda*1.03121,Res_scale_Xi*0.0978851,Res_rho*3.38142,0.308185,Res_omega*2.93795,0.520315,0.179972,0.0990989,Res_scale_Omega*0.00999385,3.59353,2.76675,0.152901,0.406583,0.0,0.348295,Res_scale*0.520315,Res_eta*1.37131,Res_scale_NStar*1.58338} // copy of 200 GeV
    };

    // Fits to preliminary BES data
    // T = 0.153 for K+ at 11.5
    // T = 0.207 for p  at 11.5

    // T = 0.2 for K+ at 62.4
    // T = 0.2 for p  at 62.4

    Double_t T_kin[N_Energies] = {0.18,0.2,0.2,0.2,0.2,0.2,0.25,0.3};

    TF1 *FlowFit         = new TF1("FlowFit",FlowFitFunc,0,3.15,5);
    TF1 *PtFit2_mod_x    = new TF1("PtFit2_mod_x",PtFitFunc2_mod_x,0.0,2.5,4);
    TF1* BlastWaveFit    = new TF1("BlastWaveFit",BlastWaveFitFunc_no_mass_array,0.0,2.5,6);
    for(Int_t i = 0; i < 6; i++)
    {
        BlastWaveFit ->SetParameter(i,0.0);
        BlastWaveFit ->SetParError(i,0.0);
    }
    for(Int_t i = 0; i < 4; i++)
    {
        PtFit2_mod_x ->SetParameter(i,0.0);
        PtFit2_mod_x ->SetParError(i,0.0);
    }

    PtFit2_mod_x ->SetParameter(2,100.0);         // amplitude
    PtFit2_mod_x ->SetParameter(3,0.0);           // shift

    for(Int_t x = 0; x < 5; x++)
    {
        FlowFit->ReleaseParameter(x);
        FlowFit->SetParError(x,0.0);
        FlowFit->SetParameter(x,0.0);
    }

    FlowFit->SetParameter(0,1.0);
    FlowFit->SetParameter(1,0.0);
    FlowFit->SetParameter(2,0.0);
    FlowFit->SetParameter(3,0.0);
    FlowFit->SetParameter(4,0.0);

    BlastWaveFit ->SetParameter(0,T_BW);
    BlastWaveFit ->SetParameter(1,rho0_BW);
    BlastWaveFit ->SetParameter(2,rhoa_BW);
    BlastWaveFit ->SetParameter(3,s2_BW);
    BlastWaveFit ->SetParameter(4,1.0);

    TGraphAsymmErrors* tgae_output_v2;
    tgae_output_v2 = (TGraphAsymmErrors*)tgae_input_v2->Clone("tgae_output_v2");
    TProfile* tp_v2_feed_down = new TProfile("tp_v2_feed_down","tp_v2_feed_down",20,0,2.5);

    // Loop over the feed down channels
    TLorentzVector pMother, pDaughterA, pDaughterB, pDaughterAB, pDaughterC;
    Double_t total_weigh_factor = 0.0;
    Double_t input_weigh_factor = particle_yields[Energy][PID_in];
    for(Int_t iFeed_down = 0; iFeed_down < N_feed_down_channels; iFeed_down++)
    {
        Int_t PID_mother    = feed_down_PID_table[PID][iFeed_down];
        Int_t PID_daughterA = decay_products[PID_mother][0];
        Int_t PID_daughterB = decay_products[PID_mother][1];
        Int_t PID_daughterC = decay_products[PID_mother][2];
        if(PID_mother == -1) break;
        Double_t Mass_mother    = mass_table[PID_mother];
        Double_t Mass_daughterA = mass_table[PID_daughterA];
        Double_t Mass_daughterB = mass_table[PID_daughterB];
        Double_t Mass_daughterC = -1.0;
        Int_t flag_three_body = 0; // flag for three body decay
        if(PID_daughterC >= 0)
        {
            Mass_daughterC = mass_table[PID_daughterC];
            flag_three_body = 1;
        }
        Double_t decay_momentum = decay_momentum_table[PID_mother];
        BlastWaveFit ->SetParameter(5,Mass_mother);
        PtFit2_mod_x ->SetParameter(1,T_kin[Energy]+Mass_mother/20.0); // Effective temperature -> small mass dependence
        PtFit2_mod_x ->SetParameter(0,Mass_mother);
        Double_t weigh_factor = particle_yields[Energy][PID_mother]; // yield from SHM

        total_weigh_factor += weigh_factor;

        //if(PID_daughterA == PID_in && PID_daughterB == PID_in) // double the weigh since both particles are used, BUT only one is analyzed...
        //{
        //     total_weigh_factor += weigh_factor;
        //}

#if 1
        cout << "iFeed_down = " << iFeed_down << ", PID_mother = " << PID_mother << ", PID_daughterA = " << PID_daughterA  << ", PID_daughterB = " << PID_daughterB
            << ", Mass_mother = " << Mass_mother << ", Mass_daughterA = " << Mass_daughterA << ", Mass_daughterB = " << Mass_daughterB << ", decay_momentum = " << decay_momentum
            << ", input_weigh_factor = " << input_weigh_factor << ", weigh_factor = " << weigh_factor
            << endl;
#endif



        //TH1F* h_Mother_pt = new TH1F("h_Mother_pt","h_Mother_pt",100,0,1.5);
        //for(Int_t iSample = 0; iSample < 1000000; iSample++)
        //{
        //    Double_t pt_val_mother = PtFit2_mod_x ->GetRandom();
        //    h_Mother_pt ->Fill(pt_val_mother);
        //}

        // Sampling
        for(Int_t iSample = 0; iSample < N_sample; iSample++)
        {
            if (iSample != 0  &&  iSample % 200 == 0)
                cout << "." << flush;
            if (iSample != 0  &&  iSample % 2000 == 0)
            {
                Float_t event_percent = 100.0*iSample/N_sample;
                cout << " " << iSample << " (" << event_percent << "%) " << "\n" << "==> Processing data, " << flush;
            }


            Double_t pt_val_mother = PtFit2_mod_x ->GetRandom();
            //Double_t pt_val_mother =  h_Mother_pt->GetRandom();
            Double_t v2_val_mother = BlastWaveFit ->Eval(pt_val_mother);
            FlowFit->SetParameter(2,v2_val_mother);
            Double_t Psi_phi_val_mother = FlowFit ->GetRandom();
            //cout << "Psi_phi_val_mother = " << Psi_phi_val_mother << endl;

            pMother.SetXYZM(pt_val_mother,0.0,0.0,Mass_mother); // Mother particle

            if(flag_three_body == 0) // two body decay
            {
                pDaughterA.SetXYZM(decay_momentum,0.0,0.0,Mass_daughterA);
                pDaughterB.SetXYZM(-decay_momentum,0.0,0.0,Mass_daughterB);

                Double_t anglex = r3b.Rndm()*TMath::Pi()*2;
                Double_t angley = r3b.Rndm()*TMath::Pi()*2;
                Double_t anglez = r3b.Rndm()*TMath::Pi()*2;

                pDaughterA.RotateZ(anglez);
                pDaughterA.RotateX(anglex);
                pDaughterA.RotateY(angley);

                pDaughterB.RotateZ(anglez);
                pDaughterB.RotateX(anglex);
                pDaughterB.RotateY(angley);

                pDaughterA.Boost(pMother.BoostVector());
                pDaughterB.Boost(pMother.BoostVector());

                Double_t AngleA = (pMother.BoostVector()).DeltaPhi(pDaughterA.BoostVector());
                Double_t AngleB = (pMother.BoostVector()).DeltaPhi(pDaughterB.BoostVector());

                Double_t PtA = pDaughterA.Pt();
                Double_t PtB = pDaughterB.Pt();

                //cout << "AngleA = " << AngleA << ", AngleB = " << AngleB << ", PtA = " << PtA << ", PtB = " << PtB << endl;

                AngleA = Psi_phi_val_mother - AngleA;
                AngleB = Psi_phi_val_mother - AngleB;

                if(AngleA < 0.0) AngleA += TMath::Pi();
                if(AngleA > TMath::Pi()) AngleA -= TMath::Pi();

                if(AngleB < 0.0) AngleB += TMath::Pi();
                if(AngleB > TMath::Pi()) AngleB -= TMath::Pi();

                Double_t cosA = TMath::Cos(2.0*AngleA);
                Double_t cosB = TMath::Cos(2.0*AngleB);

                //tp_v2_feed_down ->Fill(pt_val_mother,v2_val_mother); // test, working
                //tp_v2_feed_down ->Fill(pt_val_mother,TMath::Cos(2.0*Psi_phi_val_mother)); // test, working


                if(PID_daughterA == PID_in)
                {
                    tp_v2_feed_down ->Fill(PtA,cosA,weigh_factor);
                }
                if(PID_daughterB == PID_in)
                {
                    tp_v2_feed_down ->Fill(PtB,cosB,weigh_factor);
                }
            }

            if(flag_three_body == 1) // three body decay
            {
                // first decay
                Double_t low_limit_Mass_first = Mass_daughterA + Mass_daughterB;
                Double_t up_limit_Mass_first  = Mass_mother - Mass_daughterC;
                Double_t Mass_daughter_first  = r3b.Rndm()*(up_limit_Mass_first - low_limit_Mass_first) + low_limit_Mass_first;

                if(up_limit_Mass_first > low_limit_Mass_first)
                {
                    Double_t Decay_mom_first   = Decay_mom(Mass_mother,Mass_daughterC,Mass_daughter_first);

                    pDaughterAB.SetPxPyPzE(Decay_mom_first,0.0,0.0,TMath::Sqrt(Mass_daughterC*Mass_daughterC+Decay_mom_first*Decay_mom_first));
                    pDaughterC.SetPxPyPzE(-Decay_mom_first,0.0,0.0,TMath::Sqrt(Mass_daughter_first*Mass_daughter_first+Decay_mom_first*Decay_mom_first));

                    Double_t anglex = r3b.Rndm()*TMath::Pi()*2;
                    Double_t angley = r3b.Rndm()*TMath::Pi()*2;
                    Double_t anglez = r3b.Rndm()*TMath::Pi()*2;

                    pDaughterAB.RotateZ(anglez);
                    pDaughterAB.RotateX(anglex);
                    pDaughterAB.RotateY(angley);

                    pDaughterC.RotateZ(anglez);
                    pDaughterC.RotateX(anglex);
                    pDaughterC.RotateY(angley);

                    pDaughterAB.Boost(pMother.BoostVector());
                    pDaughterC.Boost(pMother.BoostVector());


                    // second decay
                    Double_t Decay_mom_second   = Decay_mom(Mass_daughter_first,Mass_daughterA,Mass_daughterB);


                    pDaughterA.SetXYZM(Decay_mom_second,0.0,0.0,Mass_daughterA);
                    pDaughterB.SetXYZM(-Decay_mom_second,0.0,0.0,Mass_daughterB);

                    anglex = r3b.Rndm()*TMath::Pi()*2;
                    angley = r3b.Rndm()*TMath::Pi()*2;
                    anglez = r3b.Rndm()*TMath::Pi()*2;

                    pDaughterA.RotateZ(anglez);
                    pDaughterA.RotateX(anglex);
                    pDaughterA.RotateY(angley);

                    pDaughterB.RotateZ(anglez);
                    pDaughterB.RotateX(anglex);
                    pDaughterB.RotateY(angley);

                    pDaughterA.Boost(pDaughterAB.BoostVector());
                    pDaughterB.Boost(pDaughterAB.BoostVector());

                    Double_t AngleA = (pMother.BoostVector()).DeltaPhi(pDaughterA.BoostVector());
                    Double_t AngleB = (pMother.BoostVector()).DeltaPhi(pDaughterB.BoostVector());
                    Double_t AngleC = (pMother.BoostVector()).DeltaPhi(pDaughterC.BoostVector());

                    Double_t PtA = pDaughterA.Pt();
                    Double_t PtB = pDaughterB.Pt();
                    Double_t PtC = pDaughterC.Pt();

                    //cout << "AngleA = " << AngleA << ", AngleB = " << AngleB << ", PtA = " << PtA << ", PtB = " << PtB << endl;

                    AngleA = Psi_phi_val_mother - AngleA;
                    AngleB = Psi_phi_val_mother - AngleB;
                    AngleC = Psi_phi_val_mother - AngleC;

                    if(AngleA < 0.0) AngleA += TMath::Pi();
                    if(AngleA > TMath::Pi()) AngleA -= TMath::Pi();

                    if(AngleB < 0.0) AngleB += TMath::Pi();
                    if(AngleB > TMath::Pi()) AngleB -= TMath::Pi();

                    if(AngleC < 0.0) AngleC += TMath::Pi();
                    if(AngleC > TMath::Pi()) AngleC -= TMath::Pi();

                    Double_t cosA = TMath::Cos(2.0*AngleA);
                    Double_t cosB = TMath::Cos(2.0*AngleB);
                    Double_t cosC = TMath::Cos(2.0*AngleC);

                    //tp_v2_feed_down ->Fill(pt_val_mother,v2_val_mother); // test, working
                    //tp_v2_feed_down ->Fill(pt_val_mother,TMath::Cos(2.0*Psi_phi_val_mother)); // test, working


                    if(PID_daughterA == PID_in)
                    {
                        tp_v2_feed_down ->Fill(PtA,cosA,weigh_factor);
                    }
                    if(PID_daughterB == PID_in)
                    {
                        tp_v2_feed_down ->Fill(PtB,cosB,weigh_factor);
                    }
                    if(PID_daughterC == PID_in)
                    {
                        tp_v2_feed_down ->Fill(PtC,cosC,weigh_factor);
                    }

                }
            }

        }
    }

    //Double_t diff_weigh_factor = input_weigh_factor - total_weigh_factor;
    Double_t diff_weigh_factor = input_weigh_factor;
    if(diff_weigh_factor <= 0.0) diff_weigh_factor = 1.0;

    if(PID == 0) // special treatment for pions, resonances have a 60% contribution, see http://arxiv.org/pdf/hep-ph/0407174.pdf, Fig. 1
    {
        total_weigh_factor = 0.6;
        input_weigh_factor = 0.4;
        diff_weigh_factor  = input_weigh_factor;
    }

    for(Int_t iPoint = 0; iPoint < tgae_output_v2->GetN(); iPoint++)
    {
        Double_t x_val, y_val;
        tgae_output_v2 ->GetPoint(iPoint,x_val,y_val);
        Double_t exl_ME = tgae_input_v2 ->GetErrorXlow(iPoint);
        Double_t exh_ME = tgae_input_v2 ->GetErrorXhigh(iPoint);
        Double_t eyl_ME = tgae_input_v2 ->GetErrorYlow(iPoint);
        Double_t eyh_ME = tgae_input_v2 ->GetErrorYhigh(iPoint);

        Double_t y_val_FD     = tp_v2_feed_down->GetBinContent(tp_v2_feed_down->FindBin(x_val));
        Double_t y_val_err_FD = tp_v2_feed_down->GetBinError(tp_v2_feed_down->FindBin(x_val));

        //Double_t y_val_out = (y_val*input_weigh_factor - y_val_FD*total_weigh_factor)/diff_weigh_factor;
        Double_t y_val_out = (y_val*(input_weigh_factor+total_weigh_factor) - y_val_FD*total_weigh_factor)/diff_weigh_factor;

        tgae_output_v2 ->SetPoint(iPoint,x_val,y_val_out); // original
        //tgae_output_v2 ->SetPoint(iPoint,x_val,y_val_FD);  // test
        //tgae_output_v2 ->SetPointError(iPoint,0.0,0.0,y_val_err_FD,y_val_err_FD); // test

        cout << "iPoint = " << iPoint << ", y_val_FD = " << y_val_FD << ", y_val_in = " << y_val << ", y_val_out = " << y_val_out
            << ", input_weigh_factor = " << input_weigh_factor << ", total_weigh_factor = " << total_weigh_factor << endl;
    }

    return tgae_output_v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Calcualtes the maximum/minimum envelope curve as systematic error for any number of added TGraphAsymmErrors/TH1F
// y-poins are the mean of the maximum/minimum
// Graphs need to have the same points in x
class Calc_syst_error_graph_hist
{
    TGraphAsymmErrors* tgae_stat;
    TGraphAsymmErrors* tgae_syst;
    vector<TGraphAsymmErrors*> vec_tgae;
    vector<TH1F*> vec_hist;
    Int_t flag_input, flag_add_stat_error;
public:
    void add_graph(TGraphAsymmErrors*);
    void add_hist(TH1F*);
    void add_stat_error(Int_t);
    void calculate_error();
    TGraphAsymmErrors* get_graph_stat();
    TGraphAsymmErrors* get_graph_syst();
    void clear();
};


void Calc_syst_error_graph_hist::add_graph(TGraphAsymmErrors* graph)
{
    vec_tgae.push_back(graph);
    flag_input = 0;
    flag_add_stat_error = 0;
}

void Calc_syst_error_graph_hist::add_hist(TH1F* hist)
{
    vec_hist.push_back(hist);
    flag_input = 1;
    flag_add_stat_error = 0;
}

void Calc_syst_error_graph_hist::add_stat_error(Int_t add_stat)
{
    if(add_stat == 0)
    {
        flag_add_stat_error = 0;
    }
    else
    {
        flag_add_stat_error = 1;
    }
}

void Calc_syst_error_graph_hist::calculate_error()
{
    Int_t start_point, stop_point, size;
    if(flag_input == 0) // TGraphAsymmErrors
    {
        tgae_stat   = (TGraphAsymmErrors*)vec_tgae[0]->Clone("tgae_stat");
        tgae_syst   = (TGraphAsymmErrors*)vec_tgae[0]->Clone("tgae_syst");
        start_point = 0;
        stop_point  = vec_tgae[0]->GetN();
        size        = vec_tgae.size();
    }
    if(flag_input == 1) // TH1F
    {
        tgae_stat   = new TGraphAsymmErrors(vec_hist[0]);
        tgae_stat   ->SetName("tgae_stat");
        tgae_syst   = new TGraphAsymmErrors(vec_hist[0]);
        tgae_stat   ->SetName("tgae_syst");
        start_point = 1;
        stop_point  = vec_hist[0]->GetNbinsX()+1;
        size        = vec_hist.size();
    }

    //cout << "start_point: " << start_point << ", stop_point: " << stop_point << ", flag_input: " << flag_input << ", size: " << size << endl;

    for(Int_t i_point = start_point; i_point < stop_point; i_point++)
    {
        Double_t max_y_val, min_y_val, max_y_err, min_y_err, x_val, x_err_low, x_err_high;
        Int_t i_file_use_counter = 0;
        for(Int_t i_file = 0; i_file < size; i_file++)
        {
            Double_t x,y;
            Double_t Xhigh,Xlow,Yhigh,Ylow;

            if(flag_input == 0) // TGraphAsymmErrors
            {
                vec_tgae[i_file]->GetPoint(i_point,x,y);
                Xhigh = vec_tgae[i_file]->GetErrorXhigh(i_point);
                Xlow  = vec_tgae[i_file]->GetErrorXlow(i_point);
                Yhigh = vec_tgae[i_file]->GetErrorYhigh(i_point);
                Ylow  = vec_tgae[i_file]->GetErrorYlow(i_point);
            }
            if(flag_input == 1) // TH1F
            {
                x     = vec_hist[i_file]->GetBinCenter(i_point);
                y     = vec_hist[i_file]->GetBinContent(i_point);
                Xlow  = vec_hist[i_file]->GetBinWidth(i_point)/2.0;
                Xhigh = Xlow;
                Yhigh = vec_hist[i_file]->GetBinError(i_point);
                Ylow  = Yhigh;
                //cout << "i_point: " << i_point << ", x: " << x << ", y: " << y << ", Xlow: " << Xlow << ", Ylow: " << Ylow << endl;
            }

            //if(y - Ylow < 1E-09 && x > 5.0) continue; // avoid huge negative error fluctuations
            i_file_use_counter++;

            if(i_file_use_counter == 1)
            {
                if(flag_add_stat_error)
                {
                    max_y_val  = y + Yhigh;
                    min_y_val  = y - Ylow;
                }
                else
                {
                    max_y_val  = y;
                    min_y_val  = y;
                }
                x_val      = x;
                x_err_low  = Xlow;
                x_err_high = Xhigh;
                max_y_err  = Yhigh;
                min_y_err  = Ylow;
            }
            else
            {
                if(flag_add_stat_error)
                {
                    if( (y + Yhigh) > max_y_val)
                    {
                        max_y_val = y + Yhigh;
                        max_y_err = Yhigh;
                    }
                    if( (y - Ylow) < min_y_val)
                    {
                        min_y_val = y - Ylow;
                        min_y_err = Ylow;
                    }
                }
                else
                {
                    if( y > max_y_val)
                    {
                        max_y_val = y;
                    }
                    if( y < min_y_val)
                    {
                        min_y_val = y;
                    }

                    if(y - Ylow < 1E-09) continue; // avoid large negative fluctuations in unfolding
                    if(Yhigh > max_y_err)
                    {
                        max_y_err = Yhigh;
                    }
                    if(Ylow > min_y_err)
                    {
                        min_y_err = Ylow;
                    }
                }
            }

            //if(x > 10.0 && x < 11.2) cout << "i_file: " << i_file << ", Ylow: " << Ylow << ", Yhigh: " << Yhigh << endl;
        }
        Double_t average_y_val = (max_y_val+min_y_val)/2.0;
        tgae_syst->SetPoint(i_point,x_val,average_y_val);
        tgae_stat->SetPoint(i_point,x_val,average_y_val);
        //cout << "i_point: " << i_point << ", x_val: " << x_val << ", average_y_val: " << average_y_val << ", min_y_val: " << min_y_val << ", max_y_val: " << max_y_val
        //    << ", min_y_err: " << min_y_err << ", max_y_err: " << max_y_err << endl;
        tgae_syst->SetPointError(i_point,x_err_low,x_err_high,average_y_val-min_y_val,max_y_val-average_y_val);
        if(flag_add_stat_error)
        {
            tgae_stat->SetPointError(i_point,x_err_low,x_err_high,0.0,0.0);
        }
        else
        {
            tgae_stat->SetPointError(i_point,x_err_low,x_err_high,min_y_err,max_y_err);
        }
    }
}

TGraphAsymmErrors* Calc_syst_error_graph_hist::get_graph_stat()
{
    return tgae_stat;
}

TGraphAsymmErrors* Calc_syst_error_graph_hist::get_graph_syst()
{
    return tgae_syst;
}

void Calc_syst_error_graph_hist::clear()
{
    if(flag_input == 0) // TGraphAsymmErrors
    {
        vec_tgae.clear();
    }
    if(flag_input == 1) // TH1F
    {
        vec_hist.clear();
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Calculates the errors of a ratio of numbers a and b based on Monte-Carlo
class Monte_Carlo_Error_Ratio
{
    Double_t nom, den, nom_err, den_err;
    Long64_t N_MC;
    TRandom  ran;
    Double_t RMS[2];
    Double_t ratio_nom_den;
public:
    void set_nominator_denominator(Double_t, Double_t);
    void set_nominator_denominator_error(Double_t, Double_t);
    void set_N_MC(Long64_t);
    void calc_errors();
    void truncate_RMS(Double_t);
    Double_t get_RMS_above();
    Double_t get_RMS_below();
    void clear();
    std::vector< std::vector<Double_t> > vec_RMS;
    std::vector< std::vector<Double_t> > vec_RMS_trunc;
};

void Monte_Carlo_Error_Ratio::set_nominator_denominator(Double_t a, Double_t b)
{
    nom = a;
    den = b;
}

void Monte_Carlo_Error_Ratio::set_nominator_denominator_error(Double_t a, Double_t b)
{
    nom_err = a;
    den_err = b;
}

void Monte_Carlo_Error_Ratio::set_N_MC(Long64_t i)
{
    N_MC = i;
}

void Monte_Carlo_Error_Ratio::calc_errors()
{
    vec_RMS.resize(2); // two RMS values are calculated, above and below mean value
    vec_RMS_trunc.resize(2);
    ran.SetSeed(0); // seed for random number generator changes every second
    Long64_t N_ratios = 0;
    RMS[0] = 0.0; // two RMS values are calculated, above and below mean value
    RMS[1] = 0.0;
    if(den > 0.0)
    {
        ratio_nom_den = nom/den;
    }
    else
    {
        RMS[0] = -1;
        RMS[1] = -1;
        cout << "nom: " << nom << ", den: " << den << endl;
    }
    if(RMS[0] >= 0.0)
    {
        for(Int_t i_MC = 0; i_MC < N_MC; i_MC++)
        {
            Double_t val_nom = ran.Gaus(nom,nom_err);
            Double_t val_den = ran.Gaus(den,den_err);

            if(val_den > 0.0 && val_nom > 0.0) // take only positive values (yields)
            {
                Double_t ratio_val = val_nom/val_den;
                Double_t RMS_val   = ratio_val - ratio_nom_den; // MC ratio - average
                Int_t RMS_above_below = 0; // above
                if(RMS_val < 0.0) RMS_above_below = 1; // below
                vec_RMS[RMS_above_below].push_back(RMS_val);
                RMS[RMS_above_below] += RMS_val*RMS_val;
#if 0
                cout << "i_MC: " << i_MC << ", RMS_above_below: " << RMS_above_below << ", val_nom: " << val_nom << ", nom_err: " << nom_err << ", val_den: " << val_den << ", den_err: " << den_err
                    << ", ratio_val: " << ratio_val << ", ratio_nom_den (mean): " << ratio_nom_den << ", N_ratios: " << N_ratios
                    << ", RMS_val: " << RMS[RMS_above_below]  << endl;
#endif
            }
        }
        for(Int_t i_above_below = 0; i_above_below < 2; i_above_below++)
        {
            if(vec_RMS[i_above_below].size() > 0)
            {
                RMS[i_above_below] = TMath::Sqrt((1.0/((Double_t)vec_RMS[i_above_below].size()))*RMS[i_above_below]);
            }
            else
            {
                RMS[i_above_below] = -1;
                cout << "WARNING, RMS vector size below 1, old RMS used" << endl;
            }
        }
    }
    //cout << "RMS[0]: " << RMS[0] << endl;
}

void Monte_Carlo_Error_Ratio::truncate_RMS(Double_t nSigma)
{
    for(Int_t i_above_below = 0; i_above_below < 2; i_above_below++)
    {
        //cout << "size RMS: " << vec_RMS[i_above_below].size() << endl;
        for(Int_t i_entry = 0; i_entry < vec_RMS[i_above_below].size(); i_entry++)
        {
            //cout << "i_entry: " << i_entry << ", nSigma: " << nSigma << ", RMS average: " << RMS[i_above_below]
            //    << ", diff: " << fabs(vec_RMS[i_above_below][i_entry] - ratio_nom_den) << ", RMS: "
            //    << vec_RMS[i_above_below][i_entry] << ", ratio_nom_den: " << ratio_nom_den << endl;
            //if( fabs(vec_RMS[i_above_below][i_entry] - ratio_nom_den) < nSigma*RMS[i_above_below])
            if( fabs(vec_RMS[i_above_below][i_entry]) < nSigma*RMS[i_above_below])
            {
                //cout << "ACCEPTED" << endl;
                vec_RMS_trunc[i_above_below].push_back(vec_RMS[i_above_below][i_entry]);
            }
            else
            {
                //cout << "i_entry: " << i_entry << ", REJECTED" << endl;
            }
        }
        //cout << "size RMS trunc: " << vec_RMS_trunc[i_above_below].size() << endl;
        if(vec_RMS_trunc[i_above_below].size() > 0) // make sure you always have entries
        {
            RMS[i_above_below] = 0.0;
            vec_RMS[i_above_below].clear();
            for(Int_t i_entry = 0; i_entry < vec_RMS_trunc[i_above_below].size(); i_entry++)
            {
                RMS[i_above_below] += vec_RMS_trunc[i_above_below][i_entry]*vec_RMS_trunc[i_above_below][i_entry];
                vec_RMS[i_above_below].push_back(vec_RMS_trunc[i_above_below][i_entry]);
            }
            RMS[i_above_below] = TMath::Sqrt((1.0/((Double_t)vec_RMS_trunc[i_above_below].size()))*RMS[i_above_below]);
        }
        vec_RMS_trunc[i_above_below].clear();
    }
}

Double_t Monte_Carlo_Error_Ratio::get_RMS_above()
{
    return RMS[0];
}
Double_t Monte_Carlo_Error_Ratio::get_RMS_below()
{
    return RMS[1];
}
void Monte_Carlo_Error_Ratio::clear()
{
    for(Int_t i_above_below = 0; i_above_below < 2; i_above_below++)
    {
        vec_RMS[i_above_below].clear();
        vec_RMS_trunc[i_above_below].clear();
    }
}
//----------------------------------------------------------------------------------------



