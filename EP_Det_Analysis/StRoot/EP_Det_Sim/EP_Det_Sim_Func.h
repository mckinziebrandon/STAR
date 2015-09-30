
#ifndef EP_DET_SIM_FUNC_H
#define EP_DET_SIM_FUNC_H

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

Double_t v2_pT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2 vs. pT
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT, a, b, c, d, n, f;
    pT = x_val[0];
    n  = par[0]; // number-of-constituent quarks
    a  = par[1];
    b  = par[2];
    c  = par[3];
    d  = par[4];
    f  = par[5];

    if(c != 0.0)
    {
        v2 = f*(1.0-TMath::Exp(-pT))*(a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n);
    }
    else v2 = 0.0;

    return v2;
}

Double_t calc_phi_event_plane_2nd(Double_t mQx, Double_t mQy)
{
    // calculates the angle out of a 2D vector (mQx,mQy)
    // angle is in the range from 0..pi
    if(mQx == 0.0 && mQy == 0.0) return -400.0;
    Double_t Psi = 0.0;

    Psi = TMath::ATan2(mQy,mQx);
    Psi = Psi/2.0;

    return Psi;
}

Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5, par6, par7;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    par6  = par[6];
    par7  = par[7];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x + par6*x*x*x*x*x*x + par7*x*x*x*x*x*x*x;
    return y;
}

Double_t dNdetaFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for dN/deta
    // Taken from phobos publication: http://arXiv.org/pdf/1011.1940.pdf (equation 12)
    Double_t eta, dNdeta, par0, c, alpha, beta, a;
    a      = par[0];
    alpha  = par[1];
    beta   = par[2];
    c      = par[3];
    eta    = x_val[0];
    dNdeta = c*TMath::Sqrt(1.0-1.0/TMath::Power(alpha*TMath::CosH(eta),2))/(1.0+TMath::Exp((fabs(eta)-beta)/a));
    return dNdeta;
}

Double_t dNdetaModFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for dN/deta
    // Taken from phobos publication: http://arXiv.org/pdf/1011.1940.pdf (equation 12)
    Double_t eta, dNdeta, par0, c, alpha, beta, a, par4, par5;
    a      = par[0];
    alpha  = par[1];
    beta   = par[2];
    c      = par[3];
    par4   = par[4];
    par5   = par[5];
    eta    = x_val[0];
    dNdeta = par4+par5*eta+c*TMath::Sqrt(1.0-1.0/TMath::Power(alpha*TMath::CosH(eta),2))/(1.0+TMath::Exp((fabs(eta)-beta)/a));
    return dNdeta;
}

#endif /* EP_DET_SIM_FUNC_H */