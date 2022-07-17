/// \file
/// \ingroup tutorial_fit
/// \notebook
/// Combined (simultaneous) fit of two histogram with separate functions
/// and some common parameters
///
/// See http://root.cern.ch/phpBB3//viewtopic.php?f=3&t=11740#p50908
/// for a modified version working with Fumili or GSLMultiFit
///
/// N.B. this macro must be compiled with ACliC
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"


int iparPrinc[3] = { 0,      // Rydberg Constant
                    1,    // exp common parameter
                    2
};

int iparNit[3] = { 0, // Rydberg Constant
                    2, // exp common parameter
                    1
};

int iparDif[3] = { 0, // Rydberg Constant
                2, // exp common parameter
                3
};

// Create the GlobalCHi2 structure

double Ryd(double *x, double *par){

    return par[0]*(TMath::Power(3-par[1], -2)-TMath::Power(x[0]-par[2], -2));
}

struct GlobalChi2 {
    GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction &f3) :
        fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[3];
        for (int i = 0; i < 3; ++i) p1[i] = par[iparPrinc[i] ];

        double p2[3];
        for (int i = 0; i < 3; ++i) p2[i] = par[iparNit[i] ];

        double p3[3];
        for (int i = 0; i < 3; ++i) p3[i] = par[iparDif[i] ];

        return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
    }

    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
    const  ROOT::Math::IMultiGenFunction * fChi2_3;
};

void combinedFit() {

    TH1D * hA = new TH1D("A", "Histo A", 100, 0, 100);
    TH1D * hB = new TH1D("B", "Histo B", 100, 0, 100);

    TGraphErrors *PGr = new TGraphErrors();
    TGraphErrors *NGr = new TGraphErrors();
    TGraphErrors *DGr = new TGraphErrors();

    PGr->SetPoint(0, 4, 0.000311385);
    PGr->SetPointError(0, 0, 1.06238e-6);
    PGr->SetPoint(1, 3, 0.000194211);
    PGr->SetPointError(1, 0, 6.77033e-7);

    PGr->SetMarkerStyle(kFullTriangleUp);
    PGr->SetMarkerColor(kBlack);
    PGr->SetMarkerSize(0.8);

    PGr->GetXaxis()->SetLabelFont(132);
    PGr->GetXaxis()->SetTitleFont(132);
    PGr->GetXaxis()->SetLabelSize(0.035);
    PGr->GetXaxis()->SetTitleSize(0.037);
    PGr->GetXaxis()->SetTitleOffset(1.2);
    PGr->GetXaxis()->SetTitle("n");
    PGr->GetYaxis()->SetLabelFont(132);
    PGr->GetYaxis()->SetTitleFont(132);
    PGr->GetYaxis()->SetLabelSize(0.035);
    PGr->GetYaxis()->SetTitleSize(0.037);
    PGr->GetYaxis()->SetTitleOffset(1);
    PGr->GetYaxis()->SetTitle("Inverse Wavelength #bf{[#AA^{-1}]}");
    PGr->GetYaxis()->SetMaxDigits(2);
    PGr->GetXaxis()->SetMaxDigits(2);

    NGr->SetPoint(0, 8, 0.000220115);
    NGr->SetPointError(0, 0, 7.44832e-7);
    NGr->SetPoint(1, 7, 0.000210447);
    NGr->SetPointError(1, 0, 7.37501e-7);
    NGr->SetPoint(2, 6, 0.000194211);
    NGr->SetPointError(2, 0, 6.77033e-7);
    NGr->SetPoint(3, 5, 0.000162124);
    NGr->SetPointError(3, 0, 5.52615e-7);

    NGr->SetMarkerStyle(kFullTriangleUp);
    NGr->SetMarkerColor(kBlack);
    NGr->SetMarkerSize(0.8);

    NGr->GetXaxis()->SetLabelFont(132);
    NGr->GetXaxis()->SetTitleFont(132);
    NGr->GetXaxis()->SetLabelSize(0.035);
    NGr->GetXaxis()->SetTitleSize(0.037);
    NGr->GetXaxis()->SetTitleOffset(1.2);
    NGr->GetXaxis()->SetTitle("n");
    NGr->GetYaxis()->SetLabelFont(132);
    NGr->GetYaxis()->SetTitleFont(132);
    NGr->GetYaxis()->SetLabelSize(0.035);
    NGr->GetYaxis()->SetTitleSize(0.037);
    NGr->GetYaxis()->SetTitleOffset(1);
    NGr->GetYaxis()->SetTitle("Inverse Wavelength #bf{[#AA^{-1}]}");
    NGr->GetYaxis()->SetMaxDigits(2);
    NGr->GetXaxis()->SetMaxDigits(2);

    DGr->SetPoint(0, 7, 0.000222425);
    DGr->SetPointError(0, 0, 7.77688e-7);
    DGr->SetPoint(1, 6, 0.000214289);
    DGr->SetPointError(1, 0, 7.39337e-7);
    DGr->SetPoint(2, 5, 0.000200882);
    DGr->SetPointError(2, 0, 6.93952e-7);
    DGr->SetPoint(3, 4, 0.000175773);
    DGr->SetPointError(3, 0, 5.98879e-7);

    DGr->SetMarkerStyle(kFullTriangleUp);
    DGr->SetMarkerColor(kBlack);
    DGr->SetMarkerSize(0.8);

    DGr->GetXaxis()->SetLabelFont(132);
    DGr->GetXaxis()->SetTitleFont(132);
    DGr->GetXaxis()->SetLabelSize(0.035);
    DGr->GetXaxis()->SetTitleSize(0.037);
    DGr->GetXaxis()->SetTitleOffset(1.2);
    DGr->GetXaxis()->SetTitle("n");
    DGr->GetYaxis()->SetLabelFont(132);
    DGr->GetYaxis()->SetTitleFont(132);
    DGr->GetYaxis()->SetLabelSize(0.035);
    DGr->GetYaxis()->SetTitleSize(0.037);
    DGr->GetYaxis()->SetTitleOffset(1);
    DGr->GetYaxis()->SetTitle("Inverse Wavelength #bf{[#AA^{-1}]}");
    DGr->GetYaxis()->SetMaxDigits(2);
    DGr->GetXaxis()->SetMaxDigits(2);

    TF1 * fPrin = new TF1("fPrin", Ryd, 0, 10, 3);
    fPrin->SetParameters(1.096e-3, 1.37, 0.87);

    TF1 * fNit = new TF1("fNit", Ryd, 0, 10, 3);
    fNit->SetParameters(1.096e-3, 0.87, 1.37);

    TF1 * fDif = new TF1("fDif", Ryd, 0, 10, 3);
    fDif->SetParameters(1.096e-3, 0.87, 0.01);

    // perform now global fit

    // TF1 * fSB = new TF1("fSB","expo + gaus(2)",0,100);

    ROOT::Math::WrappedMultiTF1 wfPrin(*fPrin,1);
    ROOT::Math::WrappedMultiTF1 wfNit(*fNit,1);
    ROOT::Math::WrappedMultiTF1 wfDif(*fDif,1);

    ROOT::Fit::DataOptions opt;

    ROOT::Fit::DataRange rangeP;
    rangeP.SetRange(0,10);
    ROOT::Fit::BinData dataP(opt, rangeP);
    ROOT::Fit::FillData(dataP, PGr);

    ROOT::Fit::DataRange rangeN;
    rangeN.SetRange(0,10);
    ROOT::Fit::BinData dataN(opt, rangeN);
    ROOT::Fit::FillData(dataN, NGr);

    ROOT::Fit::DataRange rangeD;
    rangeD.SetRange(0, 10);
    ROOT::Fit::BinData dataD(opt, rangeD);
    ROOT::Fit::FillData(dataD, DGr);

    ROOT::Fit::Chi2Function chi2_P(dataP, wfPrin);
    ROOT::Fit::Chi2Function chi2_N(dataN, wfNit);
    ROOT::Fit::Chi2Function chi2_D(dataD, wfDif);

    GlobalChi2 globalChi2(chi2_P, chi2_N, chi2_D);  

    ROOT::Fit::Fitter fitter;

    const int Npar = 4;
    double par0[Npar] = { 1.096e-3, 1.37, 0.87, 0.01};

    // create before the parameter settings in order to fix or set range on them
    fitter.Config().SetParamsSettings(4,par0);
    // fix 5-th parameter
    // fitter.Config().ParSettings(0).Fix();
    // fitter.Config().ParSettings(1).Fix();
    // fitter.Config().ParSettings(2).Fix();
    // set limits on the third and 4-th parameter
    // fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
    // fitter.Config().ParSettings(3).SetLimits(0,10000);
    // fitter.Config().ParSettings(3).SetStepSize(5);

    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");

    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(4,globalChi2,0,dataP.Size()+dataN.Size()+dataD.Size(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

    TCanvas * c1 = new TCanvas("","",
                                10,10,700,700);
    c1->Divide(3,1);
    c1->cd(1);
    gStyle->SetOptFit(1111);

    fPrin->SetFitResult( result, iparPrinc);
    fPrin->SetRange(rangeP().first, rangeP().second);
    fPrin->SetLineColor(kRed);
    PGr->GetListOfFunctions()->Add(fPrin);
    PGr->SetStats(true);
    PGr->Draw("AP");

    c1->cd(2);
    fNit->SetFitResult( result, iparNit);
    fNit->SetRange(rangeN().first, rangeN().second);
    fNit->SetLineColor(kRed);
    NGr->GetListOfFunctions()->Add(fNit);
    NGr->SetStats(0);
    NGr->Draw("AP");

    c1->cd(3);
    fDif->SetFitResult( result, iparDif);
    fDif->SetRange(rangeD().first, rangeD().second);
    fDif->SetLineColor(kRed);
    DGr->GetListOfFunctions()->Add(fDif);
    DGr->SetStats(0);
    DGr->Draw("AP");


}
