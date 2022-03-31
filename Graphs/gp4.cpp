#include "TCanvas.h"
#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TF1.h"

Double_t xmin = 0, xmax = 400;
double frac = 4.5;

double FindMax(double* x){
    int n = sizeof(x);
    double max = x[0];
    for(int i = 0; i<n;i++){
        if (max < x[i]) max = x[i];
    }
    return max;
}

double FindMin(double* x){
    int n = sizeof(x);
    double min = x[0];
    for(int i = 0; i<n;i++){
        if (min > x[i]) min = x[i];
    }
    return min;
}

Double_t Reta(double* x, double* par){
    return par[0] + par[1]*x[0];
}

void gp(){
    
    TCanvas *c1 = new TCanvas();    

    /*c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx();
    c1->SetGridy();*/

    c1->cd();

    TMultiGraph *MG = new TMultiGraph();

    TGraphErrors *gr = new TGraphErrors();

    /*gr->GetXaxis()->SetTitle("Photon Energy #bf{[eV]}");
    gr->GetYaxis()->SetTitle("Electron Kinetic Energy #bf{[eV]}");
    gr->GetXaxis()->SetLabelFont(132);
    gr->GetXaxis()->SetTitleFont(132);
    gr->GetYaxis()->SetLabelFont(132);
    gr->GetYaxis()->SetTitleFont(132);
    gr->GetYaxis()->SetLabelSize(0.035);
    gr->GetYaxis()->SetTitleSize(0.035);
    gr->GetXaxis()->SetLabelSize(0.035);
    gr->GetXaxis()->SetTitleSize(0.035);
    gr->GetYaxis()->SetTitleOffset(1.4);
    gr->GetXaxis()->SetMaxDigits(3);*/
    gr->SetMarkerStyle(kFullTriangleUp);
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(0.8);

    std::fstream file;
    file.open("EleMundMed.tsv", std::ios::in);

    double x1[2],ex = 0, ey = 0;

    int n = 0;

    while(1){
        file >> x1[0] >> x1[1];
    
        n = gr->GetN();

        gr->SetPoint(n, x1[1], x1[0]);
        gr->SetPointError(n, ex, ey);

        if(file.eof()) break;
    
    }

    file.close();

    gr->GetXaxis()->SetLimits(xmin, xmax);

    MG->Add(gr);

    c1->Update();

    TGraphErrors *gr2 = new TGraphErrors();

    /*gr2->GetXaxis()->SetTitle("Photon Energy #bf{[eV]}");
    gr2->GetYaxis()->SetTitle("Electron Kinetic Energy #bf{[eV]}");
    gr2->GetXaxis()->SetLabelFont(132);
    gr2->GetXaxis()->SetTitleFont(132);
    gr2->GetYaxis()->SetLabelFont(132);
    gr2->GetYaxis()->SetTitleFont(132);
    gr2->GetYaxis()->SetLabelSize(0.035);
    gr2->GetYaxis()->SetTitleSize(0.035);
    gr2->GetXaxis()->SetLabelSize(0.035);
    gr2->GetXaxis()->SetTitleSize(0.035);
    gr2->GetYaxis()->SetTitleOffset(1.4);
    gr2->GetXaxis()->SetMaxDigits(3);*/
    gr2->SetMarkerStyle(kFullTriangleUp);
    gr2->SetMarkerColor(kRed);
    gr2->SetMarkerSize(0.8);

    file.open("PhotFrac.tsv", std::ios::in);

    while(1){
        file >> x1[0] >> x1[1];
    
        n = gr2->GetN();

        gr2->SetPoint(n, x1[1], x1[0]);
        gr2->SetPointError(n, ex, ey);

        if(file.eof()) break;
    
    }

    file.close();

    gr2->GetXaxis()->SetLimits(xmin, xmax);

    float scale = frac*FindMax(gr->GetY())/(FindMax(gr2->GetY()));
    for (int i=0;i<gr2->GetN();i++) gr2->GetY()[i] *= scale;
   
    MG->Add(gr2);
    
    MG->GetXaxis()->SetLabelFont(132);
    MG->GetXaxis()->SetTitleFont(132);
    MG->GetXaxis()->SetLabelSize(0.035);
    MG->GetXaxis()->SetTitleSize(0.037);
    MG->GetXaxis()->SetTitleOffset(1.2);
    MG->GetXaxis()->SetTitle("Photon Kinetic Energy #bf{[eV]}");
    MG->GetYaxis()->SetLabelFont(132);
    MG->GetYaxis()->SetTitleFont(132);
    MG->GetYaxis()->SetLabelSize(0.035);
    MG->GetYaxis()->SetTitleSize(0.037);
    MG->GetYaxis()->SetTitleOffset(1);
    MG->GetYaxis()->SetTitle("Electron Kinetic Energy #bf{[eV]}");
    MG->Draw("AP");

    TGaxis* y2 = new TGaxis(MG->GetXaxis()->GetXmax(),0,MG->GetXaxis()->GetXmax(),
        MG->GetYaxis()->GetXmax(),0,MG->GetYaxis()->GetXmax()/scale,
        MG->GetYaxis()->GetNdivisions(),"+L");
    y2->SetLineColor(kRed);
    y2->SetLabelFont(132);
    y2->SetTitleFont(132);
    y2->SetLabelSize(0.035);
    y2->SetTitleSize(0.037);
    y2->SetTitleOffset(1);
    y2->SetMaxDigits(3);
    y2->SetTitle("Electron Emission Fraction");
    y2->SetLabelColor(kRed);
    y2->Draw();

    TF1 *ajt = new TF1("ajt",Reta,0,xmax,2);
    ajt->SetLineColor(kGreen);
    ajt->SetLineStyle(kDashDotted);

    ajt->SetParameters(-4.2,1);
    ajt->Draw("same");
}