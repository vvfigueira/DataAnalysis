#include "TCanvas.h"
#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMathText.h"
#include "TPad.h"
#include "TRatioPlot.h"
#include "TVirtualFFT.h"
#include "TSpline.h"
#include "TMultiGraph.h"

int npar = 4; // Número de parâmetros do ajuste
int intertype = 3; // Tipo de interpolação, 1 == Linear, 3 == Cúbica, 5 == Quíntica
int nmax = 48; // Número de pontos interpolados entre o ponto inicial e final do gráfico
int def = 5; // Defasagem para autocorrelação das elípses
int dx = 0, ex = 1, dy = 2, ey = 3; // Posições dos dados no arquivo
char* name = "D1.tsv"; // Nome do arquivo de dados
double confi = 0.05; // Intervalo de confiança

Double_t xmin = 3000, xmax = 6000; // Mínimo e máximo do eixo x no gráfico

Double_t Cnst(Double_t* x,Double_t* par){
    return confi;
}

Double_t Cnst2(Double_t* x,Double_t* par){
    return -confi;
}

Double_t Per(Double_t *x, Double_t *par){
    return par[0] + par[1]*TMath::Cos(x[0]*par[2]+par[3]);
}

void DFT(TGraphErrors* Graph, TGraphErrors* &FTGraph ){

    double P[Graph->GetN()];
    for(int i=0; i < Graph->GetN(); ++i) {
        P[i]=Graph->GetPointX(i);
    }

    TH1D* TempHis = new TH1D("h1","",Graph->GetN()-1,P);

    double x,y;
    for(int i=0; i < Graph->GetN()-1; ++i) {
        Graph->GetPoint(i, x, y);
        TempHis->Fill(x,y);
    }

    TH1* FTHist = 0;

    TVirtualFFT::SetTransform(0);
    FTHist = TempHis->FFT(FTHist, "MAG");

    double norm = (Graph->GetPointX(Graph->GetN()-1)-Graph->GetPointX(0));

    for(int k=0;k<Graph->GetN()/2 ;k++){
        
        FTGraph->SetPoint(k, FTHist->GetBinLowEdge(k+1)/norm, 
            FTHist->GetBinContent(k+1)/TMath::Sqrt(Graph->GetN()));
        FTGraph->SetPointError(k, 0, 0);
    
    }

    for(int k=1;k<Graph->GetN()/2 ;k++){
        
        FTGraph->SetPoint(Graph->GetN()/2 -1 + k, -FTHist->GetBinLowEdge(k+1)/norm, 
            FTHist->GetBinContent(k+1)/TMath::Sqrt(Graph->GetN()));
        FTGraph->SetPointError(Graph->GetN()/2 -1 + k, 0, 0);
    
    }

    FTGraph->GetXaxis()->SetLimits(-Graph->GetN()/(2*norm), Graph->GetN()/(2*norm));

    delete TempHis;
    delete FTHist;

}

void AutoCorr(TGraphErrors* Graph, TGraphErrors* &GCorDat){
    
    double sumsq = 0, prod = 0, autocorr = 0;
    
    for(int i=0;i<Graph->GetN();i++){
        sumsq = sumsq + pow(Graph->GetPointY(i)-Graph->GetMean(2),2);
    }

    for(int k=-Graph->GetN();k<Graph->GetN();k++){
        
        prod = 0;
        for(int j=0;j<Graph->GetN();j++){
            
            if(j-k >= 0 && j-k < Graph->GetN()){
                prod = prod + (Graph->GetPointY(j-k)-Graph->GetMean(2))*(Graph->GetPointY(j)-Graph->GetMean(2));
            }else{prod = prod + (-Graph->GetMean(2))*(Graph->GetPointY(j)-Graph->GetMean(2));}
        }

        autocorr = prod/sumsq;

        GCorDat->SetPoint(k+Graph->GetN(), k, autocorr);
        GCorDat->SetPointError(Graph->GetN()+k, 0, 0);
    
    }

}

void Format(TGraphErrors* &Graph, double Size = 0.035, int markstyle = kFullCircle, double Off = 0.65, 
    int Font = 132, double markSize = 0.6){
    Graph->GetXaxis()->SetLabelFont(Font);
    Graph->GetXaxis()->SetTitleFont(Font);
    Graph->GetYaxis()->SetLabelFont(Font);
    Graph->GetYaxis()->SetTitleFont(Font);
    Graph->GetYaxis()->SetLabelSize(Size);
    Graph->GetYaxis()->SetTitleSize(Size);
    Graph->GetXaxis()->SetLabelSize(Size);
    Graph->GetXaxis()->SetTitleSize(Size);
    Graph->GetYaxis()->SetTitleOffset(Off);
    Graph->SetMarkerStyle(markstyle);
    Graph->SetMarkerSize(markSize);
}

void Format(TMultiGraph* &Graph, double Size = 0.035, double Off = 0.65, 
    int Font = 132, double markSize = 0.6){
    Graph->GetXaxis()->SetLabelFont(Font);
    Graph->GetXaxis()->SetTitleFont(Font);
    Graph->GetYaxis()->SetLabelFont(Font);
    Graph->GetYaxis()->SetTitleFont(Font);
    Graph->GetYaxis()->SetLabelSize(Size);
    Graph->GetYaxis()->SetTitleSize(Size);
    Graph->GetXaxis()->SetLabelSize(Size);
    Graph->GetXaxis()->SetTitleSize(Size);
    Graph->GetYaxis()->SetTitleOffset(Off);

}

void gp(){

    std::fstream file;
    file.open(name, std::ios::in);

    if (!file.is_open()) {
        std::cerr << "ERRO: O arquivo não pode ser aberto!\n";
        return;
    }

    Double_t x1[4], par[npar];
    int n = 0;

    TCanvas *Princ = new TCanvas();

    Princ->SetTickx();
    Princ->SetTicky();
    Princ->SetGridx();
    Princ->SetGridy();

    TPad * GPad = new TPad("Graph","",0.,0.3,1,0.95);

    GPad->SetFillStyle(4000);
    GPad->SetMargin(0.08,0.08,0.01,0.055);
    GPad->SetGridx();
    GPad->SetGridy();
    GPad->Draw();

    TPad * RPad = new TPad("Rezid","",0.,0.,1,0.32);
    
    RPad->SetFillStyle(4000);
    RPad->SetMargin(0.08,0.08,0.21,0.01);
    RPad->SetGridx();
    RPad->SetGridy();
    RPad->Draw();

    TMultiGraph *MGraph = new TMultiGraph();

    TGraphErrors *Graph = new TGraphErrors();

    while(1){
        file >> x1[0] >> x1[1]>> x1[2] >> x1[3];
    
        n = Graph->GetN();

        Graph->SetPoint(n, x1[dx], x1[dy]);
        Graph->SetPointError(n, x1[ex], x1[ey]);

        if(file.eof()) break;
    
    }

    file.close();

    //Format(Graph, 0.05, kOpenCircle);
    Graph->SetMarkerStyle(kFullCircle);
    Graph->SetMarkerColor(kBlack);
    Graph->SetMarkerSize(0.7);
    //Graph->GetYaxis()->SetTitle("Res#acute{i}duos");
    //Graph->GetXaxis()->SetLimits(xmin, xmax);
    //Graph->GetYaxis()->SetMaxDigits(2);
    //GPad->cd();
    //Graph->Draw("AP");
    MGraph->Add(Graph);

    TF1 *Ajt = new TF1("Ajt", Per, xmin, xmax, npar);

    Ajt->SetLineColor(kRed);
    Ajt->SetParameter(0, 0);
    Ajt->SetParameter(1, 0.01);
    Ajt->SetParameter(2, TMath::TwoPi()/365.);
    Ajt->SetParameter(3, -152*TMath::TwoPi()/365.);
    Ajt->SetNpx(1.e5);

    Graph->Fit("Ajt","","",xmin,xmax);

    TGraphErrors *GInter = new TGraphErrors();

    TSpline* Spline = 0;

    double step = (Graph->GetPointX(Graph->GetN()-1)-Graph->GetPointX(0))/(nmax+1);

    switch (intertype){

    case 1 :

        int data;

        GInter->SetPointError(0,0,0);
        GInter->SetPoint(0, Graph->GetPointX(0), Graph->GetPointY(0));

        for(int i = 1; i<nmax+1; i++){
            data = 0;
            GInter->SetPointError(i,0,0);
            for(int j = 0; j < (Graph->GetN() - 1); j++){
                if( ((Graph->GetPointX(j) <= (Graph->GetPointX(0)+ step*i)) &&
                    ((Graph->GetPointX(0)+ step*i) < Graph->GetPointX(j+1)))){ data = j; j=Graph->GetN();}
            }
            GInter->SetPoint(i, Graph->GetPointX(0)+ step*i, 
                Graph->GetPointY(data)+((Graph->GetPointY(data + 1)-Graph->GetPointY(data))/(Graph->GetPointX(data + 1) - Graph->GetPointX(data)))*(Graph->GetPointX(0)+ step*i - Graph->GetPointX(data)));
        }

        GInter->SetPointError(nmax + 1,0,0);
        GInter->SetPoint(nmax + 1, Graph->GetPointX(Graph->GetN() -1), Graph->GetPointY(Graph->GetN() -1));
        break;

    case 3 :
        Spline = new TSpline3("Interpol", Graph);
        for(int i = 0; i<nmax+2; i++){
        
            GInter->SetPointError(i,0,0);
            GInter->SetPoint(i, Graph->GetPointX(0)+ step*i, 
                Spline->Eval(Graph->GetPointX(0)+ step*i));
        }
        break;

    case 5 :
        Spline = new TSpline5("Interpol", Graph);
            for(int i = 0; i<nmax+2; i++){
            
                GInter->SetPointError(i,0,0);
                GInter->SetPoint(i, Graph->GetPointX(0)+ step*i, 
                    Spline->Eval(Graph->GetPointX(0)+ step*i));
            }
        break;

    default:
        std::cerr << "ERRO: Selecione um método de Interpolação válido\n\tMétodos: {1, 3, 5}\n";
        return;
        break;

    }

    //Format(GInter, 0.05, kFullCircle);
    GInter->SetMarkerStyle(kFullTriangleUp);
    GInter->SetMarkerColor(kBlue);
    GInter->SetMarkerSize(0.7);
    //GInter->GetYaxis()->SetTitle("Res#acute{i}duos");
    //GInter->GetXaxis()->SetLimits(xmin, xmax);
    //GInter->GetYaxis()->SetMaxDigits(2);
    //GPad->cd();
    //GInter->Draw("same AP");
    MGraph->Add(GInter);
    GPad->cd();
    Format(MGraph, 0.055);
    MGraph->GetYaxis()->SetTitle("Residual");
    MGraph->GetXaxis()->SetLimits(xmin, xmax);
    MGraph->GetYaxis()->SetMaxDigits(2);
    MGraph->Draw("AP");

    Ajt->GetParameters(par);

    TMultiGraph * GRed = new TMultiGraph();

    TGraphErrors* Rezid = new TGraphErrors();    

    for(int i = 0; i< Graph->GetN(); i++){
        
        double x[1] = {Graph->GetPointX(i)};
        Rezid->SetPoint(i, Graph->GetPointX(i), -(Per(x,par) - Graph->GetPointY(i)));
        Rezid->SetPointError(i, Graph->GetErrorX(i), Graph->GetErrorY(i));

    }


    Rezid->SetMarkerStyle(kFullCircle);
    Rezid->SetMarkerColor(kBlack);
    Rezid->SetMarkerSize(0.7);
    /*Format(Rezid, 0.1, kOpenCircle, 1.3);
    Rezid->GetXaxis()->SetTitle("#bf{[Dias]}");
    Rezid->GetYaxis()->SetNdivisions(8);
    Rezid->GetXaxis()->SetLimits(xmin, xmax);
    RPad->cd();
    Rezid->Draw("AP");*/
    GRed->Add(Rezid);

    TGraphErrors* RInt = new TGraphErrors();    

    for(int i = 0; i< GInter->GetN(); i++){
        
        double x[1] = {GInter->GetPointX(i)};
        RInt->SetPoint(i, GInter->GetPointX(i), -(Per(x,par) - GInter->GetPointY(i)));
        RInt->SetPointError(i, GInter->GetErrorX(i), GInter->GetErrorY(i));

    }


    RInt->SetMarkerStyle(kFullTriangleUp);
    RInt->SetMarkerColor(kBlue);
    RInt->SetMarkerSize(0.7);
    /*Format(RInt, 0.1, kOpenCircle, 1.3);
    RInt->GetXaxis()->SetTitle("#bf{[Dias]}");
    RInt->GetYaxis()->SetNdivisions(8);
    RInt->GetXaxis()->SetLimits(xmin, xmax);
    RPad->cd();
    RInt->Draw("AP");*/
    GRed->Add(RInt);

    Format(GRed, 0.1);
    GRed->GetXaxis()->SetTitle("#bf{[Days]}");
    GRed->GetXaxis()->SetLimits(xmin, xmax);
    //GRed->GetYaxis()->SetMaxDigits(2);
    RPad->cd();
    GRed->Draw("AP");

    TCanvas* ECanv = new TCanvas();

    ECanv->SetTickx();
    ECanv->SetTicky();
    ECanv->SetGridx();
    ECanv->SetGridy();

    TMultiGraph* MEGraph = new TMultiGraph();

    TGraphErrors* IntEGraph = new TGraphErrors();

    for(int i = 0; i< GInter->GetN(); i++){
        
        if( (i+def) < GInter->GetN()){
            IntEGraph->SetPoint(i, GInter->GetPointY(i), GInter->GetPointY(i+def));
        }else{
            IntEGraph->SetPoint(i, GInter->GetPointY(i), GInter->GetPointY(i+def-GInter->GetN()));
        }
        IntEGraph->SetPointError(i, GInter->GetErrorX(i), GInter->GetErrorY(i));

    }

    IntEGraph->SetMarkerStyle(kFullTriangleUp);
    IntEGraph->SetMarkerColor(kBlue);
    IntEGraph->SetMarkerSize(0.7);
    /*Format(IntEGraph,0.035,20,1);
    IntEGraph->GetYaxis()->SetTitle("Magnitude das frequ#hat{e}ncias");
    IntEGraph->GetXaxis()->SetTitle("Frequ#hat{e}ncia #bf{[Dias^{-1}]}");
    IntEGraph->GetXaxis()->SetMaxDigits(3);
    IntEGraph->GetYaxis()->SetMaxDigits(3);
    IntEGraph->GetYaxis()->SetRangeUser(0, 0.035);
    FTCanvas->cd();
    IntEGraph->Draw("AP");*/
    MEGraph->Add(IntEGraph);

    TGraphErrors* aux = new TGraphErrors();

    for(int i = 0; i< GInter->GetN(); i++){
        
        double x[1] = {GInter->GetPointX(i)};
        aux->SetPoint(i, GInter->GetPointX(i), Per(x,par));
        aux->SetPointError(i, GInter->GetErrorX(i), GInter->GetErrorY(i));

    }

    TGraphErrors* AjEGraph = new TGraphErrors();

    for(int i = 0; i< aux->GetN(); i++){

        if( (i+def) < aux->GetN()){
            AjEGraph->SetPoint(i, aux->GetPointY(i), aux->GetPointY(i+def));
        }else{
            AjEGraph->SetPoint(i, aux->GetPointY(i), aux->GetPointY(i+def - aux->GetN()));
        }
        AjEGraph->SetPointError(i, aux->GetErrorX(i), aux->GetErrorY(i));
    }

    AjEGraph->SetMarkerStyle(kFullTriangleUp);
    AjEGraph->SetMarkerColor(kRed);
    AjEGraph->SetMarkerSize(0.7);
    /*Format(IntEGraph,0.035,20,1);
    IntEGraph->GetYaxis()->SetTitle("Magnitude das frequ#hat{e}ncias");
    IntEGraph->GetXaxis()->SetTitle("Frequ#hat{e}ncia #bf{[Dias^{-1}]}");
    IntEGraph->GetXaxis()->SetMaxDigits(3);
    IntEGraph->GetYaxis()->SetMaxDigits(3);
    IntEGraph->GetYaxis()->SetRangeUser(0, 0.035);
    FTCanvas->cd();
    IntEGraph->Draw("AP");*/
    MEGraph->Add(AjEGraph);

    Format(MEGraph, 0.035, 1, 132, 1);

    MEGraph->GetYaxis()->CenterTitle();
    MEGraph->GetYaxis()->SetTitle("#bf{N'}");
    MEGraph->GetXaxis()->CenterTitle();
    MEGraph->GetXaxis()->SetTitle("#bf{N}");
    MEGraph->GetYaxis()->SetMaxDigits(3);
    MEGraph->Draw("AP");

    TCanvas* FTCanvas = new TCanvas();

    FTCanvas->SetTickx();
    FTCanvas->SetTicky();
    FTCanvas->SetGridx();
    FTCanvas->SetGridy();

    TMultiGraph* FTG = new TMultiGraph();

    TGraphErrors* FTGraph = new TGraphErrors();

    DFT(Graph,FTGraph);

    FTGraph->SetMarkerStyle(kFullCircle);
    FTGraph->SetMarkerColor(kBlack);
    FTGraph->SetMarkerSize(0.7);
    /*Format(FTGraph,0.035,20,1);
    FTGraph->GetYaxis()->SetTitle("Magnitude das frequ#hat{e}ncias");
    FTGraph->GetXaxis()->SetTitle("Frequ#hat{e}ncia #bf{[Dias^{-1}]}");
    FTGraph->GetXaxis()->SetMaxDigits(3);
    FTGraph->GetYaxis()->SetMaxDigits(3);
    FTGraph->GetYaxis()->SetRangeUser(0, 0.035);
    FTCanvas->cd();
    FTGraph->Draw("AP");*/
    FTG->Add(FTGraph);

    TGraphErrors* FTGInt = new TGraphErrors();

    DFT(GInter,FTGInt);

    FTGInt->SetMarkerStyle(kFullTriangleUp);
    FTGInt->SetMarkerColor(kBlue);
    FTGInt->SetMarkerSize(0.7);
    /*Format(FTGInt,0.035,20,1);
    FTGInt->GetYaxis()->SetTitle("Magnitude das frequ#hat{e}ncias");
    FTGInt->GetXaxis()->SetTitle("Frequ#hat{e}ncia #bf{[Dias^{-1}]}");
    FTGInt->GetXaxis()->SetMaxDigits(3);
    FTGInt->GetYaxis()->SetMaxDigits(3);
    FTGInt->GetYaxis()->SetRangeUser(0, 0.035);
    FTCanvas->cd();
    FTGInt->Draw("AP");*/
    FTG->Add(FTGInt);

    Format(FTG, 0.035, 1, 132, 1);
    FTG->GetYaxis()->SetTitle("Frequency Magnitude");
    FTG->GetXaxis()->SetTitle("Frequency #bf{[Days^{-1}]}");
    FTG->GetYaxis()->SetMaxDigits(3);
    FTG->Draw("AP");

    TCanvas *CorDat = new TCanvas();

    CorDat->SetTickx();
    CorDat->SetTicky();
    CorDat->SetGridx();
    CorDat->SetGridy();

    TMultiGraph *GCD = new TMultiGraph(); 
    
    TGraphErrors *GCorDat = new TGraphErrors();

    AutoCorr(Graph, GCorDat);

    GCorDat->SetMarkerStyle(kFullCircle);
    GCorDat->SetMarkerColor(kBlack);
    GCorDat->SetMarkerSize(0.7);
    /*Format(GCorDat,0.035,20,1);
    GCorDat->GetYaxis()->SetTitle("Auto-Correlac#tilde{a}o");
    GCorDat->GetXaxis()->SetTitle("Defasagem");
    GCorDat->GetXaxis()->SetLimits(-Graph->GetN(),Graph->GetN());
    CorDat->cd();
    GCorDat->Draw("AL");*/
    GCD->Add(GCorDat);

    TGraphErrors *GCorInt = new TGraphErrors();

    AutoCorr(GInter, GCorInt);

    GCorInt->SetMarkerStyle(kFullTriangleUp);
    GCorInt->SetMarkerColor(kBlue);
    GCorInt->SetLineColor(kBlue);
    GCorInt->SetMarkerSize(0.7);
    /*Format(GCorInt,0.035,20,1);
    GCorInt->GetYaxis()->SetTitle("Auto-Correlac#tilde{a}o");
    GCorInt->GetXaxis()->SetTitle("Defasagem");
    GCorInt->GetXaxis()->SetLimits(-Graph->GetN(),Graph->GetN());
    CorDat->cd();
    GCorInt->Draw("AL");*/
    GCD->Add(GCorInt);

    Format(GCD, 0.035, 1, 132, 1);
    GCD->GetYaxis()->SetTitle("Auto-Correlation");
    GCD->GetXaxis()->SetTitle("Shift");
    GCD->GetYaxis()->SetMaxDigits(3);
    GCD->Draw("AL");

    TCanvas *CorRez = new TCanvas();

    CorRez->SetTickx();
    CorRez->SetTicky();
    CorRez->SetGridx();
    CorRez->SetGridy();

    TGraphErrors *GCorRez = new TGraphErrors();

    AutoCorr(Rezid, GCorRez);
    
    Format(GCorRez,0.035,20,1);
    GCorRez->GetYaxis()->SetTitle("Auto-Correlation");
    GCorRez->GetXaxis()->SetTitle("Shift");
    GCorRez->GetXaxis()->SetLimits(-Rezid->GetN(), Rezid->GetN());
    CorRez->cd();
    GCorRez->Draw("AL");

    TF1 *Confi1 = new TF1("Confi1", Cnst, -10*Rezid->GetN(), 10*Rezid->GetN(), 0);

    Confi1->SetFillColor(kRed);
    CorRez->cd();
    Confi1->Draw("same");

    TF1 *Confi2 = new TF1("Confi2", Cnst2, -10*Rezid->GetN(), 10*Rezid->GetN(), 0);
    
    Confi2->SetFillColor(kRed);
    Confi2->Draw("same");

    /*Princ->SaveAs("teste2.pdf");
    ECanv->SaveAs();
    FTCanvas->SaveAs();
    CorDat->SaveAs();
    CorRez->SaveAs();*/
}