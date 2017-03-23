#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    printf("./world_gt <opt>\n");
    return 1;
  }

  //base (blank)
  TH1D * h0 = new TH1D("h0", "", 1, -5, 5);
  h0->SetStats(0);
  h0->SetTitle("");
  h0->GetXaxis()->SetTitle("");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleFont(22);
  h0->GetXaxis()->SetTitleSize(0.055);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetTitle("");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleFont(22);
  h0->GetYaxis()->SetTitleSize(0.055);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetLabelOffset(999);
  //h0->GetYaxis()->SetNdivisions(15, 0, 0, kTRUE);

  TCanvas * c0 = new TCanvas("c0", "", 1200, 700);
  //c0->SetBottomMargin(0.15);
  //c0->SetLeftMargin(0.15);

  
  const int opt = atoi(argv[1]);

  if (opt == 1){
    gStyle->SetEndErrorSize(6);
    c0->Clear();
    h0->GetXaxis()->SetLimits(-0.95, 2.5);
    h0->GetYaxis()->SetRangeUser(-33, -4.5);
    h0->GetYaxis()->SetNdivisions(0, 0, 0, kTRUE);
    h0->DrawClone("AXIS");
    double iy = -10.0;
    double Xtext = 1.7;
    double tsize = 0.025;
    double msize = 1.3;
    double lwidth = 1.;
    TLatex * t0;
    t0 = new TLatex(0.65, -8.0, "#font[32]{g_{T}^{u}}");
    t0->DrawClone("same");
    t0 = new TLatex(-0.2, -8.0, "#font[32]{g_{T}^{d}}");
    t0->DrawClone("same");


    TGraphErrors * e0 = new TGraphErrors(2);
    e0->SetMarkerStyle(8);
    e0->SetMarkerColor(8);
    e0->SetMarkerSize(msize);
    e0->SetLineColor(8);
    e0->SetLineWidth(lwidth);
    //DSE PRD91(2015)074004 Pitschmann et al.
    iy -= 1.0;
    e0->SetPoint(0, 0.55, iy);//u
    e0->SetPointError(0, 0.08, 0);
    e0->SetPoint(1, -0.11, iy);//d
    e0->SetPointError(1, 0.02, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Pitschmann et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e0->DrawClone("pesame");
    t0->DrawClone("same");
    //DSE PRD88(2013)074036 Yamanaka et al.
    iy -= 1.0;
    e0->SetPoint(0, 0.8, iy);//u
    e0->SetPoint(1, -0.2, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Yamanaka et al. (2013)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e0->DrawClone("pxsame");
    t0->DrawClone("same");

    TGraphErrors * e1 = new TGraphErrors(2);
    e1->SetMarkerStyle(8);
    e1->SetMarkerColor(4);
    e1->SetMarkerSize(msize);
    e1->SetLineColor(4);
    e1->SetLineWidth(lwidth);
    //Lattice PRD94(2016)054508 Bhattacharya et al.
    iy -= 1.0;
    e1->SetPoint(0, 0.792, iy);//u
    e1->SetPointError(0, 0.042, 0);
    e1->SetPoint(1, -0.194, iy);//d
    e1->SetPointError(1, 0.014, iy);
    t0 = new TLatex(Xtext, iy, "#font[22]{Bhattacharya et al. (2016)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PRD92(2015)114513 Abdel-Rehim et al.
    iy -= 1.0;
    e1->SetPoint(0, 0.791, iy);//u
    e1->SetPointError(0, 0.053, 0);
    e1->SetPoint(1, -0.236, iy);//d
    e1->SetPointError(1, 0.033, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Abdel-Rehim et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PLB627(2005)113 G\"ockeler et al.
    iy -= 1.0;
    e1->SetPoint(0, 0.857, iy);//u
    e1->SetPointError(0, 0.013, 0);
    e1->SetPoint(1, -0.212, iy);//d
    e1->SetPointError(1, 0.005, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Gockeler et al. (2005)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e1->DrawClone("pesame");
    t0->DrawClone("same");

    TGraphErrors * e2 = new TGraphErrors(2);
    e2->SetMarkerStyle(8);
    e2->SetMarkerColor(6);
    e2->SetMarkerSize(msize);
    e2->SetLineColor(6);
    e2->SetLineWidth(lwidth);
    //Model PLB659(2008)214 Cloet et al.
    iy -= 1.0;
    e2->SetPoint(0, 1.04, iy);//u
    e2->SetPoint(1, -0.24, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Cloet et al. (2008)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
    //Model PLB653(2007)398 Wakamatsu
    iy -= 1.0;
    e2->SetPoint(0, 0.90, iy);//u
    e2->SetPoint(1, -0.19, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Wakamatsu (2007)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
    //Model PRD72(2005)094029 Pasquini et al.
    iy -= 1.0;
    e2->SetPoint(0, 0.97, iy);//u
    e2->SetPoint(1, -0.24, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Pasquini et al. (2005)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
    //Model PRL87(2001)242001 Gamberg and Goldstein
    iy -= 1.0;
    e2->SetPoint(0, 0.795, iy);//u
    e2->SetPointError(0, 0.294, 0);
    e2->SetPoint(1, -0.155, iy);//d
    e2->SetPointError(1, 0.205, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Gamberg, Goldstein (2001)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pesame");
    t0->DrawClone("same");
    //Model PRD64(2001)034013 Schweitzer et al.
    iy -= 1.0;
    e2->SetPoint(0, (0.63+1.06)/2.0, iy);//u
    e2->SetPoint(1, (0.63-1.06)/2.0, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Schweitzer et al. (2001)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
    //Model JPG24(1998)L71 Ma and Schmidt
    iy -= 1.0;
    e2->SetPoint(0, (0.84+1.09)/2.0, iy);//u
    e2->SetPointError(0, (1.09-0.84)/2.0, 0);
    e2->SetPoint(1, -(0.23+0.51)/2.0, iy);//d
    e2->SetPointError(1, (0.51-0.23)/2.0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Ma, Schmidt (1998)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pesame");
    t0->DrawClone("same");
    //Model PLB390(1997)287 Barone et al.
    iy -= 1.0;
    e2->SetPoint(0, 0.969, iy);//u
    e2->SetPoint(1, -0.25, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Barone et al. (1997)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
    //Model PLB407(1997)331 Schmidt and Soffer
    iy -= 1.0;
    e2->SetPoint(0, 7.0/6.0, iy);//u
    e2->SetPoint(1, -7.0/24.0, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Schmidt, Soffer (1997)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
    //Model PRD52(1995)2960  He and Ji
    if (false){
      iy -= 1.0;
      e2->SetPoint(0, 1.0, iy-0.1);//u
      e2->SetPointError(0, 0.5, 0);
      e2->SetPoint(1, 0.0, iy+0.1);//d
      e2->SetPointError(1, 0.5, 0);
      t0 = new TLatex(Xtext, iy, "#font[22]{                         He, Ji (1995)}");
      t0->SetTextAlign(12);
      t0->SetTextSize(tsize);
      e2->DrawClone("pesame");
      t0->DrawClone("same");
    }
    //Model PRD54(1996)6897 He and Ji
    iy -= 1.0;
    e2->SetPoint(0, 1.33, iy);//u
    e2->SetPointError(0, 0.53, 0);
    e2->SetPoint(1, 0.04, iy);
    e2->SetPointError(1, 0.02, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{                         He, Ji (1996)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pesame");
    t0->DrawClone("same");
    //Model PLB377(1996)577 Kim et al.
    iy -= 1.0;
    e2->SetPoint(0, 1.12, iy);//u
    e2->SetPoint(1, -0.41, iy);//d
    t0 = new TLatex(Xtext, iy, "#font[22]{Kim et al. (1996)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e2->DrawClone("pxsame");
    t0->DrawClone("same");
 
    
    TGraphAsymmErrors * e3 = new TGraphAsymmErrors(2);
    e3->SetMarkerStyle(8);
    e3->SetMarkerColor(1);
    e3->SetMarkerSize(msize);
    e3->SetLineColor(1);
    e3->SetLineWidth(lwidth);
    //Fit PRD93(2016)014009 Kang et al.
    iy -= 1.0;
    e3->SetPoint(0, 0.39, iy);//u
    e3->SetPointError(0, 0.20, 0.16, 0, 0);
    e3->SetPoint(1, -0.22, iy);//d
    e3->SetPointError(1, 0.10, 0.31, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Kang et al. (2016)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit JHEP05(2015)123 Radici et al.
    iy -= 1.0;
    e3->SetPoint(0, 0.39, iy);//u
    e3->SetPointError(0, 0.15, 0.15, 0, 0);
    e3->SetPoint(1, -0.41, iy);//d
    e3->SetPointError(1, 0.52, 0.52, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Radici et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit 1401.0483 Goldstein et al.
    iy -= 1.0;
    e3->SetPoint(0, 0.860, iy);//u
    e3->SetPointError(0, 0.248, 0.248, 0, 0);
    e3->SetPoint(1, -0.119, iy);//d
    e3->SetPointError(1, 0.060, 0.060, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Goldstein et al. (2014)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit PRD87(2013)094019 Anselmino et al.
    iy -= 1.0;
    e3->SetPoint(0, 0.39, iy);//u
    e3->SetPointError(0, 0.12, 0.18, 0, 0);
    e3->SetPoint(1, -0.25, iy);//d
    e3->SetPointError(1, 0.10, 0.30, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Anselmino et al. (2013)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit PLB767(2017)91 Ye et al.
    iy -= 1.0;
    e3->SetPoint(0, 0.413, iy);//u
    e3->SetPointError(0, 0.133, 0.133, 0, 0);
    e3->SetPoint(1, -0.229, iy);//d
    e3->SetPointError(1, 0.094, 0.094, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Ye et al. (2017)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e3->DrawClone("pesame");
    t0->DrawClone("same");


    TGraphErrors * e4 = new TGraphErrors(2);
    e4->SetMarkerStyle(8);
    e4->SetMarkerColor(2);
    e4->SetMarkerSize(msize);
    e4->SetLineColor(2);
    e4->SetLineWidth(lwidth);
    //SoLID PLB767(2017)91 Ye et al. 
    iy -= 1.0;
    e4->SetPoint(0, 0.413, iy);//u
    e4->SetPointError(0, 0.018, 0);
    e4->SetPoint(1, -0.229, iy);//d
    e4->SetPointError(1, 0.008, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{JLab12 SoLID}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    e4->DrawClone("pesame");
    t0->DrawClone("same");


    TLegend * leg = new TLegend(0.63, 0.75, 0.9, 0.9);
    leg->AddEntry(e0, "#font[22]{Dyson-Schwinger equation}", "p");
    leg->AddEntry(e1, "#font[22]{Lattice QCD}", "p");
    leg->AddEntry(e2, "#font[22]{Models}", "p");
    leg->AddEntry(e3, "#font[22]{Phenomenology}", "p");
    leg->AddEntry(e4, "#font[22]{Future experiment}", "p");


    
    leg->DrawClone("same");
    c0->Print("gt0.pdf");
    
  }

  if (opt == 2){
    c0->Clear();
    c0->SetCanvasSize(1200, 800);
    gStyle->SetEndErrorSize(6);
    h0->GetXaxis()->SetLimits(0.1, 2.3);
    h0->GetYaxis()->SetRangeUser(-26, -6.0);
    h0->GetXaxis()->SetNdivisions(5, 5, 0, kTRUE);
    h0->GetYaxis()->SetNdivisions(0, 0, 0, kTRUE);
    h0->DrawClone("AXIS");
    double iy = -10.0;
    double Xtext = 1.6;
    double tsize = 0.027;
    double msize = 1.3;
    double lwidth = 1.;
    TLatex * t0;
    t0 = new TLatex(0.75, -8.0, "#font[32]{g_{T}^{(1)} = g_{T}^{u} - g_{T}^{d}}");
    t0->SetTextAlign(22);
    t0->DrawClone("same");

    TGraphErrors * a0 = new TGraphErrors(1);
    a0->SetMarkerStyle(8);
    a0->SetMarkerColor(8);
    a0->SetMarkerSize(msize);
    a0->SetLineColor(8);
    a0->SetLineWidth(lwidth);
    //DSE PRD91(2015)074004 Pitschmann et al.
    iy -= 1.0;
    a0->SetPoint(0, 0.66, iy);
    a0->SetPointError(0, 0.10, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Pitschmann et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a0->DrawClone("pesame");
    t0->DrawClone("same");
    //DSE PRD88(2013)074036 Yamanaka et al.
    iy -= 1.0;
    a0->SetPoint(0, 1.0, iy);
    t0 = new TLatex(Xtext, iy, "#font[22]{Yamanaka et al. (2013)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a0->DrawClone("pxsame");
    t0->DrawClone("same");

    TGraphErrors * a1 = new TGraphErrors(1);
    a1->SetMarkerStyle(8);
    a1->SetMarkerColor(4);
    a1->SetMarkerSize(msize);
    a1->SetLineColor(4);
    a1->SetLineWidth(lwidth);
    //Lattice PRD94(2016)054508 Bhattacharya et al.
    iy -= 1.0;
    a1->SetPoint(0, 0.987, iy);
    a1->SetPointError(0, 0.051, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Bhattacharya et al. (2016)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PRD92(2015)114513 Abdel-Rehim et al.
    iy -= 1.0;
    a1->SetPoint(0, 1.027, iy);
    a1->SetPointError(0, 0.062, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Abdel-Rehim et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PRD91(2015)054501 Bali et al.
    iy -= 1.0;
    a1->SetPoint(0, 1.005, iy);
    a1->SetPointError(0, 0.017, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Bali et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PRD86(2012)114509 Green et al.
    iy -= 1.0;
    a1->SetPoint(0, 1.038, iy);
    a1->SetPointError(0, 0.0163, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Green et al. (2012)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PRD82(2010)014501 Aoki et al.
    iy -= 1.0;
    a1->SetPoint(0, 0.990, iy);
    a1->SetPointError(0, 0.035, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Aoki et al. (2010)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a1->DrawClone("pesame");
    t0->DrawClone("same");
    //Lattice PLB627(2005)113 Gockeler et al.
    iy -= 1.0;
    a1->SetPoint(0, 1.068, iy);
    a1->SetPointError(0, 0.016, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Gockeler et al. (2005)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a1->DrawClone("pesame");
    t0->DrawClone("same");

    TGraphAsymmErrors * a3 = new TGraphAsymmErrors(1);
    a3->SetMarkerStyle(8);
    a3->SetMarkerColor(1);
    a3->SetMarkerSize(msize);
    a3->SetLineColor(1);
    a3->SetLineWidth(lwidth);
    //Fit PRD93(2016)014009 Kang et al.
    iy -= 1.0;
    a3->SetPoint(0, 0.61, iy);
    a3->SetPointError(0, 0.25, 0.15, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Kang et al. (2016)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit JHEP05(2015)123 Radici et al.
    iy -= 1.0;
    a3->SetPoint(0, 0.81, iy);
    a3->SetPointError(0, 0.44, 0.44, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Radici et al. (2015)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit 1401.0483 Goldstein et al.
    iy -= 1.0;
    a3->SetPoint(0, 0.979, iy);
    a3->SetPointError(0, 0.25, 0.25, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Goldstein et al. (2014)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit PRD87(2013)094019 Anselmino et al.
    iy -= 1.0;
    a3->SetPoint(0, 0.64, iy);
    a3->SetPointError(0, 0.32, 0.28, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Anselmino et al. (2013)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a3->DrawClone("pesame");
    t0->DrawClone("same");
    //Fit PLB767(2017)91 Ye et al.
    iy -= 1.0;
    a3->SetPoint(0, 0.64, iy);
    a3->SetPointError(0, 0.15, 0.15, 0, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{Ye et al. (2017)}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a3->DrawClone("pesame");
    t0->DrawClone("same");

    TGraphErrors * a4 = new TGraphErrors(1);
    a4->SetMarkerStyle(8);
    a4->SetMarkerColor(2);
    a4->SetMarkerSize(msize);
    a4->SetLineColor(2);
    a4->SetLineWidth(lwidth);
    //SoLID PLB767(2017)91 Ye et al. 
    iy -= 1.0;
    a4->SetPoint(0, 0.64, iy);
    a4->SetPointError(0, 0.021, 0);
    t0 = new TLatex(Xtext, iy, "#font[22]{JLab12 SoLID}");
    t0->SetTextAlign(12);
    t0->SetTextSize(tsize);
    a4->DrawClone("pesame");
    t0->DrawClone("same");

    TLegend * leg = new TLegend(0.60, 0.75, 0.9, 0.9);
    leg->AddEntry(a0, "#font[22]{Dyson-Schwinger equation}", "p");
    leg->AddEntry(a1, "#font[22]{Lattice QCD}", "p");
    leg->AddEntry(a3, "#font[22]{Phenomenology}", "p");
    leg->AddEntry(a4, "#font[22]{Future experiment}", "p");


    leg->Draw("same");
    c0->Print("gt1.pdf");
  }



  return 0;
}
