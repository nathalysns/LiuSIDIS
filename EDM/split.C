#include "splitsusy.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"

int main(const int argc, const char * argv[]){


  double gu = 1.025 * sqrt(2.0) * MW / VEV;
  double gd = 1.025 * sqrt(2.0) * MW / VEV;
  //double gt_now[2] = {10.33, -0.33};
  double gt_now[2] = {0.413, -0.229};
  //double et_now[3] = {0.0, 0.0, 0.0};
  double et_now[3] = {0.133, 0.094, 0.002};
  double et_future[3] = {0.018, 0.008, -2.65e-5};

  double par[7] = {gu, gd, gt_now[0], gt_now[1], et_now[0], et_now[1], et_now[2]};
  double par1[7] = {gu, gd, gt_now[0], gt_now[1], et_future[0], et_future[1], et_future[2]};

  double mu[200], mu1[200], M2[200];
  for (int i = 0; i < 200; i++){
    M2[i] = 2.0 * pow(10.0, 2.0 + 0.01 * i);
    mu[i] = Solve_mu_nEDM(M2[i], par);
    mu1[i] = Solve_mu_nEDM(M2[i], par1);
  }
  
  TH1D * h0 = new TH1D("h0", "", 1, 200.0, 10000.0);
  //h0->GetXaxis()->SetRangeUser(200.0, 600.0);
  h0->SetMinimum(2.0e2);
  h0->SetMaximum(1.0e4);
  h0->SetStats(0);
  h0->GetXaxis()->SetTitle("M_{2} (GeV)");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetTitleSize(0.06);
  h0->GetYaxis()->SetTitle("#mu (GeV)");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetTitleSize(0.06);
  

  TGraph * g0 = new TGraph(200, M2, mu);
  g0->SetLineColor(4);
  g0->SetLineStyle(1);
  TGraph * g1 = new TGraph(200, M2, mu1);
  g1->SetLineColor(2);
  g1->SetLineStyle(2);
  TCanvas * c0 = new TCanvas("c0", "", 800, 800);
  gPad->SetLogx();
  gPad->SetLogy();
  c0->SetBottomMargin(0.15);
  c0->SetLeftMargin(0.15);
  h0->DrawClone("axis");
  g0->Draw("lsame");
  g1->Draw("lsame");

  c0->Print("c1.pdf");
  
  return 0;
}
