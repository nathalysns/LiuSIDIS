#include "splitsusy.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"

int main(const int argc, const char * argv[]){\

  if (argc < 2) return 0;

  const int opt = atoi(argv[1]);

  double gu = 1.025 * sqrt(2.0) * MW / VEV;
  double gd = 1.025 * sqrt(2.0) * MW / VEV;
  double gt[2] = {0.413, -0.229};
  double et_now[3] = {0.133, 0.094, 0.002};
  double et_future[3] = {0.018, 0.008, -2.65e-5};

  double par0[7] = {gu, gd, gt[0], gt[1], 0.0, 0.0, 0.0};
  double par_now[7] = {gu, gd, gt[0], gt[1], et_now[0], et_now[1], et_now[2]};
  double par_future[7] = {gu, gd, gt[0], gt[1], et_future[0], et_future[1], et_future[2]};

  if (opt == 0){//nEDM constraint
    double mu0[200], mu_now[200], mu_future[200], M2[200];
    cout << "nEDM: " << nEDMlimit <<  "   running ..." << endl;
    for (int i = 0; i < 200; i++){
      M2[i] = 1.0 * pow(10.0, 2.0 + 0.01 * i);
      mu0[i] = Solve_mu_nEDM(M2[i], par0);
      mu_now[i] = Solve_mu_nEDM(M2[i], par_now);
      mu_future[i] = Solve_mu_nEDM(M2[i], par_future);
    }
    TGraph * g0 = new TGraph(200, M2, mu0);
    g0->SetLineColor(1);
    g0->SetLineWidth(1);
    g0->SetLineStyle(2);
    TGraph * g_now = new TGraph(200, M2, mu_now);
    g_now->SetLineColor(4);
    g_now->SetLineWidth(2);
    g_now->SetLineStyle(1);
    TGraph * g_future = new TGraph(200, M2, mu_future);
    g_future->SetLineColor(2);
    g_future->SetLineWidth(2);
    g_future->SetLineStyle(7);

    TLegend * l0 = new TLegend(0.5, 0.75, 0.89, 0.89);
    l0->SetBorderSize(0);
    l0->SetFillStyle(0);
    l0->AddEntry(g_now, "#font[22]{future d_{n} + current g_{T}}", "l");
    l0->AddEntry(g_future, "#font[22]{future d_{n} + future g_{T}}", "l");

    TH1D * h0 = new TH1D("h0", "", 1, 200.0, 1.0e4);
    //h0->GetXaxis()->SetRangeUser(200.0, 600.0);
    h0->SetMinimum(2.0e2);
    h0->SetMaximum(1.0e4);
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("#font[22]{M_{2}} (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetLabelSize(0.055);
    h0->GetXaxis()->SetLabelFont(22);
    h0->GetYaxis()->SetTitle("#font[22]{#mu} (GeV)");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetLabelSize(0.055);
    h0->GetYaxis()->SetLabelFont(22);
  
    TCanvas * c0 = new TCanvas("c0", "", 800, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h0->DrawClone("axis");
    g0->DrawClone("lsame");
    g_now->DrawClone("lsame");
    g_future->DrawClone("lsame");
    

    nEDMlimit = 1.0e-27;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_nEDM(M2[i], par0);
    }
    TGraph * g27 = new TGraph(200, M2, mu0);
    g27->SetLineColor(11);
    g27->SetLineStyle(2);
    g27->DrawClone("lsame");

    nEDMlimit = 1.0e-28;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_nEDM(M2[i], par0);
    }
    TGraph * g28 = new TGraph(200, M2, mu0);
    g28->SetLineColor(11);
    g28->SetLineStyle(2);
    g28->DrawClone("lsame");

    nEDMlimit = 1.0e-29;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_nEDM(M2[i], par0);
    }
    TGraph * g29 = new TGraph(200, M2, mu0);
    g29->SetLineColor(11);
    g29->SetLineStyle(2);
    g29->DrawClone("lsame");

    TLatex * lt = new TLatex();
    lt->SetTextSize(0.03);
    lt->SetTextFont(132);
    lt->DrawLatex(230.0, 350.0, "10^{-27}");
    lt->DrawLatex(1.3e3, 1.3e3, "10^{-28}");
    lt->DrawLatex(6.0e3, 3.0e3, "10^{-29}");

    lt->SetTextFont(22);
    lt->DrawLatex(5.0e2, 7.0e2, "3.0#times10^{-28}");

    l0->DrawClone("same");
    c0->Print("nedmsplit.pdf");
  }

  if (opt == 1){//pEDM constraint
    double mu0[200], mu_now[200], mu_future[200], M2[200];
    cout << "pEDM: " << pEDMlimit <<  "   running ..." << endl;
    for (int i = 0; i < 200; i++){
      M2[i] = 1.0 * pow(10.0, 2.0 + log(5.0e4/200.0) / 199.0 * i);
      mu0[i] = Solve_mu_pEDM(M2[i], par0);
      mu_now[i] = Solve_mu_pEDM(M2[i], par_now);
      mu_future[i] = Solve_mu_pEDM(M2[i], par_future);
    }
    TGraph * g0 = new TGraph(200, M2, mu0);
    g0->SetLineColor(1);
    g0->SetLineWidth(1);
    g0->SetLineStyle(2);
    TGraph * g_now = new TGraph(200, M2, mu_now);
    g_now->SetLineColor(4);
    g_now->SetLineWidth(2);
    g_now->SetLineStyle(1);
    TGraph * g_future = new TGraph(200, M2, mu_future);
    g_future->SetLineColor(2);
    g_future->SetLineWidth(2);
    g_future->SetLineStyle(7);

    TLegend * l0 = new TLegend(0.5, 0.75, 0.89, 0.89);
    l0->SetBorderSize(0);
    l0->SetFillStyle(0);
    l0->AddEntry(g_now, "#font[22]{future d_{p} + current g_{T}}", "l");
    l0->AddEntry(g_future, "#font[22]{future d_{p} + future g_{T}}", "l");

    TH1D * h0 = new TH1D("h0", "", 1, 200.0, 5.0e4);
    //h0->GetXaxis()->SetRangeUser(200.0, 600.0);
    h0->SetMinimum(2.0e2);
    h0->SetMaximum(5.0e4);
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("#font[22]{M_{2}} (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetLabelSize(0.055);
    h0->GetXaxis()->SetLabelFont(22);
    h0->GetYaxis()->SetTitle("#font[22]{#mu} (GeV)");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetLabelSize(0.055);
    h0->GetYaxis()->SetLabelFont(22);
  
    TCanvas * c0 = new TCanvas("c0", "", 800, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h0->DrawClone("axis");
    g0->DrawClone("lsame");
    g_now->DrawClone("lsame");
    g_future->DrawClone("lsame");
    

    pEDMlimit = 1.0e-27;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_pEDM(M2[i], par0);
    }
    TGraph * g27 = new TGraph(200, M2, mu0);
    g27->SetLineColor(11);
    g27->SetLineStyle(2);
    g27->DrawClone("lsame");

    pEDMlimit = 1.0e-28;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_pEDM(M2[i], par0);
    }
    TGraph * g28 = new TGraph(200, M2, mu0);
    g28->SetLineColor(11);
    g28->SetLineStyle(2);
    g28->DrawClone("lsame");

    pEDMlimit = 1.0e-29;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_pEDM(M2[i], par0);
    }
    TGraph * g29 = new TGraph(200, M2, mu0);
    g29->SetLineColor(11);
    g29->SetLineStyle(2);
    g29->DrawClone("lsame");

    pEDMlimit = 1.0e-30;
    for (int i = 0; i < 200; i++){
      mu0[i] = Solve_mu_pEDM(M2[i], par0);
    }
    TGraph * g30 = new TGraph(200, M2, mu0);
    g30->SetLineColor(11);
    g30->SetLineStyle(2);
    g30->DrawClone("lsame");

    TLatex * lt = new TLatex();
    lt->SetTextSize(0.03);
    lt->SetTextFont(132);
    lt->DrawLatex(250.0, 350.0, "10^{-27}");
    lt->DrawLatex(1.3e3, 1.3e3, "10^{-28}");
    lt->DrawLatex(6.0e3, 3.5e3, "10^{-29}");
    lt->DrawLatex(3.0e4, 0.9e4, "10^{-30}");

    lt->SetTextFont(22);
    lt->DrawLatex(2.5e3, 2.3e3, "2.6#times10^{-29}");

    l0->DrawClone("same");
    c0->Print("pedmsplit.pdf");
  }


  return 0;
}
