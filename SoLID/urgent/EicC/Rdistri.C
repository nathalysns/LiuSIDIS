#include "sidis-region.h"
#include "eiccacceptance.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"


using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./Rdistri <opt>" << endl;
    cout << "./Rdistri <1> <hadron>" << endl;
    return 0;
  }
  
  const int opt = atoi(argv[1]);

  gRandom->SetSeed(0);

  const double Mp = 0.938272;

  TLorentzVector l(0, 0, 3.5, 3.5);
  TLorentzVector P1(0, 0, -sqrt(pow(20.0, 2) - Mp*Mp), 20.0);
  TLorentzVector P2(0, 0, -sqrt(pow(20.0*2/3, 2) - Mp*Mp), 20.0*2/3);
  TLorentzVector lp, Ph;
  double weight, acc, R1;
  double Xmin[6] = {0.001, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.999, 50.0, 0.7, 2.0, M_PI, M_PI};
  Long64_t Nsim = 10000000;
  double lumi = 1e+7 * pow(0.197327, 2);
  double time = 100.0 * 24.0 * 3600.0;

  Lsidis sidis1;
  Lsidis sidis2;

  sidis1.SetNucleus(1, 0);//proton
  sidis1.SetInitialState(l, P1);
  sidis1.SetPDFset("CJ15lo");
  sidis1.ChangeTMDpars(0.57, 0.12);
  sidis1.SetRange(Xmin, Xmax);

  sidis2.SetNucleus(0.67, 0.33);//helium-3
  sidis2.SetInitialState(l, P2);
  sidis2.SetPDFset("CJ15lo");
  sidis2.ChangeTMDpars(0.57, 0.12);
  sidis2.SetRange(Xmin, Xmax);
  
  if (opt == 0){//test
    return 0;
  }

  if (opt == 1){//R1 analysis
    sidis1.SetHadron(argv[2]);
    sidis1.SetFFset("DSSFFlo");
    sidis2.SetHadron(argv[2]);
    sidis2.SetFFset("DSSFFlo");

    TFile * fs = new TFile("R1.root", "RECREATE");
    TH1D * h1 = new TH1D("h1", "", 200, 0.0, 10.0);
    h1->SetDirectory(fs);
    TH1D * h2 = new TH1D("h2", "", 200, 0.0, 10.0);
    h2->SetDirectory(fs);
    TH2D * zPt1 = new TH2D("zPt1", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1->SetDirectory(fs);
    TH2D * zPt1cut = new TH2D("zPt1cut", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1cut->SetDirectory(fs);
    TH2D * zPt2 = new TH2D("zPt2", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1->SetDirectory(fs);
    TH2D * zPt2cut = new TH2D("zPt2cut", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1cut->SetDirectory(fs);

    int proc = 0;
    for (Long64_t i = 0; i < Nsim; i++){//proton
      if (i * 100 / Nsim > proc){
	proc = i * 100 / Nsim;
	cout << proc << "%" << endl;
      }
      weight = sidis1.GenerateEvent(0, 1);
      if (weight > 0){
	if (sidis1.GetVariable("W") < 2.3) continue;
	if (sidis1.GetVariable("Wp") < 1.6) continue;
	if (sidis1.GetVariable("Q2") < 1.0) continue;
	if (sidis1.GetVariable("z") < 0.3 || sidis1.GetVariable("z") > 0.7) continue;
	lp = sidis1.GetLorentzVector("lp");
	Ph = sidis1.GetLorentzVector("Ph");
	acc = EICC::GetAcceptance_e(lp) * EICC::GetAcceptance_hadron(Ph, argv[2]);
	if (acc > 0){
	  R1 = Calculate_R1(&sidis1);
	  h1->Fill(R1, weight * acc);
	  zPt1->Fill(sidis1.GetVariable("z"), sidis1.GetVariable("Pt"), weight * acc);
	  if (R1 < 0.4){
	    zPt1cut->Fill(sidis1.GetVariable("z"), sidis1.GetVariable("Pt"), weight * acc);
	  }
	}
      }

      weight = sidis2.GenerateEvent(0, 1);
      if (weight > 0){
	if (sidis2.GetVariable("W") < 2.3) continue;
	if (sidis2.GetVariable("Wp") < 1.6) continue;
	if (sidis2.GetVariable("Q2") < 1.0) continue;
	if (sidis2.GetVariable("z") < 0.3 || sidis2.GetVariable("z") > 0.7) continue;
	lp = sidis2.GetLorentzVector("lp");
	Ph = sidis2.GetLorentzVector("Ph");
	acc = EICC::GetAcceptance_e(lp) * EICC::GetAcceptance_hadron(Ph, argv[2]);
	if (acc > 0){
	  R1 = Calculate_R1(&sidis2);
	  h2->Fill(R1, weight * acc);
	  zPt2->Fill(sidis2.GetVariable("z"), sidis2.GetVariable("Pt"), weight * acc);
	  if (R1 < 0.4){
	    zPt2cut->Fill(sidis2.GetVariable("z"), sidis2.GetVariable("Pt"), weight * acc);
	  }
	}
      }
    }
    h1->Scale(lumi*time/Nsim);
    h2->Scale(lumi*time/Nsim);
    zPt1->Scale(lumi*time/Nsim);
    zPt1cut->Scale(lumi*time/Nsim);
    zPt2->Scale(lumi*time/Nsim);
    zPt2cut->Scale(lumi*time/Nsim);

    fs->Write();
    return 0;
  }


  if (opt == 2){//plot pi+
    TFile * fs = new TFile("R1-pip.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("h1"); h1->Scale(1 / time);
    TH1D * h2 = (TH1D *) fs->Get("h2"); h2->Scale(3 / time);
    TH2D * zPt1 = (TH2D *) fs->Get("zPt1"); zPt1->Scale(1 / time);
    TH2D * zPt2 = (TH2D *) fs->Get("zPt2"); zPt2->Scale(3 / time);
    TH2D * zPt1cut = (TH2D *) fs->Get("zPt1cut"); zPt1cut->Scale(1 / time);
    TH2D * zPt2cut = (TH2D *) fs->Get("zPt2cut"); zPt2cut->Scale(3 / time);

    TH1D * h0 = new TH1D("h0", "", 1, 0.0, 3.0);
    h0->SetTitle("");
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("#font[22]{R_{1}}");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.06);
    h0->GetXaxis()->SetRangeUser(0.0, 3.0);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("#font[22]{rate}");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetRangeUser(0.0, 5.0);
    h0->GetYaxis()->SetNdivisions(6, 5, 0);
    h0->SetMaximum(h2->GetMaximum()*1.2);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h0->DrawClone("axis");

    h1->SetStats(0);
    h1->SetLineColor(4);
    h1->SetLineWidth(2);
    h1->DrawClone("same");

    h2->SetStats(0);
    h2->SetLineColor(2);
    h2->SetLineWidth(2);
    h2->DrawClone("same");

    TLegend * l0 = new TLegend(0.75, 0.72, 0.9, 0.9);
    l0->AddEntry(h1, "#font[22]{ep}", "l");
    l0->AddEntry(h2, "#font[22]{e^{3}He}", "l");
    l0->DrawClone("same");

    TLine * V1 = new TLine(0.4, 0, 0.4, h0->GetMaximum());
    V1->SetLineColor(1);
    V1->SetLineWidth(1);
    V1->SetLineStyle(2);
    V1->DrawClone("same");

    TLine * V2 = new TLine(2.5, 0, 2.5, h0->GetMaximum());
    V2->SetLineColor(1);
    V2->SetLineWidth(1);
    V2->SetLineStyle(2);
    V2->DrawClone("same");

    c0->Print("plots/R1-pip.pdf");

    TH2D * zPt0 = new TH2D("zPt0", "", 1, 0.2, 0.8, 1, 0.0, 2.0);
    zPt0->SetTitle("");
    zPt0->SetStats(0);
    zPt0->GetXaxis()->SetTitle("#font[22]{z}");
    zPt0->GetXaxis()->CenterTitle(true);
    zPt0->GetXaxis()->SetTitleSize(0.06);
    zPt0->GetXaxis()->SetTitleOffset(1.15);
    zPt0->GetXaxis()->SetLabelSize(0.06);
    zPt0->GetXaxis()->SetRangeUser(0.2, 0.8);
    zPt0->GetXaxis()->SetNdivisions(6, 5, 0);
    zPt0->GetYaxis()->SetTitle("#font[22]{P_{hT} (GeV)}");
    zPt0->GetYaxis()->CenterTitle(true);
    zPt0->GetYaxis()->SetTitleSize(0.06);
    zPt0->GetYaxis()->SetTitleOffset(1.15);
    zPt0->GetYaxis()->SetLabelSize(0.06);
    zPt0->GetYaxis()->SetRangeUser(0.0, 2.0);
    zPt0->GetYaxis()->SetNdivisions(6, 5, 0);
    zPt0->SetMaximum(zPt1->GetMaximum());
    
    c0->Clear();
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    c0->SetLogz();

    zPt0->DrawClone("axis");
    zPt1->DrawClone("colzsame");
    c0->Print("plots/zPt-pip.pdf(");

    zPt0->DrawClone("axis");
    zPt1cut->DrawClone("colzsame");
    c0->Print("plots/zPt-pip.pdf");

    zPt0->SetMaximum(zPt2->GetMaximum());
    zPt0->DrawClone("axis");
    zPt2->DrawClone("colzsame");
    c0->Print("plots/zPt-pip.pdf");

    zPt0->DrawClone("axis");
    zPt2cut->DrawClone("colzsame");
    c0->Print("plots/zPt-pip.pdf)");
  

    return 0;
  }
      
  if (opt == 3){//plot pi-
    TFile * fs = new TFile("R1-pim.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("h1"); h1->Scale(1 / time);
    TH1D * h2 = (TH1D *) fs->Get("h2"); h2->Scale(3 / time);
    TH2D * zPt1 = (TH2D *) fs->Get("zPt1"); zPt1->Scale(1 / time);
    TH2D * zPt2 = (TH2D *) fs->Get("zPt2"); zPt2->Scale(3 / time);
    TH2D * zPt1cut = (TH2D *) fs->Get("zPt1cut"); zPt1cut->Scale(1 / time);
    TH2D * zPt2cut = (TH2D *) fs->Get("zPt2cut"); zPt2cut->Scale(3 / time);

    TH1D * h0 = new TH1D("h0", "", 1, 0.0, 3.0);
    h0->SetTitle("");
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("#font[22]{R_{1}}");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.06);
    h0->GetXaxis()->SetRangeUser(0.0, 3.0);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("#font[22]{rate}");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetRangeUser(0.0, 5.0);
    h0->GetYaxis()->SetNdivisions(6, 5, 0);
    h0->SetMaximum(h2->GetMaximum()*1.2);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h0->DrawClone("axis");

    h1->SetStats(0);
    h1->SetLineColor(4);
    h1->SetLineWidth(2);
    h1->DrawClone("same");

    h2->SetStats(0);
    h2->SetLineColor(2);
    h2->SetLineWidth(2);
    h2->DrawClone("same");

    TLegend * l0 = new TLegend(0.75, 0.72, 0.9, 0.9);
    l0->AddEntry(h1, "#font[22]{ep}", "l");
    l0->AddEntry(h2, "#font[22]{e^{3}He}", "l");
    l0->DrawClone("same");

    TLine * V1 = new TLine(0.4, 0, 0.4, h0->GetMaximum());
    V1->SetLineColor(1);
    V1->SetLineWidth(1);
    V1->SetLineStyle(2);
    V1->DrawClone("same");

    TLine * V2 = new TLine(2.5, 0, 2.5, h0->GetMaximum());
    V2->SetLineColor(1);
    V2->SetLineWidth(1);
    V2->SetLineStyle(2);
    V2->DrawClone("same");

    c0->Print("plots/R1-pim.pdf");

    TH2D * zPt0 = new TH2D("zPt0", "", 1, 0.2, 0.8, 1, 0.0, 2.0);
    zPt0->SetTitle("");
    zPt0->SetStats(0);
    zPt0->GetXaxis()->SetTitle("#font[22]{z}");
    zPt0->GetXaxis()->CenterTitle(true);
    zPt0->GetXaxis()->SetTitleSize(0.06);
    zPt0->GetXaxis()->SetTitleOffset(1.15);
    zPt0->GetXaxis()->SetLabelSize(0.06);
    zPt0->GetXaxis()->SetRangeUser(0.2, 0.8);
    zPt0->GetXaxis()->SetNdivisions(6, 5, 0);
    zPt0->GetYaxis()->SetTitle("#font[22]{P_{hT} (GeV)}");
    zPt0->GetYaxis()->CenterTitle(true);
    zPt0->GetYaxis()->SetTitleSize(0.06);
    zPt0->GetYaxis()->SetTitleOffset(1.15);
    zPt0->GetYaxis()->SetLabelSize(0.06);
    zPt0->GetYaxis()->SetRangeUser(0.0, 2.0);
    zPt0->GetYaxis()->SetNdivisions(6, 5, 0);
    zPt0->SetMaximum(zPt1->GetMaximum());
    
    c0->Clear();
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    c0->SetLogz();

    zPt0->DrawClone("axis");
    zPt1->DrawClone("colzsame");
    c0->Print("plots/zPt-pim.pdf(");

    zPt0->DrawClone("axis");
    zPt1cut->DrawClone("colzsame");
    c0->Print("plots/zPt-pim.pdf");

    zPt0->SetMaximum(zPt2->GetMaximum());
    zPt0->DrawClone("axis");
    zPt2->DrawClone("colzsame");
    c0->Print("plots/zPt-pim.pdf");

    zPt0->DrawClone("axis");
    zPt2cut->DrawClone("colzsame");
    c0->Print("plots/zPt-pim.pdf)");
  

    return 0;
  }

  if (opt == 4){//plot K+
    TFile * fs = new TFile("R1-kp.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("h1"); h1->Scale(1 / time);
    TH1D * h2 = (TH1D *) fs->Get("h2"); h2->Scale(3 / time);
    TH2D * zPt1 = (TH2D *) fs->Get("zPt1"); zPt1->Scale(1 / time);
    TH2D * zPt2 = (TH2D *) fs->Get("zPt2"); zPt2->Scale(3 / time);
    TH2D * zPt1cut = (TH2D *) fs->Get("zPt1cut"); zPt1cut->Scale(1 / time);
    TH2D * zPt2cut = (TH2D *) fs->Get("zPt2cut"); zPt2cut->Scale(3 / time);

    TH1D * h0 = new TH1D("h0", "", 1, 0.0, 3.0);
    h0->SetTitle("");
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("#font[22]{R_{1}}");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.06);
    h0->GetXaxis()->SetRangeUser(0.0, 3.0);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("#font[22]{rate}");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetRangeUser(0.0, 5.0);
    h0->GetYaxis()->SetNdivisions(6, 5, 0);
    h0->SetMaximum(h2->GetMaximum()*1.2);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h0->DrawClone("axis");

    h1->SetStats(0);
    h1->SetLineColor(4);
    h1->SetLineWidth(2);
    h1->DrawClone("same");

    h2->SetStats(0);
    h2->SetLineColor(2);
    h2->SetLineWidth(2);
    h2->DrawClone("same");

    TLegend * l0 = new TLegend(0.75, 0.72, 0.9, 0.9);
    l0->AddEntry(h1, "#font[22]{ep}", "l");
    l0->AddEntry(h2, "#font[22]{e^{3}He}", "l");
    l0->DrawClone("same");

    TLine * V1 = new TLine(0.4, 0, 0.4, h0->GetMaximum());
    V1->SetLineColor(1);
    V1->SetLineWidth(1);
    V1->SetLineStyle(2);
    V1->DrawClone("same");

    TLine * V2 = new TLine(2.5, 0, 2.5, h0->GetMaximum());
    V2->SetLineColor(1);
    V2->SetLineWidth(1);
    V2->SetLineStyle(2);
    V2->DrawClone("same");

    c0->Print("plots/R1-kp.pdf");

    TH2D * zPt0 = new TH2D("zPt0", "", 1, 0.2, 0.8, 1, 0.0, 2.0);
    zPt0->SetTitle("");
    zPt0->SetStats(0);
    zPt0->GetXaxis()->SetTitle("#font[22]{z}");
    zPt0->GetXaxis()->CenterTitle(true);
    zPt0->GetXaxis()->SetTitleSize(0.06);
    zPt0->GetXaxis()->SetTitleOffset(1.15);
    zPt0->GetXaxis()->SetLabelSize(0.06);
    zPt0->GetXaxis()->SetRangeUser(0.2, 0.8);
    zPt0->GetXaxis()->SetNdivisions(6, 5, 0);
    zPt0->GetYaxis()->SetTitle("#font[22]{P_{hT} (GeV)}");
    zPt0->GetYaxis()->CenterTitle(true);
    zPt0->GetYaxis()->SetTitleSize(0.06);
    zPt0->GetYaxis()->SetTitleOffset(1.15);
    zPt0->GetYaxis()->SetLabelSize(0.06);
    zPt0->GetYaxis()->SetRangeUser(0.0, 2.0);
    zPt0->GetYaxis()->SetNdivisions(6, 5, 0);
    zPt0->SetMaximum(zPt1->GetMaximum());
    
    c0->Clear();
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    c0->SetLogz();

    zPt0->DrawClone("axis");
    zPt1->DrawClone("colzsame");
    c0->Print("plots/zPt-kp.pdf(");

    zPt0->DrawClone("axis");
    zPt1cut->DrawClone("colzsame");
    c0->Print("plots/zPt-kp.pdf");

    zPt0->SetMaximum(zPt2->GetMaximum());
    zPt0->DrawClone("axis");
    zPt2->DrawClone("colzsame");
    c0->Print("plots/zPt-kp.pdf");

    zPt0->DrawClone("axis");
    zPt2cut->DrawClone("colzsame");
    c0->Print("plots/zPt-kp.pdf)");
  

    return 0;
  }

  if (opt == 5){//plot K-
    TFile * fs = new TFile("R1-km.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("h1"); h1->Scale(1 / time);
    TH1D * h2 = (TH1D *) fs->Get("h2"); h2->Scale(3 / time);
    TH2D * zPt1 = (TH2D *) fs->Get("zPt1"); zPt1->Scale(1 / time);
    TH2D * zPt2 = (TH2D *) fs->Get("zPt2"); zPt2->Scale(3 / time);
    TH2D * zPt1cut = (TH2D *) fs->Get("zPt1cut"); zPt1cut->Scale(1 / time);
    TH2D * zPt2cut = (TH2D *) fs->Get("zPt2cut"); zPt2cut->Scale(3 / time);

    TH1D * h0 = new TH1D("h0", "", 1, 0.0, 3.0);
    h0->SetTitle("");
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("#font[22]{R_{1}}");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.06);
    h0->GetXaxis()->SetRangeUser(0.0, 3.0);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("#font[22]{rate}");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetRangeUser(0.0, 5.0);
    h0->GetYaxis()->SetNdivisions(6, 5, 0);
    h0->SetMaximum(h2->GetMaximum()*1.2);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h0->DrawClone("axis");

    h1->SetStats(0);
    h1->SetLineColor(4);
    h1->SetLineWidth(2);
    h1->DrawClone("same");

    h2->SetStats(0);
    h2->SetLineColor(2);
    h2->SetLineWidth(2);
    h2->DrawClone("same");

    TLegend * l0 = new TLegend(0.75, 0.72, 0.9, 0.9);
    l0->AddEntry(h1, "#font[22]{ep}", "l");
    l0->AddEntry(h2, "#font[22]{e^{3}He}", "l");
    l0->DrawClone("same");

    TLine * V1 = new TLine(0.4, 0, 0.4, h0->GetMaximum());
    V1->SetLineColor(1);
    V1->SetLineWidth(1);
    V1->SetLineStyle(2);
    V1->DrawClone("same");

    TLine * V2 = new TLine(2.5, 0, 2.5, h0->GetMaximum());
    V2->SetLineColor(1);
    V2->SetLineWidth(1);
    V2->SetLineStyle(2);
    V2->DrawClone("same");

    c0->Print("plots/R1-km.pdf");

    TH2D * zPt0 = new TH2D("zPt0", "", 1, 0.2, 0.8, 1, 0.0, 2.0);
    zPt0->SetTitle("");
    zPt0->SetStats(0);
    zPt0->GetXaxis()->SetTitle("#font[22]{z}");
    zPt0->GetXaxis()->CenterTitle(true);
    zPt0->GetXaxis()->SetTitleSize(0.06);
    zPt0->GetXaxis()->SetTitleOffset(1.15);
    zPt0->GetXaxis()->SetLabelSize(0.06);
    zPt0->GetXaxis()->SetRangeUser(0.2, 0.8);
    zPt0->GetXaxis()->SetNdivisions(6, 5, 0);
    zPt0->GetYaxis()->SetTitle("#font[22]{P_{hT} (GeV)}");
    zPt0->GetYaxis()->CenterTitle(true);
    zPt0->GetYaxis()->SetTitleSize(0.06);
    zPt0->GetYaxis()->SetTitleOffset(1.15);
    zPt0->GetYaxis()->SetLabelSize(0.06);
    zPt0->GetYaxis()->SetRangeUser(0.0, 2.0);
    zPt0->GetYaxis()->SetNdivisions(6, 5, 0);
    zPt0->SetMaximum(zPt1->GetMaximum());
    
    c0->Clear();
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    c0->SetLogz();

    zPt0->DrawClone("axis");
    zPt1->DrawClone("colzsame");
    c0->Print("plots/zPt-km.pdf(");

    zPt0->DrawClone("axis");
    zPt1cut->DrawClone("colzsame");
    c0->Print("plots/zPt-km.pdf");

    zPt0->SetMaximum(zPt2->GetMaximum());
    zPt0->DrawClone("axis");
    zPt2->DrawClone("colzsame");
    c0->Print("plots/zPt-km.pdf");

    zPt0->DrawClone("axis");
    zPt2cut->DrawClone("colzsame");
    c0->Print("plots/zPt-km.pdf)");
  

    return 0;
  }


    
  

  return 0;
}
