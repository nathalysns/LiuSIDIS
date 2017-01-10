#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TChain.h"

#include "Lsidis.h"

// Acceptance 
TFile * file_e = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_electron_output.root", "r");
TFile * file_pip = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_pionp_output.root", "r");
TFile * file_pim = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_pionm_output.root", "r");
TH3F * acc_FA_e = (TH3F *) file_e->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_e = (TH3F *) file_e->Get("acceptance_ThetaPhiP_largeangle");
TH3F * acc_FA_pip = (TH3F *) file_pip->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_pip = (TH3F *) file_pip->Get("acceptance_ThetaPhiP_largeangle");
TH3F * acc_FA_pim = (TH3F *) file_pim->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_pim = (TH3F *) file_pim->Get("acceptance_ThetaPhiP_largeangle");

double Rfactor0 = 1.0e5;

double GetAcceptance_e(const TLorentzVector p){//Get electron acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 50.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(phi), acc_FA_e->GetZaxis()->FindBin(mom));
  if (mom > 3.5)
    acc += acc_LA_e->GetBinContent(acc_LA_e->GetXaxis()->FindBin(theta), acc_LA_e->GetYaxis()->FindBin(phi), acc_LA_e->GetZaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_pi(const TLorentzVector p, const char * hadron){//Get pion acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 45.0) return 0;
  double mom = p.P();
  double acc = 0;
  if (strcmp(hadron, "pi+") == 0)
    acc += acc_FA_pip->GetBinContent(acc_FA_pip->GetXaxis()->FindBin(theta), acc_FA_pip->GetYaxis()->FindBin(phi), acc_FA_pip->GetZaxis()->FindBin(mom));
  else if (strcmp(hadron, "pi-") == 0)
    acc += acc_FA_pim->GetBinContent(acc_FA_pim->GetXaxis()->FindBin(theta), acc_FA_pim->GetYaxis()->FindBin(phi), acc_FA_pim->GetZaxis()->FindBin(mom));
  else return 0;
  return acc;
}

int GetTotalRate(const double Ebeam, const char * hadron){//Estimate the total rate
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 8.0, 0.7, 1.6, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  double sum = 0.0;
  double Nsim = 1.0e7;
  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      sidis.CalculateRfactor();
      if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      sum += weight * GetAcceptance_e(lp) * GetAcceptance_pi(Ph, hadron);
    }
  }
  std::cout << "Total rate: " << sum * lumi / Nsim << std::endl;
  return 0;
}

int GenerateBinInfoFile(const char * filename, const double Ebeam, const char * hadron){//Bin the data and create the bin info file
  FILE * fp = fopen(filename, "w");
  fprintf(fp, "Q2l\t Q2u\t zl\t zu\t Ptl\t Ptu\t xl\t xu\n");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 55.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 27.5 * 24.0 * 3600.0;
  double Nsim = 1.0e6;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double acc = 0;
  int Nx = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  double Q2list[7] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
  double statlist[6] = {1.0e7, 6.4e6, 3.2e6, 1.6e6, 1.2e6, 1.0e6};
  double zlist[9] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7};
  double Ptlist[7] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.6};
  int xi = 1;
  for (int Qi = 0; Qi < 6; Qi++){//Q2 loop
    Xmin[1] = Q2list[Qi];
    Xmax[1] = Q2list[Qi+1];
    for (int zi = 0; zi < 8; zi++){//z loop
      Xmin[2] = zlist[zi];
      Xmax[2] = zlist[zi+1];
      Xmin[3] = Ptlist[0];
      for (int kj = 1; kj < 7;){//Pt loop
	Xmax[3] = Ptlist[kj];
	sidis.SetRange(Xmin, Xmax);
	TH1D * hx = new TH1D("hx", "hx", 7000, 0.0, 0.7);
	printf("Q2:%.1f-%.1f  z:%.2f-%.2f  Pt:%.1f-%.1f\n",
	       Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3]);
	for (Long64_t i = 0; i < Nsim; i++){//generate events
	  weight = sidis.GenerateEvent(0, 1);
	  if (weight > 0){
            if (sidis.GetVariable("W") < 2.3) continue;
            if (sidis.GetVariable("Wp") < 1.6) continue;
	    sidis.CalculateRfactor();
	    if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	    lp = sidis.GetLorentzVector("lp");
	    Ph = sidis.GetLorentzVector("Ph");
	    acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph, hadron);
	    if (acc > 0)
	      hx->Fill(sidis.GetVariable("x"), weight * acc);
	  }
	}
	hx->Scale(lumi * time * eff / Nsim);
	if ((hx->Integral(1, -1) < statlist[Qi] && kj < 6) || (hx->Integral(1, -1) < 0.2 * statlist[Qi] && kj == 6)){
	  hx->Delete();
	  kj++;
	  continue;
	}
	Nx = 0;
	xi = 1;
	for (int xj = 1; xj <= 7000; xj++){
	  if (hx->Integral(xi, xj) > statlist[Qi] || xj == 7000){
	    fprintf(fp, "%.1f\t %.1f\t %.2f\t %.2f\t %.1f\t %.1f\t %.4f\t %.4f\n",
		    Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3],
		    hx->GetBinLowEdge(xi), hx->GetBinLowEdge(xj+1));
	    Nx++;
	    xi = xj + 1;
	  }
	}
	std::cout << Nx << std::endl;
	hx->Delete();
	Xmin[3] = Ptlist[kj];
	kj++;
      }
    }
  }
  fclose(fp);
  return 0;
}

int AnalyzeEstatUT3(const char * readfile, const char * savefile, const double Ebeam, const char * had){//bin analysis including stat. errors
  double Hadron = 0;
  if (strcmp(had, "pi+") == 0) Hadron = 0;
  else if (strcmp(had, "pi-") == 0) Hadron = 1;
  double Nucleon = 1;
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  double Eb = Ebeam;
  double x, y, z, Q2, Pt, phih, phiS;
  double dx, dy, dz, dQ2, dPt, dphih, dphiS, dv;
  double Nacc, fp;
  double Estatraw[3], Estat[3];
  Ts->Branch("Nucleon", &Nucleon, "Nucleon/D");
  Ts->Branch("Hadron", &Hadron, "Hadron/D");
  Ts->Branch("Ebeam", &Eb, "Ebeam/D");
  Ts->Branch("x", &x, "x/D");
  Ts->Branch("y", &y, "y/D");
  Ts->Branch("z", &z, "z/D");
  Ts->Branch("Q2", &Q2, "Q2/D");
  Ts->Branch("Pt", &Pt, "Pt/D");
  Ts->Branch("dx", &dx, "dx/D");
  Ts->Branch("dy", &dy, "dy/D");
  Ts->Branch("dz", &dz, "dz/D");
  Ts->Branch("dQ2", &dQ2, "dQ2/D");
  Ts->Branch("dPt", &dPt, "dPt/D");
  Ts->Branch("dphih", &dphih, "dphih/D");
  Ts->Branch("dphiS", &dphiS, "dphiS/D");
  Ts->Branch("dv", &dv, "dv/D");
  Ts->Branch("Nacc", &Nacc, "Nacc/D");
  Ts->Branch("fp", &fp, "fp/D");
  Ts->Branch("E0statraw", &Estatraw[0], "E0statraw/D");
  Ts->Branch("E1statraw", &Estatraw[1], "E1statraw/D");
  Ts->Branch("E2statraw", &Estatraw[2], "E2statraw/D");
  Ts->Branch("E0stat", &Estat[0], "E0stat/D");
  Ts->Branch("E1stat", &Estat[1], "E1stat/D");
  Ts->Branch("E2stat", &Estat[2], "E2stat/D");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0);
  sidis.SetHadron(had);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 55.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 27.5 * 24.0 * 3600.0;
  Long64_t Nsim = 0;
  Long64_t Nrec = 0;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double weight_p = 0;
  double acc = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  Lsidis sidis_p;
  sidis_p.SetNucleus(1, 0);
  sidis_p.SetHadron(had);
  sidis_p.SetInitialState(l, P);
  sidis_p.SetPDFset("CT14lo");
  sidis_p.SetFFset("DSSFFlo");
  ifstream infile(readfile);
  char tmp[300];
  infile.getline(tmp, 256);
  int Nt = 0;
  while (infile >> Xmin[1] >> Xmax[1] >> Xmin[2] >> Xmax[2] >> Xmin[3] >> Xmax[3] >> Xmin[0] >> Xmax[0]){
    printf("%.4d  Q2[%.1f,%.1f]  z[%.2f,%.2f]  Pt[%.1f,%.1f]  x[%.4f,%.4f]\n",
	   Nt++, Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3], Xmin[0], Xmax[0]);
    sidis.SetRange(Xmin, Xmax);
    sidis_p.SetRange(Xmin, Xmax);
    TH1D * hvar = new TH1D("hvar", "hvar", 7, -0.5, 6.5);
    TH2D * hs = new TH2D("hs", "hs", 36, -M_PI, M_PI, 18, 0, M_PI);
    Nsim = 0;
    Nrec = 0;
    for (Long64_t i = 0; i < 1.0e7; i++){
      Nsim++;
      weight = sidis.GenerateEvent(0, 1);
      if (weight > 0){
        if (sidis.GetVariable("W") < 2.3) continue;
        if (sidis.GetVariable("Wp") < 1.6) continue;
	sidis.CalculateRfactor();
	if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	lp = sidis.GetLorentzVector("lp");
	Ph = sidis.GetLorentzVector("Ph");
	acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph, had);
	if (acc > 0){
	  sidis_p.SetFinalState(lp, Ph);
	  sidis_p.CalculateVariables();
	  weight_p = sidis_p.GetEventWeight(0, 1);
	  Nrec++;
	  hvar->Fill(0., weight_p * acc);
	  hvar->Fill(1., weight * acc);
	  hvar->Fill(2., weight * acc * sidis.GetVariable("x"));
	  hvar->Fill(3., weight * acc * sidis.GetVariable("y"));
	  hvar->Fill(4., weight * acc * sidis.GetVariable("z"));
	  hvar->Fill(5., weight * acc * sidis.GetVariable("Q2"));
	  hvar->Fill(6., weight * acc * sidis.GetVariable("Pt"));
	  hs->Fill(sidis.GetVariable("phih"), std::abs(sidis.GetVariable("phiS")), weight * acc);
	}
      }
      if (Nrec > 500000) break;
    }
    hvar->Scale(lumi * time * eff / Nsim);
    hs->Scale(lumi * time * eff / Nsim);
    Nacc = hvar->GetBinContent(2);
    fp = hvar->GetBinContent(1) / Nacc;
    x = hvar->GetBinContent(3) / Nacc;
    y = hvar->GetBinContent(4) / Nacc;
    z = hvar->GetBinContent(5) / Nacc;
    Q2 = hvar->GetBinContent(6) / Nacc;
    Pt = hvar->GetBinContent(7) / Nacc;
    TMatrixD MUT3(3,3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	MUT3(i,j) = 0.0;
    for (int i = 1; i <= 36; i++){
      for (int j = 1; j <= 18; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	phiS = hs->GetYaxis()->GetBinCenter(j);
	MUT3(0,0) += sin(phih - phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(0,1) += sin(phih - phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(0,2) += sin(phih - phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,0) += sin(phih + phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,1) += sin(phih + phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,2) += sin(phih + phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,0) += sin(3.0 * phih - phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,1) += sin(3.0 * phih - phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,2) += sin(3.0 * phih - phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
      }
    }
    MUT3.Invert();
    for (int i = 0; i < 3; i++){
      Estatraw[i] = sqrt(2.0 * M_PI * M_PI / Nacc * (pow(MUT3(i,0),2) + pow(MUT3(i,1), 2) + pow(MUT3(i,2), 2)) * M_PI * M_PI);
      Estat[i] = Estatraw[i] / fp / 0.7;
      if (isnan(Estat[i]))
	std::cout << "NaN warning in Estat!" << std::endl;
    }
    Ts->Fill();
    hvar->Delete();
    hs->Delete();
  }
  fs->Write();
  infile.close();
  return 0;
}
  
double CheckCurrentCut(const double Ebeam, const char * hadron, const double kT2 = 0.5, const double MiT2 = 0.5, const double MfT2 = 0.5, const char * plotname = 0){
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double Nsim = 1.0e7;
  TH2D * h0 = new TH2D("h0", "", 1, 0.2, 0.8, 1, 0.0, 1.6);
  h0->GetXaxis()->SetTitle("z");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleSize(0.05);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetTitle("P_{hT} / GeV");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleSize(0.05);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelSize(0.055);
  TH2D * hall = new TH2D("hall", "Before cut", 60, 0.2, 0.8, 160, 0.0, 1.6);
  TH2D * hcut = new TH2D("hcut", "After cut", 60, 0.2, 0.8, 160, 0.0, 1.6);
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 8.0, 0.7, 1.6, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  double weight = 0;
  double acc = 0;
  TLorentzVector lp, Ph;
  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph, hadron);
      if (acc > 0){
	hall->Fill(sidis.GetVariable("z"), sidis.GetVariable("Pt"), weight * acc);
	sidis.CalculateRfactor(kT2, MiT2, MfT2);
	if (sidis.GetVariable("Rfactor") < 0.4){
	  hcut->Fill(sidis.GetVariable("z"), sidis.GetVariable("Pt"), weight * acc);
	}
      }
    }
  }
  hall->Scale(lumi/Nsim);
  hcut->Scale(lumi/Nsim);
  double rate = hcut->Integral(1, -1);
  std::cout << "All: " << hall->Integral(1, -1) << "   Cut: " << hcut->Integral(1, -1) << std::endl;
  if (plotname != 0){
    gStyle->SetOptStat(0);
    //hall->GetZaxis()->SetRangeUser(0.01, hall->GetMaximum()/0.95);
    hcut->GetZaxis()->SetRangeUser(0.01, hall->GetMaximum());
    TCanvas * c0 = new TCanvas("c0", "", 1600, 600);
    c0->SetBorderMode(0);
    c0->SetBorderSize(2);
    c0->SetFrameBorderMode(0);
    c0->Divide(2, 1);
    c0->cd(1);
    c0->cd(1)->SetLeftMargin(0.15);
    c0->cd(1)->SetBottomMargin(0.15);
    h0->Draw();
    hall->Draw("samecolz");
    c0->cd(2);
    c0->cd(2)->SetLeftMargin(0.15);
    c0->cd(2)->SetBottomMargin(0.15);
    h0->Draw();
    hcut->Draw("samecolz");
    c0->Print(plotname);
    c0->Close();
  }
  h0->Delete();
  hall->Delete();
  hcut->Delete();
  return rate;
}

int CreateFileSivers(const char * rootfile1, const char * rootfile2, const char * csvfile){//Create file for Sivers analysis use
  TChain * Ts = new TChain("data", "data");
  Ts->Add(rootfile1);
  Ts->Add(rootfile2);
  double Nucleon, Hadron, Ebeam, x, y, z, Q2, Pt, stat, systrel, systabs, fp;
  Ts->SetBranchAddress("Nucleon", &Nucleon);
  Ts->SetBranchAddress("Hadron", &Hadron);
  Ts->SetBranchAddress("Ebeam", &Ebeam);
  Ts->SetBranchAddress("x", &x);
  Ts->SetBranchAddress("y", &y);
  Ts->SetBranchAddress("z", &z);
  Ts->SetBranchAddress("Q2", &Q2);
  Ts->SetBranchAddress("Pt", &Pt);
  Ts->SetBranchAddress("E1stat", &stat);
  Ts->SetBranchAddress("fp", &fp);
  FILE * file = fopen(csvfile, "w");
  fprintf(file, "i,Ebeam,x,y,z,Q2,pT,obs,value,stat,systrel,systabs,target,hadron,Experiment\n");
  for (int i = 0; i < Ts->GetEntries(); i++){
    std::cout << i << std::endl;
    Ts->GetEntry(i);
    systrel = 0.0;
    systabs = 0.0;
    systrel += pow(0.03, 2);//target polarization
    systrel += pow(0.05, 2);//dilution
    systrel += pow(0.025, 2);//radiative correction
    systrel += pow(0.03, 2);//diffractive meson
    systrel += pow(0.002, 2);//random coincidence
    if (Ebeam > 10.0)//raw asymmetry
      systabs += 7.78e-4 / 0.7 / fp;
    else
      systabs += 1.1e-3 / 0.7 / fp;
    systrel = sqrt(systrel);
    if (Hadron == 0)
      fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, "AUT", 0.0, stat, systrel, systabs, "proton", "pi+", "solid");
    else if (Hadron == 1)
      fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, "AUT", 0.0, stat, systrel, systabs, "proton", "pi-", "solid");
  }
  fclose(file);
  return 0;
}
