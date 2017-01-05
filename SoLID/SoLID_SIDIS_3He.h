#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"

#include "Lsidis.h"

// Acceptance 
TFile * file_e = new TFile("Acceptance/acceptance_solid_SIDIS_He3_electron_201701_1e7_output.root", "r");
TFile * file_pi = new TFile("Acceptance/acceptance_solid_SIDIS_He3_pim_201701_1e7_output.root", "r");
TH2F * acc_FA_e = (TH2F *) file_e->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_e = (TH2F *) file_e->Get("acceptance_ThetaP_largeangle");
TH2F * acc_FA_pi = (TH2F *) file_pi->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_pi = (TH2F *) file_pi->Get("acceptance_ThetaP_largeangle");

double GetAcceptance_e(const TLorentzVector p){//Get electron acceptance
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < 8.0 || theta > 30.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(mom));
  if (mom > 3.5)
    acc += acc_LA_e->GetBinContent(acc_LA_e->GetXaxis()->FindBin(theta), acc_LA_e->GetYaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_pi(const TLorentzVector p){//Get pion acceptance
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < 8.0 || theta > 18.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_pi->GetBinContent(acc_FA_pi->GetXaxis()->FindBin(theta), acc_FA_pi->GetYaxis()->FindBin(mom));
  return acc;
}

int GenerateBinInfoFile(const char * filename, const double Ebeam, const char * hadron){//Bin the data and create the bin info file
  FILE * fp = fopen(filename, "w");
  fprintf(fp, "Q2l\t Q2u\t zl\t zu\t Ptl\t Ptu\t xl\t xu\n");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2, 1);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+10 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 48.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 21.0 * 24.0 * 3600.0;
  double Nsim = 1.0e6;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double Rfactor = 0;
  double acc = 0;
  int Nx = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  double Q2list[7] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
  double statlist[6] = {2.0e7, 1.2e7, 6.0e6, 3.0e6, 2.0e6, 2.0e6};
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
	    lp = sidis.GetLorentzVector("lp");
	    Ph = sidis.GetLorentzVector("Ph");
	    acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph);
	    if (acc > 0){
	      sidis.CalculateRfactor();
	      sidis.GetVariable("Rfactor");
	      if (true)
		hx->Fill(sidis.GetVariable("x"), weight * acc);
	    }
	  }
	}
	hx->Scale(lumi * time * eff / Nsim);
	if ((hx->Integral(1, -1) < statlist[Qi] && kj < 6) || (hx->Integral(1, -1) < 0.25 * statlist[Qi] && kj == 6)){
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
