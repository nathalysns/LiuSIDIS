/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
#include <TMatrixD.h>

const double max1 = 30.;
const double min1 = 0.8;
double qmax1 = 0.17;
double qmin1 = -0.17;

/*void plot(Int_t z_flag =1, Int_t Q2_flag=1, Int_t flag_t=0, Int_t flag=0){{{*/
void plot(TString ini_ele_energy, TString ini_ion_energy, Int_t target_flag,Int_t particle_flag, Int_t sign, Int_t z_flag=1, Int_t pt_flag=1, Int_t flag_t=0, Int_t single=0, Float_t xFactor=1, Float_t yFactor=1){
	//flag_t: 1->collins, 2->sivers

	Float_t factor,shift;
	factor=(log(max1)-log(min1))/(qmax1-qmin1);
	shift=log(max1)/factor-qmax1;

//	int target_flag = 1;
//	int particle_flag =i;
//	int z_flag = i;
//	int Q2_flag = j;


	/*Input{{{*/
	TString target = "X";
	if(target_flag==1)
		target ="A1";
	else if(target_flag==2)
		target ="d2";
	else if(target_flag==3)
		target ="A3";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString particle = "X";
	if(particle_flag==1 && sign==1)
		particle ="pip";
	else if(particle_flag==1 && sign==-1)
		particle ="pim";
	else if(particle_flag==2 && sign==1)
		particle ="kp";
	else if(particle_flag==2 && sign==-1)
		particle ="km";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString filename;

	filename.Form("../results/%s/%s_%s/%s_%d_%d.dat",target.Data(),ini_ele_energy.Data(),ini_ion_energy.Data(),particle.Data(),z_flag,pt_flag);
		cerr<<Form("../results/%s/%s_%s/%s_%d_%d.dat",target.Data(),ini_ele_energy.Data(),ini_ion_energy.Data(),particle.Data(),z_flag,pt_flag)<<endl;
	ifstream infile(filename);

	Float_t ptmin,ptmax;
	if (pt_flag==1){
		ptmin = 0.; ptmax =0.2;
	}else if (pt_flag==2){
		ptmin = 0.4; ptmax =0.6;
	}else if (pt_flag==3){
		ptmin = 0.8; ptmax =1.0;
	}

	Float_t zmin,zmax;
	if (z_flag==1){
		zmin = 0.3; zmax = 0.35;
	}else if (z_flag==2){
		zmin = 0.4; zmax = 0.45;
	}else if (z_flag==3){
		zmin = 0.5; zmax = 0.55;
	}else if (z_flag==4){
		zmin = 0.6; zmax = 0.65;
	}
	/*}}}*/

	gStyle->SetOptStat(0);
	Int_t count1,count2;
	Double_t Q2[5000],x[5000],z[5000],pt[5000],y[5000],Astat[5000],coverage[5000],coef[3][5000];
	Double_t N[5000];
	Double_t temp;
	Int_t Q2_flag = -1;
	Int_t ncount=0;

	infile >> count1;
	cerr<<" Q2_total_bin="<<count1<<endl;
	Int_t Q2_temp,x_temp;
	while(infile >> Q2_flag >> count2){
		cerr<<" Q2_flag="<<Q2_flag<<", count="<<count2<<endl;
		for (Int_t j1=0;j1<count2;j1++){
			infile >> Q2_temp >> x_temp >> z[ncount] >> Q2[ncount] >> pt[ncount] >> x[ncount] >> y[ncount] >> Astat[ncount] >> N[ncount]>>
				coverage[ncount] >> coef[0][ncount] >> coef[1][ncount] >> coef[2][ncount];
			if(Q2_temp==Q2_flag && x_temp== j1){
				Astat[ncount] *= coef[2-flag_t][ncount]; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity
				if (Astat[ncount]<0.01){
					cerr<<Form("---z=%d, pt=%d, Q2 = (%d) %f , x = (%d) %f, Asys= %f",z_flag,pt_flag, Q2_flag, Q2[ncount],j1,x[ncount], Astat[ncount])<<endl;
					// 	Astat[ncount] *= (coef[0][ncount]+coef[1][ncount]+coef[2][ncount])/3.; 
					Astat[ncount] *= factor; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity
                                        Q2[ncount]=log(Q2[ncount]);
					ncount ++;
				}
			}
		}
	}

	infile.close();
	TString titlename;

	//	t1->SetNDC();


}
	/*}}}*/
