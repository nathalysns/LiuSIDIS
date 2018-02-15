#include "Lsidis.h"

using namespace std;

int main(int argc, char * argv[]){


  gRandom->SetSeed(1);

  double M_e = 0.511e-3;
  double M_p = 0.938272;
  
  Lsidis mysidis;
  TLorentzVector l;
  TLorentzVector P;
  l.SetXYZM(0.0, 0.0, -5.0, M_e);
  P.SetXYZM(100.0 * sin(-0.05), 0.0, 100.0 * cos(-0.05), M_p);

  double lumi = 1.0e+10 * pow(0.197327, 2);//10^36 cm^-2 s^-1

  mysidis.SetNucleus(1, 0);//helium-3

  mysidis.SetHadron("pi+");
  int ic = mysidis.GetHadronCharge();
  int pid = mysidis.GetHadronID();

  mysidis.SetInitialState(l, P);
  mysidis.SetPDFset("CJ15lo");
  mysidis.SetFFset("DSSFFlo");
  
  double Xmin[6] = {0.03, 0.9, 0.02, 0.0, -M_PI, -M_PI};//x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.999, 10.0, 0.999, 1.8, M_PI, M_PI};

  mysidis.SetRange(Xmin, Xmax);

  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  double vz = 0.0;
  
  double Ngen = 5.0e3;

  int Nrec = 0;
 
  FILE * file = fopen("eictest.dat", "w");

  for (int i = 0; i < Ngen; i++){
    //vz = gRandom->Uniform(-370, -330);
    //weight = mysidis.GenerateEvent(0, 1);//please check the consistency with Xmin
    weight = mysidis.GibbsSampler(0, 1);
    if (weight > 0){
      lp = mysidis.GetLorentzVector("lp");
      Ph = mysidis.GetLorentzVector("Ph");
      //number of particles, x, z, Q2, Beam polarization, Pt, phih, ...
      fprintf(file, " %d    %.1f    %.1f    %.1f    %.2E    %.1f    %.1f    %.4E    %.4E    %.4E\n",
	      2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, weight, weight*lumi, weight*lumi/Ngen);
      fprintf(file, " %d    %d    %d    %d    %d    %d    %.4E    %.4E    %.4E    %.4E    %.4E    %.4E    %.4E    %.4E\n",
	      1, -1, 1, 11, 0, 0, lp.X(), lp.Y(), lp.Z(), lp.E(), 0.000511, 0.0, 0.0, vz);
      fprintf(file, " %d    %d    %d    %d    %d    %d    %.4E    %.4E    %.4E    %.4E    %.4E    %.4E    %.4E    %.4E\n",
	      2, ic, 1, pid, 0, 0, Ph.X(), Ph.Y(), Ph.Z(), Ph.E(), Ph.M(), 0.0, 0.0, vz);
      Nrec++;
    }
  }
  fclose(file);

  cout << Nrec << endl;
    
  return 0;
}
