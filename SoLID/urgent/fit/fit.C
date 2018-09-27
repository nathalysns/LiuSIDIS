#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"
#include "TMatrixDEigen.h"
#include "Math/Integrator.h"
#include "Math/WrappedTF1.h"

#include "structure.h"

using namespace std;

double Vars[3000][5], Values[3000], Errors[3000];
TString Targets[3000], Hadrons[3000];
int Npoints = 0;

int LoadHERMES(const char * file){
  ifstream infile(file);
  char tmp[300];
  double stat, syst;
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> stat >> syst){
    Errors[Npoints] = sqrt(stat * stat + syst * syst);
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int LoadCOMPASS(const char * file){
  ifstream infile(file);
  char tmp[300];
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> Errors[Npoints]){
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int LoadSoLID(const char * file){
  ifstream infile(file);
  char tmp[300];
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> Errors[Npoints]){
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int Simulate(const double * par, const char * infile, const char * outfile){
  double var[5], value, stat, systrel, systabs, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> stat >> systrel >> systabs){
    value = AUTCollins(var, target, hadron, par);
    error = sqrt(pow(stat, 2) + pow(systabs, 2) + pow(value * systrel, 2));
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

double chi2(const double * params){
  double theory;
  double par[7] = {params[0], params[2], params[3], params[1], params[2], params[3], params[4]};
  double sum = 0.0;
  for (int i = 0; i < Npoints; i++){
    theory = AUTCollins(Vars[i], Targets[i].Data(), Hadrons[i].Data(), par);
    sum += pow( (theory - Values[i]) / Errors[i], 2);
  }
  return sum;
}

int Minimize(const double * init, double * cent, double * cov){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-6);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 5);
  min->SetFunction(f);
  min->SetVariable(0, "Nu", init[0], 1e-8);
  min->SetLimitedVariable(2, "au", init[1], 1e-8, -1.0, 5.0);
  min->SetLimitedVariable(3, "bu", init[2], 1e-8, 0.0, 10.0);
  min->SetVariable(1, "Nd", init[3], 1e-8);
  min->SetLimitedVariable(4, "ad", init[4], 1e-8, -1.0, 5.0);
  //min->SetLimitedVariable(5, "bd", init[5], 1e-8, 0.0, 10.0);
  //min->SetLimitedVariable(6, "kt2", init[6], 1e-8, 0.0, 3.0);
  min->Minimize();
  min->PrintResults();
  min->GetCovMatrix(cov);
  for (int i = 0; i < 5; i++)
    cent[i] = min->X()[i];
  return 0;
}
    

int main(const int argc, const char * argv[]){

  if (argc < 2){
    return 0;
  }

  const int opt = atoi(argv[1]);

  LoadHERMES("../expdata/urgent00.dat");
  LoadHERMES("../expdata/urgent01.dat");
  LoadHERMES("../expdata/urgent02.dat");
  LoadHERMES("../expdata/urgent03.dat");
  LoadHERMES("../expdata/urgent04.dat");
  LoadHERMES("../expdata/urgent05.dat");
  
  // LoadCOMPASS("../expdata/urgent09.dat");
  // LoadCOMPASS("../expdata/urgent10.dat");
  // LoadCOMPASS("../expdata/urgent11.dat");
  // LoadCOMPASS("../expdata/urgent12.dat");
  // LoadCOMPASS("../expdata/urgent13.dat");
  // LoadCOMPASS("../expdata/urgent14.dat");
  // LoadCOMPASS("../expdata/urgent15.dat");
  // LoadCOMPASS("../expdata/urgent16.dat");
  // LoadCOMPASS("../expdata/urgent17.dat");
  // LoadCOMPASS("../expdata/urgent18.dat");
  // LoadCOMPASS("../expdata/urgent19.dat");
  // LoadCOMPASS("../expdata/urgent20.dat");

  if (false){
    LoadCOMPASS("../expdata/urgent-solid-01.dat");
    LoadCOMPASS("../expdata/urgent-solid-02.dat");
    LoadCOMPASS("../expdata/urgent-solid-03.dat");
    LoadCOMPASS("../expdata/urgent-solid-04.dat");
  }

  if (false){
    LoadCOMPASS("../expdata/urgent-base-01.dat");
    LoadCOMPASS("../expdata/urgent-base-02.dat");
    LoadCOMPASS("../expdata/urgent-base-03.dat");
    LoadCOMPASS("../expdata/urgent-base-04.dat");
  }

  cout << Npoints << endl;


  //double init[7] = {0.5, 1.0, 2.0, -0.5, 1.0, 2.0, 0.5};
  double init[5] = {0.5, -0.5, 1.0, 2.0, 0.5};
  if (opt == 1){
    double cent[7], cov[49];
    Minimize(init, cent, cov);
  }
  
  
  return 0;
}

