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
#include "TRandom3.h"

#include "structure.h"

using namespace std;

double Vars[3000][5], Values[3000], Errors[3000];
TString Targets[3000], Hadrons[3000];
int Npoints = 0;

int LoadData(const char * file){
  ifstream infile(file);
  char tmp[300];
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> Errors[Npoints]){
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int LoadSmear(const char * file){
  ifstream infile(file);
  char tmp[300];
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> Errors[Npoints]){
    Values[Npoints] += gRandom->Gaus(0.0, Errors[Npoints]);
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int SimulateData(const double * par, const char * infile, const char * outfile){
  double var[5], value, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> error){
    value = AUTSivers(var, target, hadron, par);
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTsivers", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

int SimulateSoLID(const double * par, const char * infile, const char * outfile){
  double var[5], value, stat, systrel, systabs, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> stat >> systrel >> systabs){
    value = AUTSivers(var, target, hadron, par);
    error = sqrt(pow(stat, 2) + pow(systabs, 2) + pow(value * systrel, 2));
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

int SimulateWorld(const double * par, const char * infile, const char * outfile){
  double var[5], value, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> error){
    value = AUTSivers(var, target, hadron, par);
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

int SimulateEIC(const double * par, const char * infile, const char * outfile){
  double var[5], value, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> error){
    value = AUTSivers(var, target, hadron, par);
    if (strcmp(target, "proton") == 0) error *= sqrt(2.0);
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

double chi2(const double * par){
  double theory;
  double sum = 0.0;
  for (int i = 0; i < Npoints; i++){
    theory = AUTSivers(Vars[i], Targets[i].Data(), Hadrons[i].Data(), par);
    sum += pow( (theory - Values[i]) / Errors[i], 2);
  }
  return sum;
}

int Minimize(const double * init, double * cent){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(10000000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 11);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "Nu", init[0], 1e-8, -1.0, 1.0);
  min->SetLimitedVariable(1, "au", init[1], 1e-8, 0.0, 3.0);
  min->SetLimitedVariable(2, "bu", init[2], 1e-8, 0.0, 10.0);
  min->SetVariable(3, "cu", init[3], 1e-8);
  min->SetLimitedVariable(4, "Nd", init[4], 1e-8, -1.0, 1.0);
  min->SetLimitedVariable(5, "ad", init[5], 1e-8, 0.0, 3.0);
  min->SetLimitedVariable(6, "bd", init[6], 1e-8, 0.0, 10.0);
  min->SetVariable(7, "cd", init[7], 1e-8);
  min->SetFixedVariable(8, "Nub", init[8]);
  min->SetFixedVariable(9, "Ndb", init[9]);
  min->SetLimitedVariable(10, "kt2", init[10], 1e-8, 0.1, 1.0);
  min->Minimize();
  min->PrintResults();
  for (int i = 0; i < 11; i++)
    cent[i] = min->X()[i];
  return 0;
}

int Minimize2(const double * init, double * cent){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(10000000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 11);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "Nu", init[0], 1e-8, -1.0, 1.0);
  min->SetLimitedVariable(1, "au", init[1], 1e-8, 0.0, 3.0);
  min->SetLimitedVariable(2, "bu", init[2], 1e-8, 0.0, 10.0);
  min->SetFixedVariable(3, "cu", init[3]);
  min->SetLimitedVariable(4, "Nd", init[4], 1e-8, -1.0, 1.0);
  min->SetLimitedVariable(5, "ad", init[5], 1e-8, 0.0, 3.0);
  min->SetLimitedVariable(6, "bd", init[6], 1e-8, 0.0, 10.0);
  min->SetFixedVariable(7, "cd", init[7]);
  min->SetFixedVariable(8, "Nub", init[8]);
  min->SetFixedVariable(9, "Ndb", init[9]);
  min->SetLimitedVariable(10, "kt2", init[10], 1e-8, 0.16, 1.0);
  min->Minimize();
  min->PrintResults();
  for (int i = 0; i < 11; i++)
    cent[i] = min->X()[i];
  return 0;
}

int Simulation(const double * input){
  SimulateWorld(input, "../expdata/sivers00.dat", "simulate/world00.dat");
  SimulateWorld(input, "../expdata/sivers01.dat", "simulate/world01.dat");
  SimulateWorld(input, "../expdata/sivers02.dat", "simulate/world02.dat");
  SimulateWorld(input, "../expdata/sivers03.dat", "simulate/world03.dat");
  SimulateWorld(input, "../expdata/sivers04.dat", "simulate/world04.dat");
  SimulateWorld(input, "../expdata/sivers05.dat", "simulate/world05.dat");
  SimulateWorld(input, "../expdata/sivers06.dat", "simulate/world06.dat");
  SimulateWorld(input, "../expdata/sivers07.dat", "simulate/world07.dat");
  SimulateWorld(input, "../expdata/sivers08.dat", "simulate/world08.dat");
  SimulateWorld(input, "../expdata/sivers09.dat", "simulate/world09.dat");
  SimulateWorld(input, "../expdata/sivers10.dat", "simulate/world10.dat");
  SimulateWorld(input, "../expdata/sivers11.dat", "simulate/world11.dat");
  SimulateWorld(input, "../expdata/sivers12.dat", "simulate/world12.dat");
  SimulateWorld(input, "../expdata/sivers13.dat", "simulate/world13.dat");
  SimulateWorld(input, "../expdata/sivers14.dat", "simulate/world14.dat");
  SimulateWorld(input, "../expdata/sivers15.dat", "simulate/world15.dat");
  SimulateWorld(input, "../expdata/sivers16.dat", "simulate/world16.dat");
  SimulateWorld(input, "../expdata/sivers17.dat", "simulate/world17.dat");
  SimulateWorld(input, "../expdata/sivers18.dat", "simulate/world18.dat");
  SimulateWorld(input, "../expdata/sivers19.dat", "simulate/world19.dat");
  SimulateSoLID(input, "../expdata/solid01.dat", "simulate/solid01.dat");
  SimulateSoLID(input, "../expdata/solid02.dat", "simulate/solid02.dat");
  SimulateSoLID(input, "../expdata/solid03.dat", "simulate/solid03.dat");
  SimulateSoLID(input, "../expdata/solid04.dat", "simulate/solid04.dat");
  SimulateSoLID(input, "../expdata/base01.dat", "simulate/base01.dat");
  SimulateSoLID(input, "../expdata/base02.dat", "simulate/base02.dat");
  SimulateSoLID(input, "../expdata/base03.dat", "simulate/base03.dat");
  SimulateSoLID(input, "../expdata/base04.dat", "simulate/base04.dat");
  SimulateEIC(input, "../expdata/prot_pip.dat", "simulate/eic01.dat");
  SimulateEIC(input, "../expdata/prot_pim.dat", "simulate/eic02.dat");
  SimulateEIC(input, "../expdata/neut_pip.dat", "simulate/eic03.dat");
  SimulateEIC(input, "../expdata/neut_pim.dat", "simulate/eic04.dat");
  return 0;
}

int LoadWorld(){
  Npoints= 0;
  LoadData("../expdata/sivers00.dat");
  LoadData("../expdata/sivers01.dat");
  LoadData("../expdata/sivers02.dat");
  LoadData("../expdata/sivers03.dat");
  LoadData("../expdata/sivers04.dat");
  LoadData("../expdata/sivers05.dat");
  LoadData("../expdata/sivers06.dat");
  LoadData("../expdata/sivers07.dat");
  LoadData("../expdata/sivers08.dat");
  LoadData("../expdata/sivers09.dat");
  LoadData("../expdata/sivers10.dat");
  LoadData("../expdata/sivers11.dat");
  LoadData("../expdata/sivers12.dat");
  LoadData("../expdata/sivers13.dat");
  LoadData("../expdata/sivers14.dat");
  LoadData("../expdata/sivers15.dat");
  LoadData("../expdata/sivers16.dat");
  LoadData("../expdata/sivers17.dat");
  LoadData("../expdata/sivers18.dat");
  LoadData("../expdata/sivers19.dat");
  return 0;
}

int FitWorld(const double * init, double * cent){
  Npoints = 0;
  LoadSmear("simulate/world00.dat");
  LoadSmear("simulate/world01.dat");
  LoadSmear("simulate/world02.dat");
  LoadSmear("simulate/world03.dat");
  LoadSmear("simulate/world04.dat");
  LoadSmear("simulate/world05.dat");
  LoadSmear("simulate/world06.dat");
  LoadSmear("simulate/world07.dat");
  LoadSmear("simulate/world08.dat");
  LoadSmear("simulate/world09.dat");
  LoadSmear("simulate/world10.dat");
  LoadSmear("simulate/world11.dat");
  LoadSmear("simulate/world12.dat");
  LoadSmear("simulate/world13.dat");
  LoadSmear("simulate/world14.dat");
  LoadSmear("simulate/world15.dat");
  LoadSmear("simulate/world16.dat");
  LoadSmear("simulate/world17.dat");
  LoadSmear("simulate/world18.dat");
  LoadSmear("simulate/world19.dat");
  Minimize2(init, cent);
  return 0;
}

int FitSoLID(const double * init, double * cent){
  Npoints = 0;
  LoadSmear("simulate/world00.dat");
  LoadSmear("simulate/world01.dat");
  LoadSmear("simulate/world02.dat");
  LoadSmear("simulate/world03.dat");
  LoadSmear("simulate/world04.dat");
  LoadSmear("simulate/world05.dat");
  LoadSmear("simulate/world06.dat");
  LoadSmear("simulate/world07.dat");
  LoadSmear("simulate/world08.dat");
  LoadSmear("simulate/world09.dat");
  LoadSmear("simulate/world10.dat");
  LoadSmear("simulate/world11.dat");
  LoadSmear("simulate/world12.dat");
  LoadSmear("simulate/world13.dat");
  LoadSmear("simulate/world14.dat");
  LoadSmear("simulate/world15.dat");
  LoadSmear("simulate/world16.dat");
  LoadSmear("simulate/world17.dat");
  LoadSmear("simulate/world18.dat");
  LoadSmear("simulate/world19.dat");
  LoadSmear("simulate/solid01.dat");
  LoadSmear("simulate/solid02.dat");
  LoadSmear("simulate/solid03.dat");
  LoadSmear("simulate/solid04.dat");
  Minimize(init, cent);
  return 0;
}

int FitBase(const double * init, double * cent){
  Npoints = 0;
  LoadSmear("simulate/world00.dat");
  LoadSmear("simulate/world01.dat");
  LoadSmear("simulate/world02.dat");
  LoadSmear("simulate/world03.dat");
  LoadSmear("simulate/world04.dat");
  LoadSmear("simulate/world05.dat");
  LoadSmear("simulate/world06.dat");
  LoadSmear("simulate/world07.dat");
  LoadSmear("simulate/world08.dat");
  LoadSmear("simulate/world09.dat");
  LoadSmear("simulate/world10.dat");
  LoadSmear("simulate/world11.dat");
  LoadSmear("simulate/world12.dat");
  LoadSmear("simulate/world13.dat");
  LoadSmear("simulate/world14.dat");
  LoadSmear("simulate/world15.dat");
  LoadSmear("simulate/world16.dat");
  LoadSmear("simulate/world17.dat");
  LoadSmear("simulate/world18.dat");
  LoadSmear("simulate/world19.dat");
  LoadSmear("simulate/base01.dat");
  LoadSmear("simulate/base02.dat");
  LoadSmear("simulate/base03.dat");
  LoadSmear("simulate/base04.dat");
  Minimize(init, cent);
  return 0;
}

int FitEIC(const double * init, double * cent){
  Npoints = 0;
  LoadSmear("simulate/world00.dat");
  LoadSmear("simulate/world01.dat");
  LoadSmear("simulate/world02.dat");
  LoadSmear("simulate/world03.dat");
  LoadSmear("simulate/world04.dat");
  LoadSmear("simulate/world05.dat");
  LoadSmear("simulate/world06.dat");
  LoadSmear("simulate/world07.dat");
  LoadSmear("simulate/world08.dat");
  LoadSmear("simulate/world09.dat");
  LoadSmear("simulate/world10.dat");
  LoadSmear("simulate/world11.dat");
  LoadSmear("simulate/world12.dat");
  LoadSmear("simulate/world13.dat");
  LoadSmear("simulate/world14.dat");
  LoadSmear("simulate/world15.dat");
  LoadSmear("simulate/world16.dat");
  LoadSmear("simulate/world17.dat");
  LoadSmear("simulate/world18.dat");
  LoadSmear("simulate/world19.dat");
  LoadSmear("simulate/eic01.dat");
  LoadSmear("simulate/eic02.dat");
  LoadSmear("simulate/eic03.dat");
  LoadSmear("simulate/eic04.dat");
  Minimize(init, cent);
  return 0;
}

int GetGrid(const double * cent, const char * filename, const char * outfile){
  double par[100][11];
  double xu, dxu, xd, dxd, xub, dxub, xdb, dxdb;
  ifstream infile(filename);
  for (int i = 0; i < 100; i++){
    infile >> par[i][0] >> par[i][1] >> par[i][2] >> par[i][3] >> par[i][4] >> par[i][5] >> par[i][6] >> par[i][7] >> par[i][8] >> par[i][9] >> par[i][10];
  }
  infile.close();
  double Q2 = 2.4;
  FILE * fs = fopen(outfile, "w");
  fprintf(fs, "x\txu\tdxu\txd\tdxd\txub\tdxub\txdb\tdxdb\n");
  double x = 0;
  for (int i = 0; i < 1000; i++){
    x = pow(10.0, - 4.0 + 4.0 / 1000 * (i+1));
    xu = x * f1T1(2, x, Q2, "proton", cent);
    xd = x * f1T1(1, x, Q2, "proton", cent);
    xub = x * f1T1(-2, x, Q2, "proton", cent);
    xdb = x * f1T1(-1, x, Q2, "proton", cent);
    dxu = 0;
    dxd = 0;
    dxub = 0;
    dxdb = 0;
    for (int j = 0; j < 100; j++){
      dxu += pow(x * f1T1(2, x, Q2, "proton", par[j]) - xu, 2);
      dxd += pow(x * f1T1(1, x, Q2, "proton", par[j]) - xd, 2);
      dxub += pow(x * f1T1(-2, x, Q2, "proton", par[j]) - xub, 2);
      dxdb += pow(x * f1T1(-1, x, Q2, "proton", par[j]) - xdb, 2);
    }
    dxu = sqrt(dxu/100);
    dxd = sqrt(dxd/100);
    dxub = sqrt(dxub/100);
    dxdb = sqrt(dxdb/100);
    fprintf(fs, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	    x, xu, dxu, xd, dxd, xub, dxub, xdb, dxdb);
  }
  fclose(fs);  
  return 0;
}
    
  

int main(const int argc, const char * argv[]){

  if (argc < 2){
    return 0;
  }

  const int opt = atoi(argv[1]);
  LHAPDF::setVerbosity(0);

  gRandom->SetSeed(0);

  double input[11] = {-0.1, 1.0, 5.0, 0.0, 0.1, 1.0, 5.0, 0.0, 0.0, 0.0, 0.25};

  if (opt == -1){
    double cent[11];
    LoadWorld();
    Minimize2(input, cent);
    Simulation(cent);
  }

  if (opt == 0){
    Simulation(input);
  }

  if (opt == 1){
    double cent[11];
    FILE * file = fopen("out-world.dat", "w");
    for (int i = 0; i < 100; i++){
      FitWorld(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5], cent[6], cent[7], cent[8], cent[9], cent[10]);
    }
    fclose(file);
  }

  if (opt == 2){
    double cent[11];
    double input2[11];
    LoadWorld();
    Minimize2(input, input2);
    FILE * file = fopen("out-solid.dat", "w");
    for (int i = 0; i < 100; i++){
      FitSoLID(input2, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5], cent[6], cent[7], cent[8], cent[9], cent[10]);
    }
    fclose(file);
  }

  if (opt == 3){
    double cent[11];
    FILE * file = fopen("out-base.dat", "w");
    for (int i = 0; i < 100; i++){
      FitBase(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5], cent[6], cent[7], cent[8], cent[9], cent[10]);
    }
    fclose(file);
  }

  if (opt == 4){
    double cent[11];
    FILE * file = fopen("out-eic.dat", "w");
    for (int i = 0; i < 100; i++){
      FitEIC(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5], cent[6], cent[7], cent[8], cent[9], cent[10]);
    }
    fclose(file);
  }

  if (opt == 5){
    double cent[11];
    LoadWorld();
    Minimize2(input, cent);
    GetGrid(cent, "out-world.dat", "f1t1-world.dat");
    GetGrid(cent, "out-solid.dat", "f1t1-solid.dat");
    GetGrid(cent, "out-base.dat", "f1t1-base.dat");
    GetGrid(cent, "out-eic.dat", "f1t1-eic.dat");
  }
    
  
  
  return 0;
}

