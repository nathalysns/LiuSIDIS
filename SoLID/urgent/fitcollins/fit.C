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
  

int SimulateHERMES(const double * par, const char * infile, const char * outfile){
  double var[5], value, stat, syst, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> stat >> syst){
    value = AUTCollins(var, target, hadron, par);
    error = sqrt(pow(stat, 2) + pow(syst, 2));
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

int SimulateCOMPASS(const double * par, const char * infile, const char * outfile){
  double var[5], value, error;
  char target[10], hadron[5], tmp[10];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> error){
    value = AUTCollins(var, target, hadron, par);
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
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

double chi2(const double * par){
  double theory;
  //double par[7] = {params[0], params[2], params[3], params[1], params[2], params[3], params[4]};
  double sum = 0.0;
  for (int i = 0; i < Npoints; i++){
    theory = AUTCollins(Vars[i], Targets[i].Data(), Hadrons[i].Data(), par);
    sum += pow( (theory - Values[i]) / Errors[i], 2);
  }
  return sum;
}

int Minimize(const double * init, double * cent, double * cov = 0){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(10000000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 6);
  min->SetFunction(f);
  min->SetVariable(0, "Nu", init[0], 1e-8);
  min->SetVariable(1, "Nd", init[1], 1e-8);
  min->SetVariable(2, "a", init[2], 1e-8);
  min->SetVariable(3, "b", init[3], 1e-8);
  min->SetVariable(4, "c", init[4], 1e-8);
  min->SetVariable(5, "kt2", init[5], 1e-8);
  min->Minimize();
  min->PrintResults();
  //min->GetCovMatrix(cov);
  for (int i = 0; i < 6; i++)
    cent[i] = min->X()[i];
  return 0;
}

int Minimize2(const double * init, double * cent, double * cov = 0){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(10000000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 6);
  min->SetFunction(f);
  min->SetVariable(0, "Nu", init[0], 1e-8);
  min->SetVariable(1, "Nd", init[1], 1e-8);
  min->SetFixedVariable(2, "a", init[2]);
  min->SetVariable(3, "b", init[3], 1e-8);
  min->SetFixedVariable(4, "c", init[4]);
  min->SetVariable(5, "kt2", init[5], 1e-8);
  min->Minimize();
  min->PrintResults();
  //min->GetCovMatrix(cov);
  for (int i = 0; i < 6; i++)
    cent[i] = min->X()[i];
  return 0;
}

int Simulation(const double * input){
  SimulateHERMES(input, "../expdata/urgent00.dat", "simulate/world00.dat");
  SimulateHERMES(input, "../expdata/urgent01.dat", "simulate/world01.dat");
  SimulateHERMES(input, "../expdata/urgent02.dat", "simulate/world02.dat");
  SimulateHERMES(input, "../expdata/urgent03.dat", "simulate/world03.dat");
  SimulateHERMES(input, "../expdata/urgent04.dat", "simulate/world04.dat");
  SimulateHERMES(input, "../expdata/urgent05.dat", "simulate/world05.dat");
  SimulateCOMPASS(input, "../expdata/urgent09.dat", "simulate/world09.dat");
  SimulateCOMPASS(input, "../expdata/urgent10.dat", "simulate/world10.dat");
  SimulateCOMPASS(input, "../expdata/urgent11.dat", "simulate/world11.dat");
  SimulateCOMPASS(input, "../expdata/urgent12.dat", "simulate/world12.dat");
  SimulateCOMPASS(input, "../expdata/urgent13.dat", "simulate/world13.dat");
  SimulateCOMPASS(input, "../expdata/urgent14.dat", "simulate/world14.dat");
  SimulateCOMPASS(input, "../expdata/urgent15.dat", "simulate/world15.dat");
  SimulateCOMPASS(input, "../expdata/urgent16.dat", "simulate/world16.dat");
  SimulateCOMPASS(input, "../expdata/urgent17.dat", "simulate/world17.dat");
  SimulateCOMPASS(input, "../expdata/urgent18.dat", "simulate/world18.dat");
  SimulateCOMPASS(input, "../expdata/urgent19.dat", "simulate/world19.dat");
  SimulateCOMPASS(input, "../expdata/urgent20.dat", "simulate/world20.dat");
  SimulateSoLID(input, "../expdata/solid01.dat", "simulate/solid01.dat");
  SimulateSoLID(input, "../expdata/solid02.dat", "simulate/solid02.dat");
  SimulateSoLID(input, "../expdata/solid03.dat", "simulate/solid03.dat");
  SimulateSoLID(input, "../expdata/solid04.dat", "simulate/solid04.dat");
  SimulateSoLID(input, "../expdata/base01.dat", "simulate/base01.dat");
  SimulateSoLID(input, "../expdata/base02.dat", "simulate/base02.dat");
  SimulateSoLID(input, "../expdata/base03.dat", "simulate/base03.dat");
  SimulateSoLID(input, "../expdata/base04.dat", "simulate/base04.dat");
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
  LoadSmear("simulate/world20.dat");
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
  LoadSmear("simulate/world20.dat");
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
  LoadSmear("simulate/world20.dat");
  LoadSmear("simulate/base01.dat");
  LoadSmear("simulate/base02.dat");
  LoadSmear("simulate/base03.dat");
  LoadSmear("simulate/base04.dat");
  Minimize(init, cent);
  return 0;
}

int GetGrid(double * input, const char * filename, const char * outfile){
  double x, xhu, xhd, dxhu, dxhd;
  double par[100][6];
  ifstream file(filename);
  for (int i = 0; i < 100; i++){
    file >> par[i][0] >> par[i][1] >> par[i][2] >> par[i][3] >> par[i][4] >> par[i][5];
  }
  file.close();
  FILE * fs = fopen(outfile, "w");
  fprintf(fs, "x\txu\tdxu\txd\tdxd\n");
  for (int i = 0; i < 200; i++){
    x = 1.0 / 400 + 1.0 / 200 * i;
    xhu = h1u(x, input) * x;
    xhd = h1d(x, input) * x;
    dxhu = 0;
    dxhd = 0;
    for (int j = 0; j < 100; j++){
      dxhu += pow(h1u(x, par[j]) * x - xhu, 2);
      dxhd += pow(h1d(x, par[j]) * x - xhd, 2);
    }
    dxhu = sqrt(dxhu / 100);
    dxhd = sqrt(dxhd / 100);
    fprintf(fs, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	    x, xhu, dxhu, xhd, dxhd);
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

  double input[6] = {0.4, -0.45, 1.0, 3.0, 0.2, 0.25};

  if (opt == 0){
    Simulation(input);
  }

  if (opt == 1){
    double cent[6];
    FILE * file = fopen("out-world.dat", "w");
    for (int i = 0; i < 100; i++){
      FitWorld(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 2){
    double cent[6];
    FILE * file = fopen("out-solid.dat", "w");
    for (int i = 0; i < 100; i++){
      FitSoLID(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 3){
    double cent[6];
    FILE * file = fopen("out-base.dat", "w");
    for (int i = 0; i < 100; i++){
      FitBase(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 4){
    GetGrid(input, "out-world.dat", "h1-world.dat");
    GetGrid(input, "out-solid.dat", "h1-solid.dat");
    GetGrid(input, "out-base.dat", "h1-base.dat");
  }
   
  
  if (opt == 8){
    double u = gtu(input);
    double d = gtd(input);
    double uworld, dworld, usolid, dsolid, ubase, dbase;
    double par[6];
    ifstream fworld("out-world.dat");
    uworld = 0; dworld = 0;
    for (int i = 0; i < 100; i++){
      fworld >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] >> par[5];
      uworld += pow(gtu(par) - u, 2);
      dworld += pow(gtd(par) - d, 2);
    }
    uworld = sqrt(uworld / 100);
    dworld = sqrt(dworld / 100);
    fworld.close();
    ifstream fsolid("out-solid.dat");
    usolid = 0; dsolid = 0;
    for (int i = 0; i < 100; i++){
      fsolid >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] >> par[5];
      usolid += pow(gtu(par) - u, 2);
      dsolid += pow(gtd(par) - d, 2);
    }
    usolid = sqrt(usolid / 100);
    dsolid = sqrt(dsolid / 100);
    fsolid.close();
    ifstream fbase("out-base.dat");
    ubase = 0; dbase = 0;
    for (int i = 0; i < 100; i++){
      fbase >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] >> par[5];
      ubase += pow(gtu(par) - u, 2);
      dbase += pow(gtd(par) - d, 2);
    }
    ubase = sqrt(ubase / 100);
    dbase = sqrt(dbase / 100);
    fbase.close();
    double tol = 7.04;
    cout << u << "\t" << uworld * tol << "\t" << usolid * tol << "\t" << ubase * tol << endl;
    cout << d << "\t" << dworld * tol << "\t" << dsolid * tol << "\t" << dbase * tol << endl;     
  }
  
  return 0;
}

