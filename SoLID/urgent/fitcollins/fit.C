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
  double var[5], value, stat, systrel, systabs, error, R;
  char target[10], hadron[5], tmp[10];
  double MiT2 = 0.5;
  double MfT2 = 0.5;
  double kT2 = 0.5;
  double Mp = 0.938272;
  double Mh = 0.13957;
  double Q2, xn, z, Pt, yi, yf, yh, Rf, Ri;
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> stat >> systrel >> systabs){
    Q2 = var[0];
    z = var[3];
    xn = 2.0 * var[1] / (1.0 + sqrt(1.0 + 4.0 * pow(var[1] * Mp, 2) / Q2));
    Pt = var[4];
    yi = 0.5 * log(Q2 / MiT2);
    yf = -0.5 * log(Q2 / MfT2);
    yh = log( (sqrt(Q2) * z * (Q2 - xn * xn * Mp * Mp)) / ( 2.0 * xn * xn * Mp * Mp * sqrt(Mh * Mh + Pt * Pt)) - sqrt(Q2) / (xn * Mp) * sqrt(pow(z * (Q2 - xn * xn * Mp * Mp), 2) / (4.0 * xn * xn  * Mp * Mp * (Mh * Mh + Pt * Pt)) - 1.0));
    Rf = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MfT2) * (exp(yf - yh) + exp(yh - yf)) - sqrt(kT2) * Pt;
    Ri = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MiT2) * (exp(yi - yh) - exp(yh - yi)) - sqrt(kT2) * Pt;
    R = abs(Rf / Ri);
    //if (R > 0.4) continue;
    value = AUTCollins(var, target, hadron, par);
    error = stat;
    //error = sqrt(pow(stat, 2) + pow(systabs, 2) + pow(value * systrel, 2));
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTcollins", target, hadron, value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

int SimulateSoLIDsyst(const double * par, const char * infile, const char * outfile){
  double var[5], value, stat, systrel, systabs, error, R;
  char target[10], hadron[5], tmp[10];
  double MiT2 = 0.5;
  double MfT2 = 0.5;
  double kT2 = 0.5;
  double Mp = 0.938272;
  double Mh = 0.13957;
  double Q2, xn, z, Pt, yi, yf, yh, Rf, Ri;
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> stat >> systrel >> systabs){
    Q2 = var[0];
    z = var[3];
    xn = 2.0 * var[1] / (1.0 + sqrt(1.0 + 4.0 * pow(var[1] * Mp, 2) / Q2));
    Pt = var[4];
    yi = 0.5 * log(Q2 / MiT2);
    yf = -0.5 * log(Q2 / MfT2);
    yh = log( (sqrt(Q2) * z * (Q2 - xn * xn * Mp * Mp)) / ( 2.0 * xn * xn * Mp * Mp * sqrt(Mh * Mh + Pt * Pt)) - sqrt(Q2) / (xn * Mp) * sqrt(pow(z * (Q2 - xn * xn * Mp * Mp), 2) / (4.0 * xn * xn  * Mp * Mp * (Mh * Mh + Pt * Pt)) - 1.0));
    Rf = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MfT2) * (exp(yf - yh) + exp(yh - yf)) - sqrt(kT2) * Pt;
    Ri = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MiT2) * (exp(yi - yh) - exp(yh - yi)) - sqrt(kT2) * Pt;
    R = abs(Rf / Ri);
    //if (R > 0.4) continue;
    value = AUTCollins(var, target, hadron, par);
    //error = stat;
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
  //min->PrintResults();
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
  SimulateSoLIDsyst(input, "../expdata/solid01.dat", "simulate/solidsyst01.dat");
  SimulateSoLIDsyst(input, "../expdata/solid02.dat", "simulate/solidsyst02.dat");
  SimulateSoLIDsyst(input, "../expdata/solid03.dat", "simulate/solidsyst03.dat");
  SimulateSoLIDsyst(input, "../expdata/solid04.dat", "simulate/solidsyst04.dat");
  SimulateSoLID(input, "../expdata/base0-neutron-pip.dat", "simulate/base0-neutron-pip.dat");
  SimulateSoLID(input, "../expdata/base0-neutron-pim.dat", "simulate/base0-neutron-pim.dat");
  SimulateSoLID(input, "../expdata/base0-proton-pip.dat", "simulate/base0-proton-pip.dat");
  SimulateSoLID(input, "../expdata/base0-proton-pim.dat", "simulate/base0-proton-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base0-neutron-pip.dat", "simulate/base0syst-neutron-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base0-neutron-pim.dat", "simulate/base0syst-neutron-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base0-proton-pip.dat", "simulate/base0syst-proton-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base0-proton-pim.dat", "simulate/base0syst-proton-pim.dat");
  SimulateSoLID(input, "../expdata/base1-neutron-pip.dat", "simulate/base1-neutron-pip.dat");
  SimulateSoLID(input, "../expdata/base1-neutron-pim.dat", "simulate/base1-neutron-pim.dat");
  SimulateSoLID(input, "../expdata/base1-proton-pip.dat", "simulate/base1-proton-pip.dat");
  SimulateSoLID(input, "../expdata/base1-proton-pim.dat", "simulate/base1-proton-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base1-neutron-pip.dat", "simulate/base1syst-neutron-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base1-neutron-pim.dat", "simulate/base1syst-neutron-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base1-proton-pip.dat", "simulate/base1syst-proton-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base1-proton-pim.dat", "simulate/base1syst-proton-pim.dat");
  SimulateSoLID(input, "../expdata/base2-neutron-pip.dat", "simulate/base2-neutron-pip.dat");
  SimulateSoLID(input, "../expdata/base2-neutron-pim.dat", "simulate/base2-neutron-pim.dat");
  SimulateSoLID(input, "../expdata/base2-proton-pip.dat", "simulate/base2-proton-pip.dat");
  SimulateSoLID(input, "../expdata/base2-proton-pim.dat", "simulate/base2-proton-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base2-neutron-pip.dat", "simulate/base2syst-neutron-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base2-neutron-pim.dat", "simulate/base2syst-neutron-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base2-proton-pip.dat", "simulate/base2syst-proton-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base2-proton-pim.dat", "simulate/base2syst-proton-pim.dat");
  SimulateSoLID(input, "../expdata/base3-neutron-pip.dat", "simulate/base3-neutron-pip.dat");
  SimulateSoLID(input, "../expdata/base3-neutron-pim.dat", "simulate/base3-neutron-pim.dat");
  SimulateSoLID(input, "../expdata/base3-proton-pip.dat", "simulate/base3-proton-pip.dat");
  SimulateSoLID(input, "../expdata/base3-proton-pim.dat", "simulate/base3-proton-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base3-neutron-pip.dat", "simulate/base3syst-neutron-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base3-neutron-pim.dat", "simulate/base3syst-neutron-pim.dat");
  SimulateSoLIDsyst(input, "../expdata/base3-proton-pip.dat", "simulate/base3syst-proton-pip.dat");
  SimulateSoLIDsyst(input, "../expdata/base3-proton-pim.dat", "simulate/base3syst-proton-pim.dat");
  SimulateCOMPASS(input, "../expdata/clas01.dat", "simulate/clas01.dat");
  SimulateCOMPASS(input, "../expdata/clas02.dat", "simulate/clas02.dat");
  SimulateCOMPASS(input, "../expdata/sbs01.dat", "simulate/sbs01.dat");
  SimulateCOMPASS(input, "../expdata/sbs02.dat", "simulate/sbs02.dat");
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

int FitSoLIDsyst(const double * init, double * cent){
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
  LoadSmear("simulate/solidsyst01.dat");
  LoadSmear("simulate/solidsyst02.dat");
  LoadSmear("simulate/solidsyst03.dat");
  LoadSmear("simulate/solidsyst04.dat");
  Minimize(init, cent);
  return 0;
}

int FitBase(const double * init, double * cent, const int option = 0){
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
  if (option == 0){
    LoadSmear("simulate/base0-neutron-pip.dat");
    LoadSmear("simulate/base0-neutron-pim.dat");
    LoadSmear("simulate/base0-proton-pip.dat");
    LoadSmear("simulate/base0-proton-pim.dat");
  }
  else if (option == 1){
    LoadSmear("simulate/base1-neutron-pip.dat");
    LoadSmear("simulate/base1-neutron-pim.dat");
    LoadSmear("simulate/base1-proton-pip.dat");
    LoadSmear("simulate/base1-proton-pim.dat");
  }
  else if (option == 2){
    LoadSmear("simulate/base2-neutron-pip.dat");
    LoadSmear("simulate/base2-neutron-pim.dat");
    LoadSmear("simulate/base2-proton-pip.dat");
    LoadSmear("simulate/base2-proton-pim.dat");
  }
  else if (option == 3){
    LoadSmear("simulate/base3-neutron-pip.dat");
    LoadSmear("simulate/base3-neutron-pim.dat");
    LoadSmear("simulate/base3-proton-pip.dat");
    LoadSmear("simulate/base3-proton-pim.dat");
  }
  Minimize(init, cent);
  return 0;
}

int FitBasesyst(const double * init, double * cent, const int option = 0){
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
  if (option == 0){
    LoadSmear("simulate/base0syst-neutron-pip.dat");
    LoadSmear("simulate/base0syst-neutron-pim.dat");
    LoadSmear("simulate/base0syst-proton-pip.dat");
    LoadSmear("simulate/base0syst-proton-pim.dat");
  }
  else if (option == 1){
    LoadSmear("simulate/base1syst-neutron-pip.dat");
    LoadSmear("simulate/base1syst-neutron-pim.dat");
    LoadSmear("simulate/base1syst-proton-pip.dat");
    LoadSmear("simulate/base1syst-proton-pim.dat");
  }
  else if (option == 2){
    LoadSmear("simulate/base2syst-neutron-pip.dat");
    LoadSmear("simulate/base2syst-neutron-pim.dat");
    LoadSmear("simulate/base2syst-proton-pip.dat");
    LoadSmear("simulate/base2syst-proton-pim.dat");
  }
  else if (option == 3){
    LoadSmear("simulate/base3syst-neutron-pip.dat");
    LoadSmear("simulate/base3syst-neutron-pim.dat");
    LoadSmear("simulate/base3syst-proton-pip.dat");
    LoadSmear("simulate/base3syst-proton-pim.dat");
  }
  Minimize(init, cent);
  return 0;
}

int FitOthers(const double * init, double * cent){
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
  LoadSmear("simulate/clas01.dat");
  LoadSmear("simulate/clas02.dat");
  LoadSmear("simulate/sbs01.dat");
  LoadSmear("simulate/sbs02.dat");
  Minimize(init, cent);
  return 0;
}

int FitAll(const double * init, double * cent){
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
  LoadSmear("simulate/solidsyst01.dat");
  LoadSmear("simulate/solidsyst02.dat");
  LoadSmear("simulate/solidsyst03.dat");
  LoadSmear("simulate/solidsyst04.dat");
  LoadSmear("simulate/clas01.dat");
  LoadSmear("simulate/clas02.dat");
  LoadSmear("simulate/sbs01.dat");
  LoadSmear("simulate/sbs02.dat");
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

int GetTensorCharge(double * input, const char * filename){
  const int Npt = 1000;
  double X[Npt];
  double step = 1.0 / Npt;
  for (int i = 0; i < Npt; i++)
    X[i] = (i + 0.5) * step;

  return 0;
}


  

  

int main(const int argc, const char * argv[]){

  if (argc < 2){
    return 0;
  }

  const int opt = atoi(argv[1]);
  LHAPDF::setVerbosity(0);

  gRandom->SetSeed(0);

  double input[6] = {0.4, -0.45, 1.0, 3.0, 0., 0.25};

  if (opt == 0){
    Simulation(input);
  }

  if (opt == 1){//fit world
    double cent[6];
    FILE * file = fopen("out-world.dat", "w");
    for (int i = 0; i < 100; i++){
      FitWorld(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 20){//fit solid 
    double cent[6];
    FILE * file = fopen("out-solid.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitSoLID(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 21){//fit solidsyst
    double cent[6];
    FILE * file = fopen("out-solidsyst.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitSoLIDsyst(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 300){//fit base
    double cent[6];
    FILE * file = fopen("out-base0.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBase(input, cent, 0);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 301){//fit base
    double cent[6];
    FILE * file = fopen("out-base1.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBase(input, cent, 1);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 302){//fit base
    double cent[6];
    FILE * file = fopen("out-base2.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBase(input, cent, 2);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 303){//fit base
    double cent[6];
    FILE * file = fopen("out-base3.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBase(input, cent, 3);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 310){//fit basesyst
    double cent[6];
    FILE * file = fopen("out-base0syst.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBasesyst(input, cent, 0);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 311){//fit basesyst
    double cent[6];
    FILE * file = fopen("out-base1syst.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBasesyst(input, cent, 1);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 312){//fit basesyst
    double cent[6];
    FILE * file = fopen("out-base2syst.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBasesyst(input, cent, 2);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  if (opt == 313){//fit basesyst
    double cent[6];
    FILE * file = fopen("out-base3syst.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitBasesyst(input, cent, 3);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  

  if (opt == 4){// fit others
    double cent[6];
    FILE * file = fopen("out-others.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitOthers(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }

  if (opt == 5){// fit all
    double cent[6];
    FILE * file = fopen("out-all.dat", "w");
    for (int i = 0; i < 100; i++){
      cout << i << endl;
      FitAll(input, cent);
      fprintf(file, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	      cent[0], cent[1], cent[2], cent[3], cent[4], cent[5]);
    }
    fclose(file);
  }
  
  if (opt == 6){//solid h1
    GetGrid(input, "out-world.dat", "h1-world.dat");
    GetGrid(input, "out-solid.dat", "h1-solid.dat");
    GetGrid(input, "out-solidsyst.dat", "h1-solidsyst.dat");
  }

  if (opt == 65){
    GetGrid(input, "out-others.dat", "h1-others.dat");
    GetGrid(input, "out-all.dat", "h1-all.dat");
  }

  if (opt == 60){//base h1
    GetGrid(input, "out-world.dat", "h1-world.dat");
    GetGrid(input, "out-base0.dat", "h1-base0.dat");
    GetGrid(input, "out-base1.dat", "h1-base1.dat");
    GetGrid(input, "out-base2.dat", "h1-base2.dat");
    GetGrid(input, "out-base3.dat", "h1-base3.dat");
  }
  if (opt == 61){//basesyst h1
    GetGrid(input, "out-world.dat", "h1-world.dat");
    GetGrid(input, "out-base0syst.dat", "h1-base0syst.dat");
    GetGrid(input, "out-base1syst.dat", "h1-base1syst.dat");
    GetGrid(input, "out-base2syst.dat", "h1-base2syst.dat");
    GetGrid(input, "out-base3syst.dat", "h1-base3syst.dat");
  }
  
  if (opt == 8){
    double u = gtu(input);
    double d = gtd(input);
    double uworld, dworld, usolid, dsolid, ubase, dbase, uothers, dothers, uall, dall;
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
    ifstream fothers("out-others.dat");
    uothers = 0; dothers = 0;
    for (int i = 0; i < 100; i++){
      fothers >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] >> par[5];
      uothers += pow(gtu(par) - u, 2);
      dothers += pow(gtd(par) - d, 2);
    }
    uothers = sqrt(uothers / 100);
    dothers = sqrt(dothers / 100);
    fothers.close();
    ifstream fall("out-all.dat");
    uall = 0; dall = 0;
    for (int i = 0; i < 100; i++){
      fall >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] >> par[5];
      uall += pow(gtu(par) - u, 2);
      dall += pow(gtd(par) - d, 2);
    }
    uall = sqrt(uall / 100);
    dall = sqrt(dall / 100);
    fall.close();
    double tol = 7.04;
    cout << u << "\t" << uworld * tol << "\t" << usolid * tol << "\t" << ubase * tol << "\t" << uothers * tol << "\t" << uall * tol << endl;
    cout << d << "\t" << dworld * tol << "\t" << dsolid * tol << "\t" << dbase * tol << "\t" << dothers * tol << "\t" << dall * tol << endl;     
  }
  
  return 0;
}

