#include "structure.h"

#include "TString.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TRandom3.h"

double Vars[3000][5], Values[3000], Errors[3000];
TString Targets[3000], Hadrons[3000];
int Npoints = 0;

int LoadData(const char * filename){
  ifstream infile(filename);
  char tmp[300];
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> Errors[Npoints]){
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int LoadSmear(const char * filename){
  ifstream infile(filename);
  char tmp[300];
  infile.getline(tmp, 300);
  while (infile >> Vars[Npoints][0] >> Vars[Npoints][1] >> Vars[Npoints][2] >> Vars[Npoints][3] >> Vars[Npoints][4] >> tmp >> Targets[Npoints] >> Hadrons[Npoints] >> Values[Npoints] >> Errors[Npoints]){
    Values[Npoints] += gRandom->Gaus(0.0, Errors[Npoints]);
    Npoints = Npoints + 1;
  }
  infile.close();
  return 0;
}

int SetParameters(const double * par){
  Nuv = par[0];
  auv = par[1];
  buv = par[2];
  Ndv = par[3];
  adv = par[4];
  bdv = par[5];
  Nub = par[6];
  Ndb = par[7];
  M1 = par[8];
  return 0;
}

int SimulateData(const double * par, const char * infile, const char * outfile){
  SetParameters(par);
  double var[5], value, error;
  TString target, hadron;
  char tmp[300];
  ifstream fs(infile);
  fs.getline(tmp, 300);
  FILE * file = fopen(outfile, "w");
  fprintf(file, "# Q2\tx\ty\tz\tpT\tobs\ttarget\thadron\tvalue\terror\n");
  while (fs >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> tmp >> target >> hadron >> value >> error){
    value = AUTSivers(var, target.Data(), hadron.Data());
    fprintf(file, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%s\t%s\t%s\t%.4E\t%.4E\n",
	    var[0], var[1], var[2], var[3], var[4],
	    "AUTsivers", target.Data(), hadron.Data(), value, error);
  }
  fs.close();
  fclose(file);
  return 0;
}

double chi2(const double * par){
  SetParameters(par);
  double theory;
  double sum = 0;
  for (int i = 0; i < Npoints; i++){
    theory = AUTSivers(Vars[i], Targets[i].Data(), Hadrons[i].Data());
    sum += pow(theory - Values[i], 2) / pow(Errors[i], 2);
  }
  return sum;
}

int Minimize(const double * init, double * cent){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(10000000);
  min->SetTolerance(1e-6);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 9);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "Nuv", init[0], 1e-6, -1.0, 1.0);
  min->SetLimitedVariable(1, "auv", init[1], 1e-6, 0.0, 5.0);
  min->SetLimitedVariable(2, "buv", init[2], 1e-6, 0.0, 25.0);
  min->SetLimitedVariable(3, "Ndv", init[3], 1e-6, -1.0, 1.0);
  min->SetLimitedVariable(4, "adv", init[4], 1e-6, 0.0, 5.0);
  min->SetLimitedVariable(5, "bdv", init[5], 1e-6, 0.0, 25.0);
  min->SetLimitedVariable(6, "Nub", init[6], 1e-6, -1.0, 1.0);
  min->SetLimitedVariable(7, "Ndb", init[7], 1e-6, -1.0, 1.0);
  min->SetLimitedVariable(8, "M1", init[8], 1e-6, 0.0, 3.0);
  min->Minimize();
  min->PrintResults();
  for (int i = 0; i < 9; i++)
    cent[i] = min->X()[i];
  return 0;
}

int LoadWorld(){
  Npoints = 0;
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
  LoadData("../expdata/sivers20.dat");
  LoadData("../expdata/sivers21.dat");
  LoadData("../expdata/sivers22.dat");
  LoadData("../expdata/sivers23.dat");
  LoadData("../expdata/sivers24.dat");
  LoadData("../expdata/sivers25.dat");
  LoadData("../expdata/sivers26.dat");
  LoadData("../expdata/sivers27.dat");
  LoadData("../expdata/sivers28.dat");
  LoadData("../expdata/sivers29.dat");
  LoadData("../expdata/sivers30.dat");
  LoadData("../expdata/sivers31.dat");
  LoadData("../expdata/sivers32.dat");
  LoadData("../expdata/sivers33.dat");
  return Npoints;
}

int Simulation(const double * par){
  SimulateData(par, "../expdata/prot_pip.dat", "simulate/prop_pip.dat");
  SimulateData(par, "../expdata/prot_pim.dat", "simulate/prop_pim.dat");
  SimulateData(par, "../expdata/prot_kp.dat", "simulate/prop_kp.dat");
  SimulateData(par, "../expdata/prot_km.dat", "simulate/prop_km.dat");
  SimulateData(par, "../expdata/neut_pip.dat", "simulate/neut_pip.dat");
  SimulateData(par, "../expdata/neut_pim,.dat", "simulate/neut_pim.dat");
  SimulateData(par, "../expdata/neut_kp.dat", "simulate/neut_kp.dat");
  SimulateData(par, "../expdata/neut_km.dat", "simulate/neut_km.dat");
  return 0;
}

int LoadWorldSmear(){
  Npoints = 0;
  LoadSmear("../expdata/sivers00.dat");
  LoadSmear("../expdata/sivers01.dat");
  LoadSmear("../expdata/sivers02.dat");
  LoadSmear("../expdata/sivers03.dat");
  LoadSmear("../expdata/sivers04.dat");
  LoadSmear("../expdata/sivers05.dat");
  LoadSmear("../expdata/sivers06.dat");
  LoadSmear("../expdata/sivers07.dat");
  LoadSmear("../expdata/sivers08.dat");
  LoadSmear("../expdata/sivers09.dat");
  LoadSmear("../expdata/sivers10.dat");
  LoadSmear("../expdata/sivers11.dat");
  LoadSmear("../expdata/sivers12.dat");
  LoadSmear("../expdata/sivers13.dat");
  LoadSmear("../expdata/sivers20.dat");
  LoadSmear("../expdata/sivers21.dat");
  LoadSmear("../expdata/sivers22.dat");
  LoadSmear("../expdata/sivers23.dat");
  LoadSmear("../expdata/sivers24.dat");
  LoadSmear("../expdata/sivers25.dat");
  LoadSmear("../expdata/sivers26.dat");
  LoadSmear("../expdata/sivers27.dat");
  LoadSmear("../expdata/sivers28.dat");
  LoadSmear("../expdata/sivers29.dat");
  LoadSmear("../expdata/sivers30.dat");
  LoadSmear("../expdata/sivers31.dat");
  LoadSmear("../expdata/sivers32.dat");
  LoadSmear("../expdata/sivers33.dat");
  return Npoints;
}

int LoadEICSmear(){
  Npoints = 0;
  LoadSmear("../expdata/sivers00.dat");
  LoadSmear("../expdata/sivers01.dat");
  LoadSmear("../expdata/sivers02.dat");
  LoadSmear("../expdata/sivers03.dat");
  LoadSmear("../expdata/sivers04.dat");
  LoadSmear("../expdata/sivers05.dat");
  LoadSmear("../expdata/sivers06.dat");
  LoadSmear("../expdata/sivers07.dat");
  LoadSmear("../expdata/sivers08.dat");
  LoadSmear("../expdata/sivers09.dat");
  LoadSmear("../expdata/sivers10.dat");
  LoadSmear("../expdata/sivers11.dat");
  LoadSmear("../expdata/sivers12.dat");
  LoadSmear("../expdata/sivers13.dat");
  LoadSmear("../expdata/sivers20.dat");
  LoadSmear("../expdata/sivers21.dat");
  LoadSmear("../expdata/sivers22.dat");
  LoadSmear("../expdata/sivers23.dat");
  LoadSmear("../expdata/sivers24.dat");
  LoadSmear("../expdata/sivers25.dat");
  LoadSmear("../expdata/sivers26.dat");
  LoadSmear("../expdata/sivers27.dat");
  LoadSmear("../expdata/sivers28.dat");
  LoadSmear("../expdata/sivers29.dat");
  LoadSmear("../expdata/sivers30.dat");
  LoadSmear("../expdata/sivers31.dat");
  LoadSmear("../expdata/sivers32.dat");
  LoadSmear("../expdata/sivers33.dat");
  LoadSmear("simulate/prot_pip.dat");
  LoadSmear("simulate/prot_pim.dat");
  LoadSmear("simulate/prot_kp.dat");
  LoadSmear("simulate/prot_km.dat");
  LoadSmear("simulate/neut_pip.dat");
  LoadSmear("simulate/neut_pim.dat");
  LoadSmear("simulate/neut_kp.dat");
  LoadSmear("simulate/neut_km.dat");
  return Npoints;
}


int FitWorld(const int Nsample = 200){
  double init[9] = {0.172, 0.599, 5.07, -0.41, 2.36, 15.9, 0.059, -0.128, 1.32};
  double cent[9];
  FILE * fs = fopen("out-world.dat", "w");
  fprintf(fs, "idx\tNuv\tauv\tbuv\tNdv\tadv\tbdv\tNub\tNdb\tM1\n");
  for (int i = 0; i < Nsample; i++){
    LoadWorldSmear();
    Minimize(init, cent);
    fprintf(fs, "%d\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\n",
	    i, cent[0], cent[1], cent[2], cent[3], cent[4], cent[5], cent[6], cent[7], cent[8]);
  }
  fclose(fs);
  return 0;
}

int FitEIC(const int Nsample = 200){
  double init[9] = {0.172, 0.599, 5.07, -0.41, 2.36, 15.9, 0.059, -0.128, 1.32};
  double cent[9];
  FILE * fs = fopen("out-eic.dat", "w");
  fprintf(fs, "idx\tNuv\tauv\tbuv\tNdv\tadv\tbdv\tNub\tNdb\tM1\n");
  for (int i = 0; i < Nsample; i++){
    LoadEICSmear();
    Minimize(init, cent);
    fprintf(fs, "%d\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\n",
	    i, cent[0], cent[1], cent[2], cent[3], cent[4], cent[5], cent[6], cent[7], cent[8]);
  }
  fclose(fs);
  return 0;
}

int GetKT(const char * filename, const char * savefile, const double x, const double Q2){
  double KT[200];
  for (int i = 0; i < 200; i++){
    KT[i] = i * 0.01;
  }
  double u, d, ub, db, uv, dv;
  double du, dd, dub, ddb, duv, ddv;
  int idx;
  ifstream infile(filename);
  char tmp[300];
  infile.getline(tmp, 300);
  double kt;
  double par[200][9];
  for (int i = 0; i < 200; i++){
    infile >> idx >> par[i][0] >> par[i][1] >> par[i][2] >> par[i][3] >> par[i][4] >> par[i][5] >> par[i][6] >> par[i][7] >> par[i][8];
  }
  infile.close();
  FILE * fs = fopen(savefile, "w");
  fprintf(fs, "kt\tu\td\tub\tdb\tuv\tdv\tdu\tdd\tdub\tddb\tduv\tddv\n");
  for (int i = 0; i < 200; i++){
    kt = KT[i];
    u = 0; d = 0; ub = 0; db = 0;
    for (int j = 0; j <= idx; j++){
      SetParameters(par[j]);
      u += f1T(2, x, kt, Q2);
      d += f1T(1, x, kt, Q2);
      ub += f1T(-2, x, kt, Q2);
      db += f1T(-1, x, kt, Q2);
      uv += f1T(2, x, kt, Q2) - f1T(-2, x, kt, Q2);
      dv += f1T(1, x, kt, Q2) - f1T(-1, x, kt, Q2);
    }
    u = u / (idx + 1);
    d = d / (idx + 1);
    ub = ub / (idx + 1);
    db = db / (idx + 1);
    uv = uv / (idx + 1);
    dv = dv / (idx + 1);
    du = 0; dd = 0; dub = 0; ddb = 0; duv = 0; ddv = 0;
    for (int j = 0; j <= idx; j++){
      SetParameters(par[j]);
      du += pow(f1T(2, x, kt, Q2) - u, 2);
      dd += pow(f1T(1, x, kt, Q2) - d, 2);
      dub += pow(f1T(-2, x, kt, Q2) - ub, 2);
      ddb += pow(f1T(-1, x, kt, Q2) - db, 2);
      duv += pow(f1T(2, x, kt, Q2) - f1T(-2, x, kt, Q2) - uv, 2);
      ddv += pow(f1T(1, x, kt, Q2) - f1T(-1, x, kt, Q2) - dv, 2);
    }
    du = sqrt(du / idx);
    dd = sqrt(dd / idx);
    dub = sqrt(dub / idx);
    ddb = sqrt(ddb / idx);
    duv = sqrt(duv / idx);
    ddv = sqrt(ddv / idx);
    fprintf(fs, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n",
	    kt, u, d, ub, db, uv, dv, du, dd, dub, ddb, duv, ddv);
  }
  fclose(fs);
  return 0;
}

int main(const int argc, const char * argv[]){

  const int opt = atoi(argv[1]);

  double init[9] = {0.18,
		    1.0,
		    6.6,
		    -0.52,
		    1.9,
		    10.0,
		    -0.01,
		    -0.06,
		    0.9};
  double cent[9];
  SetParameters(init);
  if (opt == 0){//fit and simulate
    cout << LoadWorld() << endl;
    Minimize(init, cent);
    Simulation(cent);
  }

  if (opt == 1){//world fit samples
    FitWorld(5000);
  }

  if (opt == 2){//eic fit samples
    FitEIC(5000);
  }

  if (opt == 3){//generate kt files
    GetKT("out-world.dat", "f1t-world-x=0.005.dat", 0.005, 2.4);
    GetKT("out-world.dat", "f1t-world-x=0.05.dat", 0.05, 2.4);
    GetKT("out-world.dat", "f1t-world-x=0.1.dat", 0.1, 2.4);
    GetKT("out-world.dat", "f1t-world-x=0.2.dat", 0.2, 2.4);
    GetKT("out-eic.dat", "f1t-eic-x=0.005.dat", 0.005, 2.4);
    GetKT("out-eic.dat", "f1t-eic-x=0.05.dat", 0.05, 2.4);
    GetKT("out-eic.dat", "f1t-eic-x=0.1.dat", 0.1, 2.4);
    GetKT("out-eic.dat", "f1t-eic-x=0.2.dat", 0.2, 2.4);
  }
    

  return 0;
}
