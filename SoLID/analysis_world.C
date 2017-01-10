#include <cstring>

#include "Lsidis.h"

using namespace std;

int main(int argc, char * argv[]){
  if (argc < 4){
    cout << "./analysis_world <R0> <readfile> <savefile> " << endl;
    return 0;
  }

  double R0 = atof(argv[1]);
	
  ifstream infile(argv[2]);
  string tmp;
  getline(infile, tmp, '\n');
  
  double i, Ebeam, x, y, z, Q2, pT, value, stat, systrel, systabs;
  string obs, target, hadron, experiment;
   
  FILE * file = fopen(argv[3], "w");
  
  Lsidis mysidis;
  mysidis.SetHadron("pi+");
  mysidis.SetPDFset("CT14lo");
  mysidis.SetFFset("DSSFFlo");
  TLorentzVector l, P;
  int j = 0;
  while (infile.good()){
    getline(infile, tmp, ',');
    i = atof(tmp.data());
    getline(infile, tmp, ',');
    Ebeam = atof(tmp.data());
    getline(infile, tmp, ',');
    x = atof(tmp.data());
    getline(infile, tmp, ',');
    y = atof(tmp.data());
    getline(infile, tmp, ',');
    z = atof(tmp.data());
    getline(infile, tmp, ',');
    Q2 = atof(tmp.data());
    getline(infile, tmp, ',');
    pT = atof(tmp.data());
    getline(infile, obs, ',');
    getline(infile, tmp, ',');
    value = atof(tmp.data());
    getline(infile, tmp, ',');
    stat = atof(tmp.data());
    getline(infile, tmp, ',');
    systrel = atof(tmp.data());
    getline(infile, tmp, ',');
    systabs = atof(tmp.data());
    getline(infile, target, ',');
    getline(infile, hadron, ',');
    getline(infile, experiment, '\n');
    if (!infile.good()) break;
	
    if (target == "proton") mysidis.SetNucleus(1, 0);
    else if (target == "neutron") mysidis.SetNucleus(0, 1);
    else if (target == "deuteron") mysidis.SetNucleus(1, 1);
	
    if (hadron == "pi+" || hadron == "h+") mysidis.SetHadron("pi+");
    else if (hadron == "pi-" || hadron == "h-") mysidis.SetHadron("pi-");
    else if (hadron == "pi0") mysidis.SetHadron("pi0");
    else if (hadron == "k+") mysidis.SetHadron("K+");
    else if (hadron == "k-") mysidis.SetHadron("K-");
    else if (hadron == "k0") mysidis.SetHadron("K0");
	
    l.SetXYZT(0, 0, Ebeam, Ebeam);
    P.SetXYZT(0, 0, 0, 0.938272);
    mysidis.SetInitialState(l, P);
    
    mysidis.SetVariables(x, y, z, pT, 0.0, 0.0);
    mysidis.CalculateFinalState();

    if (mysidis.GetVariable("W") < 2.3) continue;
    if (mysidis.GetVariable("Wp") < 1.6) continue;
	  
    mysidis.CalculateRfactor();
    if (mysidis.GetVariable("Rfactor") > R0) continue;
  
    fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      j++, Ebeam, x, y, z, Q2, pT, obs.data(), value, stat, systrel, systabs, target.data(), hadron.data(), experiment.data());
  }
  
  infile.close();
  fclose(file);

  return 0;
}
