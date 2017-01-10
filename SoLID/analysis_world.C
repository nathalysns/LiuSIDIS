#include <cstring>

#include "Lsidis.h"

using namespace std;

int main(int argc, char * argv[]){
  if (argc < 3){
    cout << "./analysis_world <readfile> <savefile>" << endl;
    return 0;
  }

  ifstream infile(argv[1]);
  string tmp;
  getline(infile, tmp, '\n');
  
  double i, Ebeam, x, y, z, Q2, pT, value, stat, systrel, systabs;
  string obs, target, hadron, experiment;
  double R;
 
  FILE * file = fopen(argv[2], "w");
  
  Lsidis mysidis;
  TLorentzVector l, P;
  int j = 0;
  //while (infile.good()){
    getline(infile, tmp, ',');
    i = atof(tmp);
    getline(infile, tmp, ',');
    Ebeam = atof(tmp);
    getline(infile, tmp, ',');
    x = atof(tmp);
    getline(infile, tmp, ',');
    y = atof(tmp);
    getline(infile, tmp, ',');
    z = atof(tmp);
    getline(infile, tmp, ',');
    Q2 = atof(tmp);
    getline(infile, tmp, ',');
    pT = atof(tmp);
    getline(infile, obs, ',');
    getline(infile, tmp, ',');
    value = atof(tmp);
    getline(infile, tmp, ',');
    stat = atof(tmp);
    getline(infile, tmp, ',');
    systrel = atof(tmp);
    getline(infile, tmp, ',');
    systabs = atof(tmp);
    getline(infile, target, ',');
    getline(infile, hadron, ',');
    getline(infile, experiment, '\n');
  //}
  printf("%.0f,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, obs.data(), value, stat, systrel, systabs, target.data(), hadron.data(), experiment.data());
  
  infile.close();

  return 0;
}
