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
 
  FILE * file = fopen(argv[2], "w");
  
  Lsidis mysidis;
  TLorentzVector l, P;
  int j = 0;
  //while (infile >> i >> Ebeam >> x >> y >> z >> Q2 >> pT >> obs >> value >> stat >> systrel >> systabs >> target >> hadron >> experiment){
  //}
  
  printf("%s", tmp.data());
  
  infile.close();

  return 0;
}
