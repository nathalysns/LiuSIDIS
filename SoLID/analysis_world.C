#include <cstring>

#include "Lsidis.h"

using namespace std;

int main(int argc, char * argv[]){
  if (argc < 3){
    cout << "./analysis_world <readfile> <savefile>" << endl;
    return 0;
  }

  string header;
  
  ifstream infile(argv[1]);
  getline(infile, header);
  cout << header << endl;
  /*
  
  getline(infile, header, '\r');
  int Ncol = 1;
  for (unsigned int i = 0; i < header.length(); i++){
    if (header[i] == ',') Ncol++;
  }
  infile.clear();
  infile.seekg(0, ios::beg);

  string column[Ncol];
  for (int i = 0; i < Ncol - 1; i++){
    getline(infile, column[i], ',');
  }
  getline(infile, column[Ncol-1], '\r');
  
  cout << header << endl;
  cout << Ncol << endl;
  for (int i = 0; i < Ncol; i++)
    cout << column[i] << endl;
*/
  infile.close();

  return 0;
}
