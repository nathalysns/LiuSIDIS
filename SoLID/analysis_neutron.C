#include "SoLID_SIDIS_3He.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "./analysis_neutron <opt>" << endl;
    cout << "opt = 0: total rate" << endl;
    cout << "opt = 1: binning data and create bin info file" << endl;
    cout << "opt = 2: bin analysis including Estat" << endl;   
    return 0;
  }

  int opt = atoi(argv[1]);
  
  gRandom->SetSeed(2);

  if (opt == 0){
    GetTotalRate(11.0, "pi+");
    GetTotalRate(8.8, "pi+");
    GetTotalRate(11.0, "pi-");
    GetTotalRate(8.8, "pi-");
  }

  if (opt == 1){
    GenerateBinInfoFile("NeutronResults/bin_info_N11p.dat", 11.0, "pi+");
    GenerateBinInfoFile("NeutronResults/bin_info_N8p.dat", 8.8, "pi+");
    GenerateBinInfoFile("NeutronResults/bin_info_N11m.dat", 11.0, "pi-");
    GenerateBinInfoFile("NeutronResults/bin_info_N8m.dat", 8.8, "pi-");
  }

  if (opt == 2){
    AnalyzeEstatUT3("NeutronResults/bin_info_N11p.dat", "NeutronResults/binsN11p.root", 11.0, "pi+");
    AnalyzeEstatUT3("NeutronResults/bin_info_N8p.dat", "NeutronResults/binsN8p.root", 8.8, "pi+");
    AnalyzeEstatUT3("NeutronResults/bin_info_N11m.dat", "NeutronResults/binsN11m.root", 11.0, "pi-");
    AnalyzeEstatUT3("NeutronResults/bin_info_N8m.dat", "NeutronResults/binsN8m.root", 8.8, "pi-");
  }

  return 0;
}
	

	



      

  
  

