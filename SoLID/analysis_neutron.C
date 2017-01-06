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
  
  gRandom->SetSeed(23);

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

  if (opt == 3){
    CheckCurrentCut(11.0, "pi+", 0.5, 0.5, 0.5, "NeutronResults/cutN11p.pdf");
    CheckCurrentCut(8.8, "pi+", 0.5, 0.5, 0.5, "NeutronResults/cutN8p.pdf");
    CheckCurrentCut(11.0, "pi-", 0.5, 0.5, 0.5, "NeutronResults/cutN11m.pdf");
    CheckCurrentCut(8.8, "pi-", 0.5, 0.5, 0.5, "NeutronResults/cutN8m.pdf");
  }

  if (opt == 4){
    double kT2 = atof(argv[2]);
    double MiT2 = atof(argv[3]);
    double MfT2 = atof(argv[4]);
    CheckCurrentCut(11.0, "pi+", kT2, MiT2, MfT2);
  }


  return 0;
}
	

	



      

  
  

