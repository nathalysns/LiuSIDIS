#include "SoLID_SIDIS_3He.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "./analysis_neutron <opt>" << endl;
    cout << "opt = 0: binning data and create bin info file" << endl;
    return 0;
  }

  int opt = atoi(argv[1]);

  if (opt == 0){
    gRandom->SetSeed(2);
    GenerateBinInfoFile("NeutronResults/bin_info_N11p.dat", 11.0, "pi+");
    GenerateBinInfoFile("NeutronResults/bin_info_N8p.dat", 8.8, "pi+");
    GenerateBinInfoFile("NeutronResults/bin_info_N11m.dat", 11.0, "pi-");
    GenerateBinInfoFile("NeutronResults/bin_info_N8m.dat", 8.8, "pi-");
  }



  return 0;
}
	

	



      

  
  

