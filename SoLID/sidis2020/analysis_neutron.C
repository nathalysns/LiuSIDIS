#include "SoLID_SIDIS_3He.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "./analysis_neutron <opt>" << endl;
    cout << "opt = 0: total rate" << endl;
    cout << "     ./analysis 0" << endl;
    cout << "opt = 1: binning data and create bin info file" << endl;
    cout << "     ./analysis 1" << endl;
    cout << "opt = 2: bin analysis including Estat" << endl;
    cout << "     ./analysis 2" << endl;
    cout << "opt = 3: output file for Sivers analysis" << endl;
    cout << "     ./analysis 5 <rootfile1> <rootfile2> <csvfile>" << endl;
    return 0;
  }

  int opt = atoi(argv[1]);
  
  gRandom->SetSeed(2);
  LHAPDF::setVerbosity(0);

  if (opt == 0){
    gRandom->SetSeed(0);
    cout << "Enhanced:" << endl;
    GetTotalRate(11.0, "pi+");
    GetTotalRate(8.8, "pi+");
    GetTotalRate(11.0, "pi-");
    GetTotalRate(8.8, "pi-");
    pimin = 2.5;
    cout << "Baseline:" << endl;
    GetTotalRate(11.0, "pi+");
    GetTotalRate(8.8, "pi+");
    GetTotalRate(11.0, "pi-");
    GetTotalRate(8.8, "pi-");
  }

  if (opt == 1){
    pimin = 0;
    GenerateBinInfoFile("Projections/bin_enhanced_N11p.dat", 11.0, "pi+");
    GenerateBinInfoFile("Projections/bin_enhanced_N8p.dat", 8.8, "pi+");
    GenerateBinInfoFile("Projections/bin_enhanced_N11m.dat", 11.0, "pi-");
    GenerateBinInfoFile("Projections/bin_enhanced_N8m.dat", 8.8, "pi-");
    pimin = 2.5;
    GenerateBinInfoFile("Projections/bin_base_N11p.dat", 11.0, "pi+");
    GenerateBinInfoFile("Projections/bin_base_N8p.dat", 8.8, "pi+");
    GenerateBinInfoFile("Projections/bin_base_N11m.dat", 11.0, "pi-");
    GenerateBinInfoFile("Projections/bin_base_N8m.dat", 8.8, "pi-");
  }


  if (opt == 2){
    pimin = 0;
    AnalyzeEstatUT3("Projections/bin_enhanced_N11p.dat", "Projections/enhancedN11p.root", 11.0, "pi+");
    AnalyzeEstatUT3("Projections/bin_enhanced_N8p.dat", "Projections/enhancedN8p.root", 8.8, "pi+");
    AnalyzeEstatUT3("Projections/bin_enhanced_N11m.dat", "Projections/enhancedN11m.root", 11.0, "pi-");
    AnalyzeEstatUT3("Projections/bin_enhanced_N8m.dat", "Projections/enhancedN8m.root", 8.8, "pi-");
    pimin = 2.5;
    AnalyzeEstatUT3("Projections/bin_base_N11p.dat", "Projections/baseN11p.root", 11.0, "pi+");
    AnalyzeEstatUT3("Projections/bin_base_N8p.dat", "Projections/baseN8p.root", 8.8, "pi+");
    AnalyzeEstatUT3("Projections/bin_base_N11m.dat", "Projections/baseN11m.root", 11.0, "pi-");
    AnalyzeEstatUT3("Projections/bin_base_N8m.dat", "Projections/baseN8m.root", 8.8, "pi-");
  }

  if (opt == 3){
    CreateFileSivers("Projections/enhancedN11p.root", "Projections/enhancedN8p.root", "enhancedNpip.csv");
    CreateFileSivers("Projections/enhancedN11m.root", "Projections/enhancedN8m.root", "enhancedNpim.csv");
    CreateFileSivers("Projections/baseN11p.root", "Projections/baseN8p.root", "baseNpip.csv");
    CreateFileSivers("Projections/baseN11m.root", "Projections/baseN8m.root", "baseNpim.csv");
  }

  return 0;
}
	

	



      

  
  

