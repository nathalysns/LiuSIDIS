#include "SoLID_SIDIS_NH3.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "./analysis_proton <opt>" << endl;
    cout << "opt = 0: total rate" << endl;
    cout << "     ./analysis 0" << endl;
    cout << "opt = 1: binning data and create bin info file" << endl;
    cout << "     ./analysis 1" << endl;
    cout << "opt = 2: bin analysis including Estat" << endl;
    cout << "     ./analysis 2" << endl;
    cout << "opt = 3: output file for Sivers analysis" << endl;
    cout << "     ./analysis 3 <rootfile1> <rootfile2> <csvfile>" << endl;
    return 0;
  }

  int opt = atoi(argv[1]);
  LHAPDF::setVerbosity(0);
  
  gRandom->SetSeed(2);

  if (opt == 0){
    gRandom->SetSeed(0);
    cout << "Enhanced:" << endl;
    pimin = 0;
    GetTotalRate(11.0, "pi+");
    GetTotalRate(8.8, "pi+");
    GetTotalRate(11.0, "pi-");
    GetTotalRate(8.8, "pi-");
    cout << "Baseline:" << endl;
    pimin = 2.5;
    GetTotalRate(11.0, "pi+");
    GetTotalRate(8.8, "pi+");
    GetTotalRate(11.0, "pi-");
    GetTotalRate(8.8, "pi-");

  }

  if (opt == 1){
    pimin = 0;
    GenerateBinInfoFile("Projections/bin_enhanced_P11p.dat", 11.0, "pi+");
    GenerateBinInfoFile("Projections/bin_enhanced_P8p.dat", 8.8, "pi+");
    GenerateBinInfoFile("Projections/bin_enhanced_P11m.dat", 11.0, "pi-");
    GenerateBinInfoFile("Projections/bin_enhanced_P8m.dat", 8.8, "pi-");
    pimin = 2.5;
    GenerateBinInfoFile("Projections/bin_base_P11p.dat", 11.0, "pi+");
    GenerateBinInfoFile("Projections/bin_base_P8p.dat", 8.8, "pi+");
    GenerateBinInfoFile("Projections/bin_base_P11m.dat", 11.0, "pi-");
    GenerateBinInfoFile("Projections/bin_base_P8m.dat", 8.8, "pi-");
  }

  if (opt == -1){
    Rfactor0 = 0.4;
    pimin = 0;
    GenerateBinInfoFile("Projections/bin_enhanced_P11p_cut.dat", 11.0, "pi+");
    GenerateBinInfoFile("Projections/bin_enhanced_P8p_cut.dat", 8.8, "pi+");
    GenerateBinInfoFile("Projections/bin_enhanced_P11m_cut.dat", 11.0, "pi-");
    GenerateBinInfoFile("Projections/bin_enhanced_P8m_cut.dat", 8.8, "pi-");
    pimin = 2.5;
    GenerateBinInfoFile("Projections/bin_base_P11p_cut.dat", 11.0, "pi+");
    GenerateBinInfoFile("Projections/bin_base_P8p_cut.dat", 8.8, "pi+");
    GenerateBinInfoFile("Projections/bin_base_P11m_cut.dat", 11.0, "pi-");
    GenerateBinInfoFile("Projections/bin_base_P8m_cut.dat", 8.8, "pi-");
  }

  if (opt == 2){
    pimin = 0;
    AnalyzeEstatUT3("Projections/bin_enhanced_P11p.dat", "Projections/enhancedP11p.root", 11.0, "pi+");
    AnalyzeEstatUT3("Projections/bin_enhanced_P8p.dat", "Projections/enhancedP8p.root", 8.8, "pi+");
    AnalyzeEstatUT3("Projections/bin_enhanced_P11m.dat", "Projections/enhancedP11m.root", 11.0, "pi-");
    AnalyzeEstatUT3("Projections/bin_enhanced_P8m.dat", "Projections/enhancedP8m.root", 8.8, "pi-");
    pimin = 2.5;
    AnalyzeEstatUT3("Projections/bin_base_P11p.dat", "Projections/baseP11p.root", 11.0, "pi+");
    AnalyzeEstatUT3("Projections/bin_base_P8p.dat", "Projections/baseP8p.root", 8.8, "pi+");
    AnalyzeEstatUT3("Projections/bin_base_P11m.dat", "Projections/baseP11m.root", 11.0, "pi-");
    AnalyzeEstatUT3("Projections/bin_base_P8m.dat", "Projections/baseP8m.root", 8.8, "pi-");
  }

  if (opt == 3){
    CreateFileSivers("Projections/enhancedP11p.root", "Projections/enhancedP8p.root", "enhancedPpip.csv");
    CreateFileSivers("Projections/enhancedP11m.root", "Projections/enhancedP8m.root", "enhancedPpim.csv");
    CreateFileSivers("Projections/baseP11p.root", "Projections/baseP8p.root", "basePpip.csv");
    CreateFileSivers("Projections/baseP11m.root", "Projections/baseP8m.root", "basePpim.csv");
  }


  return 0;
}
	
