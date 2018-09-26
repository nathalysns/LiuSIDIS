#include "SoLID_SIDIS_NH3.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "./analysis_proton <opt>" << endl;
    cout << "opt = 0: total rate" << endl;
    cout << "     ./analysis 0" << endl;
    cout << "     ./analysis 0 <R0>" << endl;
    cout << "opt = 1: binning data and create bin info file" << endl;
    cout << "     ./analysis 1" << endl;
    cout << "     ./analysis 1 <R0>" << endl;
    cout << "opt = 2: bin analysis including Estat" << endl;
    cout << "     ./analysis 2" << endl;
    cout << "     ./analysis 2 <R0>" << endl;
    cout << "opt = 3: compare current cut effect" << endl;
    cout << "opt = 4: see parameter sensitivity of the current cut" << endl;
    cout << "     ./analysis 4 <kT2> <MiT2> <MfT2>" << endl;
    cout << "opt = 5: output file for Sivers analysis" << endl;
    cout << "     ./analysis 5 <rootfile1> <rootfile2> <csvfile>" << endl;
    cout << "opt = 6: make kinematic coverage file" << endl;
    cout << "     ./analysis 6 " << endl;
    cout << "opt = 7: make rate distribution file" << endl;
    cout << "     ./analysis 7 " << endl;
    cout << "opt = 8: make rate distribution file with R < 0.4" << endl;
    cout << "     ./analysis 8 " << endl;
    return 0;
  }

  int opt = atoi(argv[1]);
  
  gRandom->SetSeed(2);

  if (opt == 0){
    gRandom->SetSeed(0);
    if (argc > 2) Rfactor0 = atof(argv[2]);
    GetTotalRate(11.0, "pi+");
    GetTotalRate(8.8, "pi+");
    GetTotalRate(11.0, "pi-");
    GetTotalRate(8.8, "pi-");
    PKmax = 7.5;
    GetTotalRate(11.0, "K+");
    GetTotalRate(11.0, "K-");
    PKmax = 5.0;
    GetTotalRate(11.0, "K+");
    GetTotalRate(11.0, "K-");
    PKmax = 2.5;
    GetTotalRate(11.0, "K+");
    GetTotalRate(11.0, "K-");
  }

  if (opt == 1){
    if (argc == 2){
      GenerateBinInfoFile("Projections/bin_info_P11p.dat", 11.0, "pi+");
      GenerateBinInfoFile("Projections/bin_info_P8p.dat", 8.8, "pi+");
      GenerateBinInfoFile("Projections/bin_info_P11m.dat", 11.0, "pi-");
      GenerateBinInfoFile("Projections/bin_info_P8m.dat", 8.8, "pi-");
    }
    else if (argc == 3){
      Rfactor0 = atof(argv[2]);
      GenerateBinInfoFile("Projections/bin_info_P11p_cut.dat", 11.0, "pi+");
      GenerateBinInfoFile("Projections/bin_info_P8p_cut.dat", 8.8, "pi+");
      GenerateBinInfoFile("Projections/bin_info_P11m_cut.dat", 11.0, "pi-");
      GenerateBinInfoFile("Projections/bin_info_P8m_cut.dat", 8.8, "pi-");
    }
  }

  if (opt == 2){
    if (argc == 2){
      AnalyzeEstatUT3("Projections/bin_info_P11p.dat", "Projections/binsP11p.root", 11.0, "pi+");
      AnalyzeEstatUT3("Projections/bin_info_P8p.dat", "Projections/binsP8p.root", 8.8, "pi+");
      AnalyzeEstatUT3("Projections/bin_info_P11m.dat", "Projections/binsP11m.root", 11.0, "pi-");
      AnalyzeEstatUT3("Projections/bin_info_P8m.dat", "Projections/binsP8m.root", 8.8, "pi-");
    }
    else if (argc == 3){
      Rfactor0 = atof(argv[2]);
      AnalyzeEstatUT3("Projections/bin_info_P11p_cut.dat", "Projections/binsP11p_cut.root", 11.0, "pi+");
      AnalyzeEstatUT3("Projections/bin_info_P8p_cut.dat", "Projections/binsP8p_cut.root", 8.8, "pi+");
      AnalyzeEstatUT3("Projections/bin_info_P11m_cut.dat", "Projections/binsP11m_cut.root", 11.0, "pi-");
      AnalyzeEstatUT3("Projections/bin_info_P8m_cut.dat", "Projections/binsP8m_cut.root", 8.8, "pi-");
    }
  }

  if (opt == 5){
    CreateFileSivers("Projections/binsP11p.root", "Projections/binsP8p.root", "base03.csv");
    CreateFileSivers("Projections/binsP11m.root", "Projections/binsP8m.root", "base04.csv");
  }


  return 0;
}
	

	



      

  
  

