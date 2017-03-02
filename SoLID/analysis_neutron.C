#include "SoLID_SIDIS_3He.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "./analysis_neutron <opt>" << endl;
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
    cout << "opt = 6: make kinematic covarage file" << endl;
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
  }

  if (opt == 1){
    if (argc == 2){
      GenerateBinInfoFile("NeutronResults/bin_info_N11p.dat", 11.0, "pi+");
      GenerateBinInfoFile("NeutronResults/bin_info_N8p.dat", 8.8, "pi+");
      GenerateBinInfoFile("NeutronResults/bin_info_N11m.dat", 11.0, "pi-");
      GenerateBinInfoFile("NeutronResults/bin_info_N8m.dat", 8.8, "pi-");
    }
    else if (argc == 3){
      Rfactor0 = atof(argv[2]);
      GenerateBinInfoFile("NeutronResults/bin_info_N11p_cut.dat", 11.0, "pi+");
      GenerateBinInfoFile("NeutronResults/bin_info_N8p_cut.dat", 8.8, "pi+");
      GenerateBinInfoFile("NeutronResults/bin_info_N11m_cut.dat", 11.0, "pi-");
      GenerateBinInfoFile("NeutronResults/bin_info_N8m_cut.dat", 8.8, "pi-");
    }
  }


  if (opt == 2){
    if (argc == 2){
      AnalyzeEstatUT3("NeutronResults/bin_info_N11p.dat", "NeutronResults/binsN11p.root", 11.0, "pi+");
      AnalyzeEstatUT3("NeutronResults/bin_info_N8p.dat", "NeutronResults/binsN8p.root", 8.8, "pi+");
      AnalyzeEstatUT3("NeutronResults/bin_info_N11m.dat", "NeutronResults/binsN11m.root", 11.0, "pi-");
      AnalyzeEstatUT3("NeutronResults/bin_info_N8m.dat", "NeutronResults/binsN8m.root", 8.8, "pi-");
    }
    else if (argc == 3){
      Rfactor0 = atof(argv[2]);
      AnalyzeEstatUT3("NeutronResults/bin_info_N11p_cut.dat", "NeutronResults/binsN11p_cut.root", 11.0, "pi+");
      AnalyzeEstatUT3("NeutronResults/bin_info_N8p_cut.dat", "NeutronResults/binsN8p_cut.root", 8.8, "pi+");
      AnalyzeEstatUT3("NeutronResults/bin_info_N11m_cut.dat", "NeutronResults/binsN11m_cut.root", 11.0, "pi-");
      AnalyzeEstatUT3("NeutronResults/bin_info_N8m_cut.dat", "NeutronResults/binsN8m_cut.root", 8.8, "pi-");
    }
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

  if (opt == 5){
    CreateFileSivers(argv[2], argv[3], argv[4]);
  }

  if (opt == 6){
    MakeKinematicCoveragePlots(11.0, "NeutronResults/coverage_11.root");
    MakeKinematicCoveragePlots(8.8, "NeutronResults/coverage_8.8.root");
  }	

  if (opt == 7){
    MakeRateDistributionPlots(11.0, "NeutronResults/ratedistri_11.root");
    MakeRateDistributionPlots(8.8, "NeutronResults/ratedistri_8.8.root");
  }

  if (opt == 8){
    Rfactor0 = 0.4;
    MakeRateDistributionPlots(11.0, "NeutronResults/ratedistri_11_R0.4.root");
    MakeRateDistributionPlots(8.8, "NeutronResults/ratedistri_8.8_R0.4.root");
  }

  return 0;
}
	

	



      

  
  

