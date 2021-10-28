/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* StructureFunctions.cc                                                                                                   */
/* modified by: Alexander Epping: a_eppi01@uni-muenster.de                                                                 */
/* GitHub: https://github.com/alexanderepping/apfelxx_photon_mod                                                           */
/* 25 Oct 2021                                                                                                             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* original file: https://github.com/vbertone/APFEL_Examples                                                               */
/* original author: Valerio Bertone: valerio.bertone@cern.ch                                                               */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Program to calculate the Structure Functions from a given LHAPDF set.                                                   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program takes an evolved LHAPDF set (produced by EvolutionStructureFunctions.cc) as input and calculates           */
/* the Structure Functions from that. Their values are then output to a file and the terminal.                             */
/*                                                                                                                         */
/* The output file is saved in plottoingPython/ and can be used by plottingStructureFunctions.                             */
/*                                                                                                                         */
/* This program and the EvolutionStructureFunctions.cc program can be run by using the run_StructureFunctions.sh file      */
/* in the bashFiles directory.                                                                                             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



///////////////////////////////////////
// imports
///////////////////////////////////////

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include "apfel/apfelxx.h" 

// used to get the used perturbative order of the pointlike contributions
#include "apfel/pointlikecontributions.h"

// used to write to files
#include<fstream>



///////////////////////////////////////
// definitions; can be changed by user
///////////////////////////////////////

// Name of the output file
const std::string OutputFile = "/home/alexander/Documents/apfelxx_photon_mod/plottingPython/data_StructureFunctions.txt";

// names of the used LHAPDF sets
const std::string NameLHAPDFSet = "GRVCustomSetLO";
const std::string NameLHAPDFSetEvolved = "Evolved" + NameLHAPDFSet;

// array of final scale values for which data should be output (mu^2)
//                 |  0,1,2:ALEPH  |  3,4,5:AMY
const double arr_mu2[] = {9.9, 20.7, 284, 6.8, 73, 390};   

// Vector of test values of x (xlha used for console output, xlha2 used for file output and plotting)
std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
//std::vector<double> xlha2{0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
std::vector<double> xlha2{0.0001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110, 0.115, 0.120, 0.125, 0.130, 0.135, 0.140, 0.145, 0.150, 0.155, 0.160, 0.165, 0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.200, 0.205, 0.210, 0.215, 0.220, 0.225, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.300, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355, 0.360, 0.365, 0.370, 0.375, 0.380, 0.385, 0.390, 0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425, 0.430, 0.435, 0.440, 0.445, 0.450, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505, 0.510, 0.515, 0.520, 0.525, 0.530, 0.535, 0.540, 0.545, 0.550, 0.555, 0.560, 0.565, 0.570, 0.575, 0.580, 0.585, 0.590, 0.595, 0.600, 0.605, 0.610, 0.615, 0.620, 0.625, 0.630, 0.635, 0.640, 0.645, 0.650, 0.655, 0.660, 0.665, 0.670, 0.675, 0.680, 0.685, 0.690, 0.695, 0.700, 0.705, 0.710, 0.715, 0.720, 0.725, 0.730, 0.735, 0.740, 0.745, 0.750, 0.755, 0.760, 0.765, 0.770, 0.775, 0.780, 0.785, 0.790, 0.795, 0.800, 0.805, 0.810, 0.815, 0.820, 0.825, 0.830, 0.835, 0.840, 0.845, 0.850, 0.855, 0.860, 0.865, 0.870, 0.875, 0.880, 0.885, 0.890, 0.895, 0.900};//, 0.905, 0.910, 0.915, 0.920, 0.925, 0.930, 0.935, 0.940, 0.945, 0.950, 0.955, 0.960, 0.965, 0.970, 0.975, 0.980, 0.985, 0.990, 0.995}; 


///////////////////////////////////////
// main program
///////////////////////////////////////

int main()
{
  ///////////////////////////////////////
  // opening and preparing output file
  ///////////////////////////////////////

  // opening output file
  std::ofstream file;
  file.open(OutputFile);

  // print basic information () on the following data
  file << "# mu2 values:" << std::endl;

  for (double mu2 : arr_mu2) 
  {
    if (mu2 == arr_mu2[0])
      file << std::to_string(mu2);
    else
      file << ", " << std::to_string(mu2);
  }
  
  file << "\n# num_x_vals\n" << xlha2.size() << "\n# x, APFEL++, LHAPDF, ratio" << std::endl;



  /////////////////////////////////////
  // definitions
  /////////////////////////////////////

  // flavor number with corresponding flavor
  std::map<int, std::string> mapFlavors {{-6, "tbar"}, {-5, "bbar"}, {-4, "cbar"}, {-3, "sbar"}, {-2, "ubar"}, {-1, "dbar"},
                                         {0, "gluon"}, {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}};

  
  // Open LHAPDF sets
  LHAPDF::PDF* dist = LHAPDF::mkPDF(NameLHAPDFSet);
  LHAPDF::PDF* distEvolved = LHAPDF::mkPDF(NameLHAPDFSetEvolved);
  

  // Retrieve evolution parameters from the LHAPDF set
  const int    pto          = dist->orderQCD();
  const double Qref         = 91.1876; // mass of the z-boson
  const double asref        = dist->alphasQ(Qref);
  const double mc           = dist->quarkThreshold(4);
  const double mb           = dist->quarkThreshold(5);
  const double mt           = dist->quarkThreshold(6);
  const double Qin          = dist->qMin();
  

 
  /////////////////////////////////////
  // Initialise APFEL++
  /////////////////////////////////////

  // Define grid in terms of subgrids. Each subgrid is constructed
  // through apfel::SubGrid{nx, xmin, id} and has "nx" intervals, is
  // logarithmically distributed between "xmin" and 1, and has "id"
  // interpolation degree.
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 6e-1, 3}, apfel::SubGrid{50, 8e-1, 5}}};

  // Vector of thresholds. Not including the top thesholds
  // automatically means that the code will work in a nfmax = 5
  // scheme.
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb};

  // Use alpha_s from LHAPDF
  const auto as = [&] (double const& mu) -> double{ return dist->alphasQ(mu); };

  // Effective EW charges for a space-like process (i.e. DIS)
  const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

  // Define input PDFs as a lambda function using the LHAPDF
  // set. Notice that APFEL++ takes PDFs in the QCD evolution basis
  // (g, Sigma, V, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35). The
  // function "PhysToQCDEv" rotates PDFs from the physical basis
  // (tbar, bbar, ..., g, ..., b, t) into the physical basis.
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); };
  const auto PDFsEvolved = [&] (double const& x, double const& Q) -> std::map<int, double>{return apfel::PhysToQCDEv(distEvolved->xfxQ(x, Q));};

  // Initialise F2 NC structure function in the ZM-VFNS. Functions
  // exist for FL and xF3 (btw, APFEL produces results for xF3, not
  // F3). Also CC structure functions are implemented in APFEL++. Also
  // the massive scheme and its massless limit are supported. At the
  // moment CC massive structure functions aren't there but can be
  // easily implemented if needed.
  const auto F2 = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, pto, as, fBq);
  const auto F2Evolved = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolved, pto, as, fBq);

  // Tabulate the ftructure function over 50 intervals in the range
  // [1: 1000] GeV using interpolation degree 3. Take into account
  // that at thresholds straucture functions may not be continuos
  // feeding the contructor with the thresholds.
  const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2totalEvolved {[&] (double const& Q) -> apfel::Distribution{ return F2Evolved.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};



  /////////////////////////////////////
  // output results
  /////////////////////////////////////

  for (int i_mu2=0; i_mu2<sizeof(arr_mu2)/sizeof(arr_mu2[0]); i_mu2++)
  {
    double mu = pow(arr_mu2[i_mu2], 0.5);

    std::cout << std::scientific;

    // terminal output:
    // Print down PDF at "mu"
    std::cout << "____________________________________________________________" << std::endl << std::endl;
    std::cout << "Used Mu                   : " << mu <<  " GeV" << std::endl << std::endl;
    std::cout << "     x            APFEL++/           LHAPDF           ratio" << std::endl;
    std::cout << "                  Evolved                                  " << std::endl;
    for (double x : xlha)
      std::cout << x << "\t" << F2totalEvolved.EvaluatexQ(x, mu) << "\t" << F2total.EvaluatexQ(x, mu) << "\t" << F2totalEvolved.EvaluatexQ(x, mu)/F2total.EvaluatexQ(x, mu) << std::endl;


    // file output:
    // Print down PDF at "mu"
    for (double x : xlha2)
      file << x << ", " << F2totalEvolved.EvaluatexQ(x, mu) << ", " << F2total.EvaluatexQ(x, mu) << ", " << F2totalEvolved.EvaluatexQ(x, mu)/F2total.EvaluatexQ(x, mu) << std::endl;
  }

  // terminal output:
  // print perturbative order and name of LHAPDF-Set
  std::cout << "____________________________________________________________" << std::endl;
  std::cout << "____________________________________________________________\n" << std::endl;
  std::cout << std::defaultfloat;
  std::cout << "Used Perturbative Order   : " << std::to_string(pto) << std::endl;
  std::cout << "Used Perturbative Order PL: " << std::to_string(apfel::ptoPL) << "\n" << std::endl;
  std::cout << "Used ZBoson Mass          : " << Qref << " GeV" << std::endl;
  std::cout << "Used alphas @ ZBoson Mass : " << asref << "\n" << std::endl;
  std::cout << "Used Charm Quark Mass     : " << mc << " GeV" << std::endl;
  std::cout << "Used Charm Quark Mass     : " << mc << " GeV" << std::endl;
  std::cout << "Used Bottom Quark Mass    : " << mb << " GeV" << std::endl;
  std::cout << "Used Top Quark Mass       : " << mt << " GeV\n" << std::endl;
  std::cout << "Used LHAPDF Set           : " << NameLHAPDFSet << std::endl;
  
  // close output file
  file.close();

  return 0;
}
