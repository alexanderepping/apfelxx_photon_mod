/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* EvolutionFlavors.cc                                                                                                     */
/* modified by: Alexander Epping: a_eppi01@uni-muenster.de                                                                 */
/* GitHub: https://github.com/alexanderepping/apfelxx_photon_mod                                                           */
/* 22 Oct 2021                                                                                                             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* original file: https://github.com/vbertone/APFEL_Examples                                                               */
/* original author: Valerio Bertone: valerio.bertone@cern.ch                                                               */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Program to evolve the PDFs with given initial conditions and print out the results for multiple energies and particles. */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program takes a LHAPDF set as input and evolves the PDFs using the DGLAP evolution, starting at an initial         */
/* energy. The evolved PDFs are then output to a file and the terminal at multiple final energies (mu) and for different   */
/* quark flavors as well as for the gluon.                                                                                 */
/*                                                                                                                         */
/* The output file can be read by plottingEvolutionFlavors.py in the plottingPython directory.                             */
/*                                                                                                                         */
/* This program and the plottingEvolutionFlavors program can be run by using the run_EvolutionFlavors.sh file              */
/* in the bashFiles directory.                                                                                             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



///////////////////////////////////////
// imports
///////////////////////////////////////

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include "apfel/apfelxx.h" 

// used to write to files
#include<fstream>



///////////////////////////////////////
// definitions; can be changed by user
///////////////////////////////////////

// Name of the output file
const std::string OutputFile = "/home/alexander/Documents/apfelxx_photon_mod/plottingPython/data_EvolutionFlavors.txt";

// name of the used LHAPDF set
const std::string NameLHAPDFSet = "GRVCustomSet";

// array of final scale values for which data should be output
double arr_mu[] = {1.295, 4, 10, 20}; 

// Vector of test values of x (xlha used for console output, xlha2 used for file output and plotting)
std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
std::vector<double> xlha2{0.0001, 0.0002, 0.00030000000000000003, 0.0004, 0.0005, 0.0006000000000000001, 0.0007, 0.0008, 0.0009000000000000001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009000000000000001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6000000000000001, 0.7000000000000001, 0.8, 0.9};



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
  file << "# mu values:" << std::endl;

  for (double mu : arr_mu) 
  {
    if (mu == arr_mu[0])
      file << std::to_string(mu);
    else
      file << ", " << std::to_string(mu);
  }
  
  file << "\n# num_x_vals\n" << xlha2.size() << "\n# x, apfel++, lhapdf, ratio" << std::endl;



  /////////////////////////////////////
  // definitions
  /////////////////////////////////////

  // flavor number with corresponding flavor
  std::map<int, std::string> mapFlavors {{-6, "tbar"}, {-5, "bbar"}, {-4, "cbar"}, {-3, "sbar"}, {-2, "ubar"}, {-1, "dbar"},
                                         {0, "gluon"}, {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}};

  
  // Open LHAPDF set
  LHAPDF::PDF* dist = LHAPDF::mkPDF(NameLHAPDFSet);
  

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

  // Vector of threholds. Not including the top thesholds
  // automatically means that the code will work in a nfmax = 5
  // scheme.
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb};

  // Vector of masses (set equal to thresholds but can be different)
  const std::vector<double> Masses = Thresholds;

  // Construct alpha_s object
  apfel::AlphaQCD a{asref, Qref, Masses, Thresholds, pto};
  //AlphaQCD a{asref, Qref, Thresholds, pto}; // This constructor can also be used if thresholds and masses are equal

  // Tabulate alpha_s over an interval between 0.9 and 1001 GeV on a
  // grid with 100 intervals distributed as log(log(Q/Lambda)) and set
  // interpolation degree to 3.
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Define input PDFs as a lambda function using the LHAPDF
  // set. Notice that APFEL++ takes PDFs in the QCD evolution basis
  // (g, Sigma, V, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35). The
  // function "PhysToQCDEv" rotates PDFs from the physical basis
  // (tbar, bbar, ..., g, ..., b, t) into the physical basis. One can
  // use this function to define and set of functions to be evolved.
  const auto InPDFs = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); };

  // Initialise evolution setting "Qin" as an initial scale. This
  // means that the function "InPDFs" will be called at "Qin". It uses
  // alpha_s defined above but any other function can be used.
  auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs, Qin, pto, as);
  //auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), InPDFs, Qin, pto, as); // This constructor can also be used if thresholds and masses are equal

  // Finally tabulate PDFs. The constructor is the same as that of
  // "as".
  //const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};



  /////////////////////////////////////
  // output results
  /////////////////////////////////////

  for (int i_mu=0; i_mu<sizeof(arr_mu)/sizeof(arr_mu[0]); i_mu++)
  {
    double mu = arr_mu[i_mu];

    std::cout << std::scientific;

    // terminal output:
    // print alphas from APFEL++, and LHAPDF only once per mu
    std::cout << "____________________________________________________________" << std::endl;
    std::cout << "____________________________________________________________" << std::endl;
    std::cout << "\nmu = " << mu << " GeV\n" << std::endl;
    std::cout << "APFEL++: AlphaQCD(Q) = " << Alphas.Evaluate(mu) << std::endl;
    std::cout << "LHAPDF:  AlphaQCD(Q) = " << dist->alphasQ(mu) << std::endl;

    for (int flavor=0; flavor<=6; flavor++)
    {

      // Call evolved PDFs at "mu". Notice that APFEL++ returns PDFs in
      // the QCD evolution basis, therefore one needs to rotate them back
      // to the physical basis which is done through the function
      // "QCDEvToPhys".
      const std::map<int, apfel::Distribution> tpdfs = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects());


      // terminal output:
      // Print down PDF at "mu"
      std::cout << "____________________________________________________________" << std::endl << std::endl;
      std::cout << "Shown particle:           : " << mapFlavors.at(flavor) << std::endl;
      std::cout << "Used Mu                   : " << mu <<  " GeV" << std::endl << std::endl;
      std::cout << "     x             APFEL++           LHAPDF           ratio" << std::endl;
      for (double x : xlha)
        std::cout << x << "\t" << tpdfs.at(flavor).Evaluate(x) << "\t" << dist->xfxQ(flavor, x, mu) << "\t" << tpdfs.at(flavor).Evaluate(x)/dist->xfxQ(flavor, x, mu) << std::endl;


      // file output:
      // Print down PDF at "mu"
      for (double x : xlha2)
        file << x << ", " << tpdfs.at(flavor).Evaluate(x) << ", " << dist->xfxQ(flavor, x, mu) << ", " << tpdfs.at(flavor).Evaluate(x)/dist->xfxQ(flavor, x, mu) << std::endl;
    }
  }

  // terminal output:
  // print perturbative order and name of LHAPDF-Set
  std::cout << "____________________________________________________________" << std::endl;
  std::cout << "____________________________________________________________\n" << std::endl;
  std::cout << std::defaultfloat;
  std::cout << "Used Perturbative Order   : " << std::to_string(pto) << "\n" << std::endl;
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

