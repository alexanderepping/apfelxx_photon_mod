// used to write files
#include<fstream>

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include "apfel/apfelxx.h"

// Name of the used LHAPDF set
const std::string NameLHAPDFSet = "GRVCustomSet";

// Open LHAPDF set
LHAPDF::PDF* dist = LHAPDF::mkPDF(NameLHAPDFSet);

int main()
{
  // Final scale
  //double mu = 100;
  double mu = 10; 

  // used particle for output
  int flavor = 0;
  
  // flavor number with corresponding flavor
  std::map<int, std::string> mapFlavors {{-6, "tbar"}, {-5, "bbar"}, {-4, "cbar"}, {-3, "sbar"}, {-2, "ubar"}, {-1, "dbar"},
                                      {0, "gluon"}, {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}};
  
  // include Apfel or just show LHAPDF results?
  const bool   includeApfel = true;

  // Retrieve evolution parameters from the LHAPDF set
  const int    pto          = dist->orderQCD();
  const double Qref         = 91.1876;
  const double asref        = dist->alphasQ(Qref);
  const double mc           = dist->quarkThreshold(4);
  const double mb           = dist->quarkThreshold(5);
  const double mt           = dist->quarkThreshold(6);
  const double Qin          = dist->qMin();


  if ( includeApfel ) {
    /////////////////////
    // Initialise APFEL++
    /////////////////////

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

    // Call evolved PDFs at "mu". Notice that APFEL++ returns PDFs in
    // the QCD evolution basis, therefore one needs to rotate them back
    // to the physical basis which is done through the function
    // "QCDEvToPhys".
    const std::map<int, apfel::Distribution> tpdfs = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects());

    // Print results
    std::cout << std::scientific;

    // Print alphas from APFEL++, and LHAPDF
    std::cout << "\nmu = " << mu << " GeV\n" << std::endl;
    std::cout << "APFEL++: AlphaQCD(Q) = " << Alphas.Evaluate(mu) << std::endl;
    std::cout << "LHAPDF:  AlphaQCD(Q) = " << dist->alphasQ(mu) << std::endl;

    // Vector of test values of x
    std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
    std::vector<double> xlha2{1e-09, 2e-09, 3.0000000000000004e-09, 4e-09, 5e-09, 6.000000000000001e-09, 7.000000000000001e-09, 8e-09, 9.000000000000001e-09, 1e-08, 2e-08, 3.0000000000000004e-08, 4e-08, 5e-08, 6.000000000000001e-08, 7e-08, 8e-08, 9e-08, 1e-07, 2e-07, 3e-07, 4e-07, 5e-07, 6e-07, 7e-07, 8e-07, 9e-07, 1e-06, 2e-06, 3e-06, 4e-06, 4.9999999999999996e-06, 6e-06, 7e-06, 8e-06, 9e-06, 1e-05, 2e-05, 3.0000000000000004e-05, 4e-05, 5e-05, 6.000000000000001e-05, 7.000000000000001e-05, 8e-05, 9e-05, 0.0001, 0.0002, 0.00030000000000000003, 0.0004, 0.0005, 0.0006000000000000001, 0.0007, 0.0008, 0.0009000000000000001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009000000000000001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6000000000000001, 0.7000000000000001, 0.8, 0.9, };

    //addition ->
    std::ofstream file;
    file.open("/local0/a_eppi01/apfelxx_photon_mod/plottingPython/data_Evolution.txt");//addition 
    file << mapFlavors.at(flavor) << std::endl;
    for (double x : xlha2)
    {
      file << x << "\t" << tpdfs.at(flavor).Evaluate(x) << "\t" << dist->xfxQ(flavor, x, mu) << std::endl;
    }
    file.close();
    //addition <-a

    // Print down PDF at "mu"
    std::cout << "\n     x             APFEL++           LHAPDF" << std::endl;
    for (double x : xlha)
      std::cout << x << "\t" << tpdfs.at(flavor).Evaluate(x) << "\t" << dist->xfxQ(flavor, x, mu) << std::endl;
        
    std::cout << "\n";
    std::cout << "\nUsed LHAPDF Set           : " << NameLHAPDFSet;
    std::cout << "\nUsed Perturbative Order   : " << std::to_string(pto);
    std::cout << "\nUsed Mu                   : " << std::to_string(mu);
    std::cout << "\nShown particle:           : " << mapFlavors.at(flavor) << "\n";
    

    

  }
  else {
    // Print results
    std::cout << std::scientific;

    // Print alphas from LHAPDF
    std::cout << "\nmu = " << mu << " GeV\n" << std::endl;
    std::cout << "LHAPDF:  AlphaQCD(Q) = " << dist->alphasQ(mu) << std::endl;

    // Vector of test values of x
    //std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
    std::vector<double> xlha{1e-09, 2e-09, 3.0000000000000004e-09, 4e-09, 5e-09, 6.000000000000001e-09, 7.000000000000001e-09, 8e-09, 9.000000000000001e-09, 1e-08, 2e-08, 3.0000000000000004e-08, 4e-08, 5e-08, 6.000000000000001e-08, 7e-08, 8e-08, 9e-08, 1e-07, 2e-07, 3e-07, 4e-07, 5e-07, 6e-07, 7e-07, 8e-07, 9e-07, 1e-06, 2e-06, 3e-06, 4e-06, 4.9999999999999996e-06, 6e-06, 7e-06, 8e-06, 9e-06, 1e-05, 2e-05, 3.0000000000000004e-05, 4e-05, 5e-05, 6.000000000000001e-05, 7.000000000000001e-05, 8e-05, 9e-05, 0.0001, 0.0002, 0.00030000000000000003, 0.0004, 0.0005, 0.0006000000000000001, 0.0007, 0.0008, 0.0009000000000000001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009000000000000001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6000000000000001, 0.7000000000000001, 0.8, 0.9, };

    // Print down PDF at "mu"
    std::cout << "\n     x           LHAPDF" << std::endl;
    for (double x : xlha)
        std::cout << x << "\t" << dist->xfxQ(flavor, x, mu) << std::endl;
    std::cout << "\n";
    std::cout << "\nUsed LHAPDF Set           : " << NameLHAPDFSet;
    std::cout << "\nUsed Perturbative Order   : " << std::to_string(pto);
    std::cout << "\nUsed Mu                   : " << std::to_string(mu);
    std::cout << "\nShown particle:           : " << mapFlavors.at(flavor) << "\n";
  }

  return 0;
}

