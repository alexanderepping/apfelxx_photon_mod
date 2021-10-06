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
  double mu = 100;

  // Retrieve evolution parameters from the LHAPDF set
  const int    pto          = dist->orderQCD();
  const double Qref         = 91.1876;
  const double asref        = dist->alphasQ(Qref);
  const double mc           = dist->quarkThreshold(4);
  const double mb           = dist->quarkThreshold(5);
  const double mt           = dist->quarkThreshold(6);
  const double Qin          = dist->qMin();
  const bool   includeApfel = true;


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
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3};

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

    // Print down PDF at "mu"
    std::cout << "\n     x             APFEL++           LHAPDF" << std::endl;
    for (double x : xlha)
        std::cout << x << "\t" << tpdfs.at(1).Evaluate(x) << "\t" << dist->xfxQ(1, x, mu) << std::endl;
    std::cout << "\n";
    std::cout <<  "\nUsed LHAPDF Set           : "+NameLHAPDFSet+"\nUsed Perturbative Order   : "+std::to_string(pto)+"\n";
  }
  else {
    // Print results
    std::cout << std::scientific;

    // Print alphas from LHAPDF
    std::cout << "\nmu = " << mu << " GeV\n" << std::endl;
    std::cout << "LHAPDF:  AlphaQCD(Q) = " << dist->alphasQ(mu) << std::endl;

    // Vector of test values of x
    std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

    // Print down PDF at "mu"
    std::cout << "\n     x           LHAPDF" << std::endl;
    for (double x : xlha)
        std::cout << x << "\t" << dist->xfxQ(1, x, mu) << std::endl;
    std::cout << "\n";
    std::cout <<  "\nUsed LHAPDF Set           : "+NameLHAPDFSet+"\nUsed Perturbative Order   : "+std::to_string(pto)+"\n";
  }

  return 0;
}

