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
  const int    pto   = dist->orderQCD();
  const double Qref  = 91.1876;
  const double asref = dist->alphasQ(Qref);
  const double mc    = dist->quarkThreshold(4);
  const double mb    = dist->quarkThreshold(5);
  const double mt    = dist->quarkThreshold(6);
  const double Qin   = dist->qMin();

  /////////////////////
  // Initialise APFEL++
  /////////////////////

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

  // Initialise F2 NC structure function in the ZM-VFNS. Functions
  // exist for FL and xF3 (btw, APFEL produces results for xF3, not
  // F3). Also CC structure functions are implemented in APFEL++. Also
  // the massive scheme and its massless limit are supported. At the
  // moment CC massive structure functions aren't there but can be
  // easily implemented if needed.
  const auto F2 = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, pto, as, fBq);

  // Tabulate the ftructure function over 50 intervals in the range
  // [1: 1000] GeV using interpolation degree 3. Take into account
  // that at thresholds straucture functions may not be continuos
  // feeding the contructor with the thresholds.
  const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  // Print results
  std::cout << std::scientific;

  // Vector of test values of x
  std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  // Print structure function
  std::cout << "     x             APFEL++" << std::endl;
  for (double x : xlha)
      std::cout << x << "\t" << F2total.EvaluatexQ(x, mu) << std::endl;
  std::cout << "\n";
    std::cout <<  "\nUsed LHAPDF Set           : "+NameLHAPDFSet+"\nUsed Perturbative Order   : "+std::to_string(pto)+"\n";

  return 0;
}
