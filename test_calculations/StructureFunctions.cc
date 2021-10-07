// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include "apfel/apfelxx.h"

// APFEL libs
#include "apfel/APFEL.h"

// Open LHAPDF set
LHAPDF::PDF* dist = LHAPDF::mkPDF("CT14nnlo");

// Function used by APFEL as an input. It can be used to implement any
// input set of PDFs.
extern "C" void externalsetapfel_(double const& x, double const& Q, double* xf)
{
  // Call PDFs from LHAPDF
  const std::map<int, double> xpdf = dist->xfxQ(x, Q);

  // The array "xf" is defined of the range [0,...,12] corresponding to
  // tbar, bbar, ..., g, ..., b, t.
  for (auto f : xpdf)
    xf[(f.first == 21 ? 6 : f.first + 6)] = f.second;
}

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

  // Vector of threholds. Not including the top thesholds
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

  /////////////////////
  // Initialise APFEL
  /////////////////////

  // I hope that the name of the setting functions are
  // self-explanatory. One remark on the "SetPDFSet" function: the
  // particular argument "external" instructs the code to call the
  // function "externalsetapfel_" (with that precise signature)
  // defined above. That function can thus be used to implement any
  // set of functions to be evolved without the need to hack APFEL.
  APFEL::SetMassScheme("ZM-VFNS");
  APFEL::SetProcessDIS("NC");
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetPDFSet("external");
  //APFEL::SetPDFSet("CT14nnlo"); // This should give identical results
  APFEL::SetPerturbativeOrder(pto);
  APFEL::SetAlphaQCDRef(dist->alphasQ(mu), mu); // Make sure to be using the same value of alphas at "mu"
  APFEL::SetPoleMasses(mc, mb, mt);

  // Intialise structure functions and compute them at "mu". Notice
  // that the function "ComputeStructureFunctionsAPFEL" takes two
  // arguments corresponding to initial and final scales. If they are
  // different, APFEL performs the evolution between them and finally
  // computes the structure functions at the final scale. If they are
  // set equal no evolution is performed and PDFs are callred directly
  // at the scale "mu".
  APFEL::InitializeAPFEL_DIS();
  APFEL::ComputeStructureFunctionsAPFEL(mu, mu);

  // Print results
  std::cout << std::scientific;

  // Vector of test values of x
  std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  // Print structure function
  std::cout << "     x             APFEL++         APFEL" << std::endl;
  for (double x : xlha)
      std::cout << x << "\t" << F2total.EvaluatexQ(x, mu) << "\t" << APFEL::F2total(x) << std::endl;
  std::cout << "\n";

  return 0;
}
