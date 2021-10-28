/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* EvolutionStructureFunctions.cc                                                                                          */
/* modified by: Alexander Epping: a_eppi01@uni-muenster.de                                                                 */
/* GitHub: https://github.com/alexanderepping/apfelxx_photon_mod                                                           */
/* 25 Oct 2021                                                                                                             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* original file: https://github.com/vbertone/APFEL_Examples                                                               */
/* original author: Valerio Bertone: valerio.bertone@cern.ch                                                               */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Program to evolve the PDFs with given initial conditions and print out the results in a LHAPDF format to be used by     */
/* the StructureFunctions.cc program.                                                                                      */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program takes a LHAPDF set as input and evolves the PDFs using the DGLAP evolution, starting at an initial         */
/* energy. The evolved PDFs are then output to a filein the LHAPDF format.                                                 */
/*                                                                                                                         */
/* The output file is then used by StructureFunctions.cc                                                                   */
/*                                                                                                                         */
/* This program and the StructureFunctions.cc program can be run by using the run_StructureFunctions.sh file               */
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
// definitions, can be changed by user
///////////////////////////////////////

// name of the used LHAPDF set
const std::string NameLHAPDFSet = "GRVCustomSetLO";

// Name of the output file
const std::string OutputFile = "/home/alexander/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/EvolvedGRVCustomSetLO/EvolvedGRVCustomSetLO_0000.dat";



///////////////////////////////////////
// functions
///////////////////////////////////////

// define cutoff function
double Cutoff (double const& val)
{
  if (val > 1.0e-08)
    return val;
  else
    return 0.000000e-00;
}



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

  // print introduction data to file
  std::string lhapdfFormalia = "PdfType: central\nFormat: lhagrid1\n---\n1.000000e-09 1.131190e-09 1.381630e-09 1.687530e-09 2.061150e-09 2.517500e-09 3.074880e-09 3.755670e-09 4.587180e-09 5.602800e-09 6.843270e-09 8.358390e-09 1.020900e-08 1.246930e-08 1.523000e-08 1.860190e-08 2.272050e-08 2.775080e-08 3.389490e-08 4.139940e-08 5.056530e-08 6.176060e-08 7.543460e-08 9.213600e-08 1.125350e-07 1.374510e-07 1.678830e-07 2.050520e-07 2.504520e-07 3.059020e-07 3.736300e-07 4.563530e-07 5.573900e-07 6.807980e-07 8.315290e-07 1.015630e-06 1.240500e-06 1.515140e-06 1.850600e-06 2.260330e-06 2.760770e-06 3.372020e-06 4.118590e-06 5.030460e-06 6.144210e-06 7.504560e-06 9.166090e-06 1.119550e-05 1.367420e-05 1.670170e-05 2.039950e-05 2.491600e-05 3.043250e-05 3.717030e-05 4.539990e-05 5.545160e-05 6.772870e-05 8.272410e-05 1.010390e-04 1.234100e-04 1.507330e-04 1.841060e-04 2.248670e-04 2.746540e-04 3.354630e-04 4.097350e-04 5.004510e-04 6.112530e-04 7.465860e-04 9.118820e-04 1.113780e-03 1.360370e-03 1.661560e-03 2.029430e-03 2.478750e-03 3.027550e-03 3.697860e-03 4.516580e-03 5.516560e-03 6.737950e-03 8.229750e-03 1.005180e-02 1.116455e-02 1.227730e-02 1.363645e-02 1.499560e-02 1.665560e-02 1.831560e-02 2.034320e-02 2.237080e-02 2.484725e-02 2.732370e-02 3.034850e-02 3.337330e-02 3.706775e-02 4.076220e-02 4.527465e-02 4.978710e-02 5.529860e-02 6.081010e-02 6.754185e-02 7.427360e-02 8.249580e-02 9.071800e-02 1.007605e-01 1.108030e-01 1.230690e-01 1.353350e-01 1.400000e-01 1.446650e-01 1.496515e-01 1.546380e-01 1.599685e-01 1.652990e-01 1.709965e-01 1.766940e-01 1.827850e-01 1.888760e-01 1.953865e-01 2.018970e-01 2.088560e-01 2.158150e-01 2.232540e-01 2.306930e-01 2.386450e-01 2.465970e-01 2.550970e-01 2.635970e-01 2.726830e-01 2.817690e-01 2.914815e-01 3.011940e-01 3.115760e-01 3.219580e-01 3.330560e-01 3.441540e-01 3.560165e-01 3.678790e-01 3.805600e-01 3.932410e-01 4.067955e-01 4.203500e-01 4.348395e-01 4.493290e-01 4.648170e-01 4.803050e-01 4.968610e-01 5.134170e-01 5.311145e-01 5.488120e-01 5.549780e-01 5.611440e-01 5.674485e-01 5.737530e-01 5.801995e-01 5.866460e-01 5.932375e-01 5.998290e-01 6.065680e-01 6.133070e-01 6.201980e-01 6.270890e-01 6.341345e-01 6.411800e-01 6.483840e-01 6.555880e-01 6.629540e-01 6.703200e-01 6.778515e-01 6.853830e-01 6.930835e-01 7.007840e-01 7.086575e-01 7.165310e-01 7.245815e-01 7.326320e-01 7.408635e-01 7.490950e-01 7.575115e-01 7.659280e-01 7.745335e-01 7.831390e-01 7.919380e-01 8.007370e-01 8.097340e-01 8.187310e-01 8.217745e-01 8.248180e-01 8.278840e-01 8.309500e-01 8.340390e-01 8.371280e-01 8.402400e-01 8.433520e-01 8.464875e-01 8.496230e-01 8.527815e-01 8.559400e-01 8.591215e-01 8.623030e-01 8.655090e-01 8.687150e-01 8.719440e-01 8.751730e-01 8.784265e-01 8.816800e-01 8.849575e-01 8.882350e-01 8.915370e-01 8.948390e-01 8.981655e-01 9.014920e-01 9.048435e-01 9.081950e-01 9.115710e-01 9.149470e-01 9.183485e-01 9.217500e-01 9.251765e-01 9.286030e-01 9.320550e-01 9.355070e-01 9.389845e-01 9.424620e-01 9.459660e-01 9.494700e-01 9.529995e-01 9.565290e-01 9.600845e-01 9.636400e-01 9.672225e-01 9.708050e-01 9.744140e-01 9.780230e-01 9.816585e-01 9.852940e-01 9.889570e-01 9.926200e-01 9.963100e-01 1.000000e+00\n1.295000e+00 1.298750e+00 1.464610e+00 1.659240e+00 1.890670e+00 2.167490e+00 2.500670e+00 2.904300e+00 3.396640e+00 4.001450e+00 4.750000e+00 5.767150e+00 7.070720e+00 8.758190e+00 1.096570e+01 1.388560e+01 1.779290e+01 2.308550e+01 3.034710e+01 4.044480e+01 5.468640e+01 7.507240e+01 1.047120e+02 1.485170e+02 2.143800e+02 3.152120e+02 4.725370e+02 7.229460e+02 1.129950e+03 1.806160e+03 2.955930e+03 4.958860e+03 8.538140e+03 1.510790e+04 2.751070e+04 5.162750e+04 1.000000e+05\n-5 -4 -3 -2 -1 1 2 3 4 5 21";
  file << lhapdfFormalia << std::endl;



  /////////////////////////////////////
  // definitions
  /////////////////////////////////////
  
  // Vector of test values of x
  std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  // array of final scale values for which data should be output; taken from makeLHAPDF/LHAPDF_template.dat
  double muLHAPDF[] = {1.295000e+00, 1.298750e+00, 1.464610e+00, 1.659240e+00, 1.890670e+00, 2.167490e+00, 2.500670e+00, 2.904300e+00, 3.396640e+00, 4.001450e+00, 4.750000e+00, 5.767150e+00, 7.070720e+00, 8.758190e+00, 1.096570e+01, 1.388560e+01, 1.779290e+01, 2.308550e+01, 3.034710e+01, 4.044480e+01, 5.468640e+01, 7.507240e+01, 1.047120e+02, 1.485170e+02, 2.143800e+02, 3.152120e+02, 4.725370e+02, 7.229460e+02, 1.129950e+03, 1.806160e+03, 2.955930e+03, 4.958860e+03, 8.538140e+03, 1.510790e+04, 2.751070e+04, 5.162750e+04, 1.000000e+05}; 


  // array of x values for which data should be output; taken from makeLHAPDF/LHAPDF_template.dat
  std::vector<double> xLHAPDF{1.000000e-09, 1.131190e-09, 1.381630e-09, 1.687530e-09, 2.061150e-09, 2.517500e-09, 3.074880e-09, 3.755670e-09, 4.587180e-09, 5.602800e-09, 6.843270e-09, 8.358390e-09, 1.020900e-08, 1.246930e-08, 1.523000e-08, 1.860190e-08, 2.272050e-08, 2.775080e-08, 3.389490e-08, 4.139940e-08, 5.056530e-08, 6.176060e-08, 7.543460e-08, 9.213600e-08, 1.125350e-07, 1.374510e-07, 1.678830e-07, 2.050520e-07, 2.504520e-07, 3.059020e-07, 3.736300e-07, 4.563530e-07, 5.573900e-07, 6.807980e-07, 8.315290e-07, 1.015630e-06, 1.240500e-06, 1.515140e-06, 1.850600e-06, 2.260330e-06, 2.760770e-06, 3.372020e-06, 4.118590e-06, 5.030460e-06, 6.144210e-06, 7.504560e-06, 9.166090e-06, 1.119550e-05, 1.367420e-05, 1.670170e-05, 2.039950e-05, 2.491600e-05, 3.043250e-05, 3.717030e-05, 4.539990e-05, 5.545160e-05, 6.772870e-05, 8.272410e-05, 1.010390e-04, 1.234100e-04, 1.507330e-04, 1.841060e-04, 2.248670e-04, 2.746540e-04, 3.354630e-04, 4.097350e-04, 5.004510e-04, 6.112530e-04, 7.465860e-04, 9.118820e-04, 1.113780e-03, 1.360370e-03, 1.661560e-03, 2.029430e-03, 2.478750e-03, 3.027550e-03, 3.697860e-03, 4.516580e-03, 5.516560e-03, 6.737950e-03, 8.229750e-03, 1.005180e-02, 1.116455e-02, 1.227730e-02, 1.363645e-02, 1.499560e-02, 1.665560e-02, 1.831560e-02, 2.034320e-02, 2.237080e-02, 2.484725e-02, 2.732370e-02, 3.034850e-02, 3.337330e-02, 3.706775e-02, 4.076220e-02, 4.527465e-02, 4.978710e-02, 5.529860e-02, 6.081010e-02, 6.754185e-02, 7.427360e-02, 8.249580e-02, 9.071800e-02, 1.007605e-01, 1.108030e-01, 1.230690e-01, 1.353350e-01, 1.400000e-01, 1.446650e-01, 1.496515e-01, 1.546380e-01, 1.599685e-01, 1.652990e-01, 1.709965e-01, 1.766940e-01, 1.827850e-01, 1.888760e-01, 1.953865e-01, 2.018970e-01, 2.088560e-01, 2.158150e-01, 2.232540e-01, 2.306930e-01, 2.386450e-01, 2.465970e-01, 2.550970e-01, 2.635970e-01, 2.726830e-01, 2.817690e-01, 2.914815e-01, 3.011940e-01, 3.115760e-01, 3.219580e-01, 3.330560e-01, 3.441540e-01, 3.560165e-01, 3.678790e-01, 3.805600e-01, 3.932410e-01, 4.067955e-01, 4.203500e-01, 4.348395e-01, 4.493290e-01, 4.648170e-01, 4.803050e-01, 4.968610e-01, 5.134170e-01, 5.311145e-01, 5.488120e-01, 5.549780e-01, 5.611440e-01, 5.674485e-01, 5.737530e-01, 5.801995e-01, 5.866460e-01, 5.932375e-01, 5.998290e-01, 6.065680e-01, 6.133070e-01, 6.201980e-01, 6.270890e-01, 6.341345e-01, 6.411800e-01, 6.483840e-01, 6.555880e-01, 6.629540e-01, 6.703200e-01, 6.778515e-01, 6.853830e-01, 6.930835e-01, 7.007840e-01, 7.086575e-01, 7.165310e-01, 7.245815e-01, 7.326320e-01, 7.408635e-01, 7.490950e-01, 7.575115e-01, 7.659280e-01, 7.745335e-01, 7.831390e-01, 7.919380e-01, 8.007370e-01, 8.097340e-01, 8.187310e-01, 8.217745e-01, 8.248180e-01, 8.278840e-01, 8.309500e-01, 8.340390e-01, 8.371280e-01, 8.402400e-01, 8.433520e-01, 8.464875e-01, 8.496230e-01, 8.527815e-01, 8.559400e-01, 8.591215e-01, 8.623030e-01, 8.655090e-01, 8.687150e-01, 8.719440e-01, 8.751730e-01, 8.784265e-01, 8.816800e-01, 8.849575e-01, 8.882350e-01, 8.915370e-01, 8.948390e-01, 8.981655e-01, 9.014920e-01, 9.048435e-01, 9.081950e-01, 9.115710e-01, 9.149470e-01, 9.183485e-01, 9.217500e-01, 9.251765e-01, 9.286030e-01, 9.320550e-01, 9.355070e-01, 9.389845e-01, 9.424620e-01, 9.459660e-01, 9.494700e-01, 9.529995e-01, 9.565290e-01, 9.600845e-01, 9.636400e-01, 9.672225e-01, 9.708050e-01, 9.744140e-01, 9.780230e-01, 9.816585e-01, 9.852940e-01, 9.889570e-01, 9.926200e-01, 9.963100e-01, 1.000000e+00};

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

  std::cout << std::scientific;
  for (double x : xLHAPDF)
  {
    //file << "\nx=" << x << std::endl; //debug

    // file output:
    // Print down PDF at "mu"
    for (double mu : muLHAPDF)
    {
      std::cout << std::scientific;
      // Call evolved PDFs at "mu". Notice that APFEL++ returns PDFs in
      // the QCD evolution basis, therefore one needs to rotate them back
      // to the physical basis which is done through the function
      // "QCDEvToPhys".
      const std::map<int, apfel::Distribution> tpdfs = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects());



      //file << "mu=" << mu << "    " << Cutoff(tpdfs.at(-5).Evaluate(x)) << " " << Cutoff(tpdfs.at(-4).Evaluate(x)) << " " << Cutoff(tpdfs.at(-3).Evaluate(x)) << " " << Cutoff(tpdfs.at(-2).Evaluate(x)) << " " << Cutoff(tpdfs.at(-1).Evaluate(x)) << " " << Cutoff(tpdfs.at(1).Evaluate(x)) << " " << Cutoff(tpdfs.at(2).Evaluate(x)) << " " << Cutoff(tpdfs.at(3).Evaluate(x)) << " " << Cutoff(tpdfs.at(4).Evaluate(x)) << " " << Cutoff(tpdfs.at(5).Evaluate(x)) << " " << Cutoff(tpdfs.at(0).Evaluate(x)) << std::endl; //debug

      file << Cutoff(tpdfs.at(-5).Evaluate(x)) << " " << Cutoff(tpdfs.at(-4).Evaluate(x)) << " " << Cutoff(tpdfs.at(-3).Evaluate(x)) << " " << Cutoff(tpdfs.at(-2).Evaluate(x)) << " " << Cutoff(tpdfs.at(-1).Evaluate(x)) << " " << Cutoff(tpdfs.at(1).Evaluate(x)) << " " << Cutoff(tpdfs.at(2).Evaluate(x)) << " " << Cutoff(tpdfs.at(3).Evaluate(x)) << " " << Cutoff(tpdfs.at(4).Evaluate(x)) << " " << Cutoff(tpdfs.at(5).Evaluate(x)) << " " << Cutoff(tpdfs.at(0).Evaluate(x)) << std::endl;
    }
  }
  
  // close output file
  file << "---" << std::endl;
  file.close();

  return 0;
}

