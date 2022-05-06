/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* StructureFunctions.cc                                                                                                   */
/* modified by: Alexander Epping: a_eppi01@uni-muenster.de                                                                 */
/* GitHub: https://github.com/alexanderepping/apfelxx_photon_mod                                                           */
/* 11 Nov 2021                                                                                                             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Program to calculate the Structure Functions from a given LHAPDF set.                                                   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program takes an LHAPDF set as input, evolves the PDFs, prints them to a file and then calculates the Structure    */
/* Functions from the evolved and not-evolved PDFs. Their values are then output to a file and the terminal.               */
/*                                                                                                                         */
/* The output file for the LHAPDF data is saved in the respective LHAPDF share directory                                   */
/* The output file for the Structure Functions is saved in plottingPython/ and can be used by plottingStructureFunctions.  */
/*                                                                                                                         */
/* This program and the plotting program can be run by using the run_StructureFunctions.sh file in the bashFiles directory.*/
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



////////////////////////////////////////////////////////////
// imports
////////////////////////////////////////////////////////////

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include "apfel/apfelxx.h" 

// used to write to files
#include<fstream>



////////////////////////////////////////////////////////////
// definitions, can be changed by user
////////////////////////////////////////////////////////////

// name of the used LHAPDF set
const std::string NameLHAPDFSet = "GRVCustomSetHO";

#define asGRV

// Name of the output file for the Evolved PDFs
const std::string OutputFileEvolved = "/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/Evolved"+NameLHAPDFSet+"/Evolved"+NameLHAPDFSet+"_0000.dat";

// Name of the output file for the Structure Functions of the Evolved PDFs
const std::string OutputFileStructureFunctions = "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/data_StructureFunctions.txt";



////////////////////////////////////////////////////////////
// defintions
////////////////////////////////////////////////////////////

// array of final scale values for which data should be output (mu^2)
//                 |  0,1,2:ALEPH  |  3,4,5:AMY
const double arr_mu2[] = {9.9, 20.7, 284, 6.8, 73, 390};   

// xlha is used for console output of the structure function data
std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

// xlha2 is used for file output and plotting of the structure function data
std::vector<double> xlha2{0.0001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110, 0.115, 0.120, 0.125, 0.130, 0.135, 0.140, 0.145, 0.150, 0.155, 0.160, 0.165, 0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.200, 0.205, 0.210, 0.215, 0.220, 0.225, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.300, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355, 0.360, 0.365, 0.370, 0.375, 0.380, 0.385, 0.390, 0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425, 0.430, 0.435, 0.440, 0.445, 0.450, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505, 0.510, 0.515, 0.520, 0.525, 0.530, 0.535, 0.540, 0.545, 0.550, 0.555, 0.560, 0.565, 0.570, 0.575, 0.580, 0.585, 0.590, 0.595, 0.600, 0.605, 0.610, 0.615, 0.620, 0.625, 0.630, 0.635, 0.640, 0.645, 0.650, 0.655, 0.660, 0.665, 0.670, 0.675, 0.680, 0.685, 0.690, 0.695, 0.700, 0.705, 0.710, 0.715, 0.720, 0.725, 0.730, 0.735, 0.740, 0.745, 0.750, 0.755, 0.760, 0.765, 0.770, 0.775, 0.780, 0.785, 0.790, 0.795, 0.800, 0.805, 0.810, 0.815, 0.820, 0.825, 0.830, 0.835, 0.840, 0.845, 0.850, 0.855, 0.860, 0.865, 0.870, 0.875, 0.880, 0.885, 0.890, 0.895, 0.900};//, 0.905, 0.910, 0.915, 0.920, 0.925, 0.930, 0.935, 0.940, 0.945, 0.950, 0.955, 0.960, 0.965, 0.970, 0.975, 0.980, 0.985, 0.990, 0.995}; 

// array of final scale values for which data of the Evolved PDFs should be output; taken from makeLHAPDF/LHAPDF_template.dat
double muLHAPDF[] = {1.295000e+00, 1.298750e+00, 1.464610e+00, 1.659240e+00, 1.890670e+00, 2.167490e+00, 2.500670e+00, 2.904300e+00, 3.396640e+00, 4.001450e+00, 4.750000e+00, 5.767150e+00, 7.070720e+00, 8.758190e+00, 1.096570e+01, 1.388560e+01, 1.779290e+01, 2.308550e+01, 3.034710e+01, 4.044480e+01, 5.468640e+01, 7.507240e+01, 1.047120e+02, 1.485170e+02, 2.143800e+02, 3.152120e+02, 4.725370e+02, 7.229460e+02, 1.129950e+03, 1.806160e+03, 2.955930e+03, 4.958860e+03, 8.538140e+03, 1.510790e+04, 2.751070e+04, 5.162750e+04, 1.000000e+05}; 

// array of x values for which data of the Evolved PDFs should be output; taken from makeLHAPDF/LHAPDF_template.dat
std::vector<double> xLHAPDF{1.000000e-09, 1.131190e-09, 1.381630e-09, 1.687530e-09, 2.061150e-09, 2.517500e-09, 3.074880e-09, 3.755670e-09, 4.587180e-09, 5.602800e-09, 6.843270e-09, 8.358390e-09, 1.020900e-08, 1.246930e-08, 1.523000e-08, 1.860190e-08, 2.272050e-08, 2.775080e-08, 3.389490e-08, 4.139940e-08, 5.056530e-08, 6.176060e-08, 7.543460e-08, 9.213600e-08, 1.125350e-07, 1.374510e-07, 1.678830e-07, 2.050520e-07, 2.504520e-07, 3.059020e-07, 3.736300e-07, 4.563530e-07, 5.573900e-07, 6.807980e-07, 8.315290e-07, 1.015630e-06, 1.240500e-06, 1.515140e-06, 1.850600e-06, 2.260330e-06, 2.760770e-06, 3.372020e-06, 4.118590e-06, 5.030460e-06, 6.144210e-06, 7.504560e-06, 9.166090e-06, 1.119550e-05, 1.367420e-05, 1.670170e-05, 2.039950e-05, 2.491600e-05, 3.043250e-05, 3.717030e-05, 4.539990e-05, 5.545160e-05, 6.772870e-05, 8.272410e-05, 1.010390e-04, 1.234100e-04, 1.507330e-04, 1.841060e-04, 2.248670e-04, 2.746540e-04, 3.354630e-04, 4.097350e-04, 5.004510e-04, 6.112530e-04, 7.465860e-04, 9.118820e-04, 1.113780e-03, 1.360370e-03, 1.661560e-03, 2.029430e-03, 2.478750e-03, 3.027550e-03, 3.697860e-03, 4.516580e-03, 5.516560e-03, 6.737950e-03, 8.229750e-03, 1.005180e-02, 1.116455e-02, 1.227730e-02, 1.363645e-02, 1.499560e-02, 1.665560e-02, 1.831560e-02, 2.034320e-02, 2.237080e-02, 2.484725e-02, 2.732370e-02, 3.034850e-02, 3.337330e-02, 3.706775e-02, 4.076220e-02, 4.527465e-02, 4.978710e-02, 5.529860e-02, 6.081010e-02, 6.754185e-02, 7.427360e-02, 8.249580e-02, 9.071800e-02, 1.007605e-01, 1.108030e-01, 1.230690e-01, 1.353350e-01, 1.400000e-01, 1.446650e-01, 1.496515e-01, 1.546380e-01, 1.599685e-01, 1.652990e-01, 1.709965e-01, 1.766940e-01, 1.827850e-01, 1.888760e-01, 1.953865e-01, 2.018970e-01, 2.088560e-01, 2.158150e-01, 2.232540e-01, 2.306930e-01, 2.386450e-01, 2.465970e-01, 2.550970e-01, 2.635970e-01, 2.726830e-01, 2.817690e-01, 2.914815e-01, 3.011940e-01, 3.115760e-01, 3.219580e-01, 3.330560e-01, 3.441540e-01, 3.560165e-01, 3.678790e-01, 3.805600e-01, 3.932410e-01, 4.067955e-01, 4.203500e-01, 4.348395e-01, 4.493290e-01, 4.648170e-01, 4.803050e-01, 4.968610e-01, 5.134170e-01, 5.311145e-01, 5.488120e-01, 5.549780e-01, 5.611440e-01, 5.674485e-01, 5.737530e-01, 5.801995e-01, 5.866460e-01, 5.932375e-01, 5.998290e-01, 6.065680e-01, 6.133070e-01, 6.201980e-01, 6.270890e-01, 6.341345e-01, 6.411800e-01, 6.483840e-01, 6.555880e-01, 6.629540e-01, 6.703200e-01, 6.778515e-01, 6.853830e-01, 6.930835e-01, 7.007840e-01, 7.086575e-01, 7.165310e-01, 7.245815e-01, 7.326320e-01, 7.408635e-01, 7.490950e-01, 7.575115e-01, 7.659280e-01, 7.745335e-01, 7.831390e-01, 7.919380e-01, 8.007370e-01, 8.097340e-01, 8.187310e-01, 8.217745e-01, 8.248180e-01, 8.278840e-01, 8.309500e-01, 8.340390e-01, 8.371280e-01, 8.402400e-01, 8.433520e-01, 8.464875e-01, 8.496230e-01, 8.527815e-01, 8.559400e-01, 8.591215e-01, 8.623030e-01, 8.655090e-01, 8.687150e-01, 8.719440e-01, 8.751730e-01, 8.784265e-01, 8.816800e-01, 8.849575e-01, 8.882350e-01, 8.915370e-01, 8.948390e-01, 8.981655e-01, 9.014920e-01, 9.048435e-01, 9.081950e-01, 9.115710e-01, 9.149470e-01, 9.183485e-01, 9.217500e-01, 9.251765e-01, 9.286030e-01, 9.320550e-01, 9.355070e-01, 9.389845e-01, 9.424620e-01, 9.459660e-01, 9.494700e-01, 9.529995e-01, 9.565290e-01, 9.600845e-01, 9.636400e-01, 9.672225e-01, 9.708050e-01, 9.744140e-01, 9.780230e-01, 9.816585e-01, 9.852940e-01, 9.889570e-01, 9.926200e-01, 9.963100e-01, 1.000000e+00};

// flavor number with corresponding flavor
std::map<int, std::string> mapFlavors {{-6, "tbar"}, {-5, "bbar"}, {-4, "cbar"}, {-3, "sbar"}, {-2, "ubar"}, {-1, "dbar"},
                                      {0, "gluon"}, {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}};


////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////

// define cutoff function
double Cutoff (double const& val)
{
  if (val > 1.0e-08)
    return val;
  else
    return 0.000000e-00;
}



////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////
int main()
{     
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



    //////////////////////////////////////////////////////////
    // Evolve PDFs
    //////////////////////////////////////////////////////////

    // Define Grid
    const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 6e-1, 3}, apfel::SubGrid{50, 8e-1, 5}}};


    // Define Thresholds and Masses
    const std::vector<double> Thresholds = {0, 0, 0, mc, mb};
    const std::vector<double> Masses = Thresholds;


    // Calculate AlphaS(Q)
    apfel::AlphaQCD a{asref, Qref, Masses, Thresholds, pto};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    // const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
//----calculate alphas(Q) using the evolution calculated by apfel----
  #ifdef asApfel
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
  #endif

  //----calculate Alphas(Q) using equation(2.11) given in                  ----
  //----Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01)----
  #ifdef asGRV
  const auto as = [&] (double const& mu) -> double
  {
    int nf;

    if (mu < 1.5)       nf = 3;
    else if (mu < 4.5)  nf = 4;
    else if (mu < 100.) nf = 5;
    else                nf = 5;


    double beta0 =  11 -  2. * nf / 3.; 
    double beta1 = 102 - 38. * nf / 3.;

    if (pto == 0) beta1 = 0.;


    std::map<int,std::vector<double>> lambdas = {{0, {0, 0, 0, 0.232, 0.200, 0.153, 0.082}}, {1, {0, 0, 0, 0.248, 0.200, 0.131, 0.053}}}; 

    double lnQ2Lambda2 = log((mu * mu) / (lambdas.at(pto)[nf]*lambdas.at(pto)[nf])); 

    return 4*M_PI/(beta0 * lnQ2Lambda2 + beta1 / beta0 * log(lnQ2Lambda2)); 
  };
  #endif


    // Define Input PDFs
    const auto InPDFs = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); };


    // Initialise and Build DglapObjects
    auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs, Qin, pto, as);


    // Tabulate PDFs
    //const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3, pto, as};



    //////////////////////////////////////////////////////////
    // print results of Evolution to file
    //////////////////////////////////////////////////////////

    // opening output file
    std::ofstream file;
    file.open(OutputFileEvolved);

    // print introduction data to file
    std::string lhapdfFormalia = "PdfType: central\nFormat: lhagrid1\n---\n1.000000e-09 1.131190e-09 1.381630e-09 1.687530e-09 2.061150e-09 2.517500e-09 3.074880e-09 3.755670e-09 4.587180e-09 5.602800e-09 6.843270e-09 8.358390e-09 1.020900e-08 1.246930e-08 1.523000e-08 1.860190e-08 2.272050e-08 2.775080e-08 3.389490e-08 4.139940e-08 5.056530e-08 6.176060e-08 7.543460e-08 9.213600e-08 1.125350e-07 1.374510e-07 1.678830e-07 2.050520e-07 2.504520e-07 3.059020e-07 3.736300e-07 4.563530e-07 5.573900e-07 6.807980e-07 8.315290e-07 1.015630e-06 1.240500e-06 1.515140e-06 1.850600e-06 2.260330e-06 2.760770e-06 3.372020e-06 4.118590e-06 5.030460e-06 6.144210e-06 7.504560e-06 9.166090e-06 1.119550e-05 1.367420e-05 1.670170e-05 2.039950e-05 2.491600e-05 3.043250e-05 3.717030e-05 4.539990e-05 5.545160e-05 6.772870e-05 8.272410e-05 1.010390e-04 1.234100e-04 1.507330e-04 1.841060e-04 2.248670e-04 2.746540e-04 3.354630e-04 4.097350e-04 5.004510e-04 6.112530e-04 7.465860e-04 9.118820e-04 1.113780e-03 1.360370e-03 1.661560e-03 2.029430e-03 2.478750e-03 3.027550e-03 3.697860e-03 4.516580e-03 5.516560e-03 6.737950e-03 8.229750e-03 1.005180e-02 1.116455e-02 1.227730e-02 1.363645e-02 1.499560e-02 1.665560e-02 1.831560e-02 2.034320e-02 2.237080e-02 2.484725e-02 2.732370e-02 3.034850e-02 3.337330e-02 3.706775e-02 4.076220e-02 4.527465e-02 4.978710e-02 5.529860e-02 6.081010e-02 6.754185e-02 7.427360e-02 8.249580e-02 9.071800e-02 1.007605e-01 1.108030e-01 1.230690e-01 1.353350e-01 1.400000e-01 1.446650e-01 1.496515e-01 1.546380e-01 1.599685e-01 1.652990e-01 1.709965e-01 1.766940e-01 1.827850e-01 1.888760e-01 1.953865e-01 2.018970e-01 2.088560e-01 2.158150e-01 2.232540e-01 2.306930e-01 2.386450e-01 2.465970e-01 2.550970e-01 2.635970e-01 2.726830e-01 2.817690e-01 2.914815e-01 3.011940e-01 3.115760e-01 3.219580e-01 3.330560e-01 3.441540e-01 3.560165e-01 3.678790e-01 3.805600e-01 3.932410e-01 4.067955e-01 4.203500e-01 4.348395e-01 4.493290e-01 4.648170e-01 4.803050e-01 4.968610e-01 5.134170e-01 5.311145e-01 5.488120e-01 5.549780e-01 5.611440e-01 5.674485e-01 5.737530e-01 5.801995e-01 5.866460e-01 5.932375e-01 5.998290e-01 6.065680e-01 6.133070e-01 6.201980e-01 6.270890e-01 6.341345e-01 6.411800e-01 6.483840e-01 6.555880e-01 6.629540e-01 6.703200e-01 6.778515e-01 6.853830e-01 6.930835e-01 7.007840e-01 7.086575e-01 7.165310e-01 7.245815e-01 7.326320e-01 7.408635e-01 7.490950e-01 7.575115e-01 7.659280e-01 7.745335e-01 7.831390e-01 7.919380e-01 8.007370e-01 8.097340e-01 8.187310e-01 8.217745e-01 8.248180e-01 8.278840e-01 8.309500e-01 8.340390e-01 8.371280e-01 8.402400e-01 8.433520e-01 8.464875e-01 8.496230e-01 8.527815e-01 8.559400e-01 8.591215e-01 8.623030e-01 8.655090e-01 8.687150e-01 8.719440e-01 8.751730e-01 8.784265e-01 8.816800e-01 8.849575e-01 8.882350e-01 8.915370e-01 8.948390e-01 8.981655e-01 9.014920e-01 9.048435e-01 9.081950e-01 9.115710e-01 9.149470e-01 9.183485e-01 9.217500e-01 9.251765e-01 9.286030e-01 9.320550e-01 9.355070e-01 9.389845e-01 9.424620e-01 9.459660e-01 9.494700e-01 9.529995e-01 9.565290e-01 9.600845e-01 9.636400e-01 9.672225e-01 9.708050e-01 9.744140e-01 9.780230e-01 9.816585e-01 9.852940e-01 9.889570e-01 9.926200e-01 9.963100e-01 1.000000e+00\n1.295000e+00 1.298750e+00 1.464610e+00 1.659240e+00 1.890670e+00 2.167490e+00 2.500670e+00 2.904300e+00 3.396640e+00 4.001450e+00 4.750000e+00 5.767150e+00 7.070720e+00 8.758190e+00 1.096570e+01 1.388560e+01 1.779290e+01 2.308550e+01 3.034710e+01 4.044480e+01 5.468640e+01 7.507240e+01 1.047120e+02 1.485170e+02 2.143800e+02 3.152120e+02 4.725370e+02 7.229460e+02 1.129950e+03 1.806160e+03 2.955930e+03 4.958860e+03 8.538140e+03 1.510790e+04 2.751070e+04 5.162750e+04 1.000000e+05\n-5 -4 -3 -2 -1 1 2 3 4 5 21";
    file << lhapdfFormalia << std::endl;

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

        file << Cutoff(tpdfs.at(-5).Evaluate(x)) << " " << Cutoff(tpdfs.at(-4).Evaluate(x)) << " " << Cutoff(tpdfs.at(-3).Evaluate(x)) << " " << Cutoff(tpdfs.at(-2).Evaluate(x)) << " " << Cutoff(tpdfs.at(-1).Evaluate(x)) << " " << Cutoff(tpdfs.at(1).Evaluate(x)) << " " << Cutoff(tpdfs.at(2).Evaluate(x)) << " " << Cutoff(tpdfs.at(3).Evaluate(x)) << " " << Cutoff(tpdfs.at(4).Evaluate(x)) << " " << Cutoff(tpdfs.at(5).Evaluate(x)) << " " << Cutoff(tpdfs.at(0).Evaluate(x)) << std::endl;
      }
    }
    
    // close output file
    file << "---" << std::endl;
    file.close();


    //////////////////////////////////////////////////////////
    // Calculate Structure Functions
    //////////////////////////////////////////////////////////

    // Effective EW charges for a space-like process (i.e. DIS)
    const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };


    // Define PDFs functions 
    // Cutoff not used!!
    const auto PDFs        = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); };
    const auto PDFsEvolved = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };


    // Initialise and Build Structure Functions
    const auto F2        = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, pto, as, fBq);
    const auto F2Evolved = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolved, pto, as, fBq);

    std::cout << std::to_string(PDFs(1,0).size()) << std::endl;
    std::cout << std::to_string(PDFsEvolved(1,0).size()) << std::endl;

    // Tabulate Structure Functions
    const apfel::TabulateObject<apfel::Distribution> F2total        {[=] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2totalEvolved {[=] (double const& Q) -> apfel::Distribution{ return F2Evolved.at(0).Evaluate(Q);}, 50, 1, 1000, 3, Thresholds};


    ///////////////////////////////////////////////////////////
    // print data of Structure Functions to file and terminal
    //////////////////////////////////////////////////////////

    // opening output file
    std::ofstream file2;
    file2.open(OutputFileStructureFunctions);

    // print basic information () on the following data
    file2 << "# mu2 values:" << std::endl;

    for (double mu2 : arr_mu2) 
    {
      if (mu2 == arr_mu2[0])
        file2 << std::to_string(mu2);
      else
        file2 << ", " << std::to_string(mu2);
    }
    
    file2 << "\n# num_x_vals\n" << xlha2.size() << "\n# x, APFEL++, LHAPDF, ratio" << std::endl;

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


      const std::function<double(double const&, double const&)> F2GRVManual = [&] (double const& x, double const& mu) -> double
      {
        int nf;

        if (mu < 1.5)       nf = 3;
        else if (mu < 4.5)  nf = 4;
        else if (mu < 100.) nf = 5;
        else                nf = 5;
        
        const std::map<int, double> quarkCharges2 = {{1, 1./9.}, {2, 4./9.}, 
                                                     {3, 1./9.}, {4, 4./9.}, 
                                                     {5, 1./9.}, {6, 4./9.}};

        double result = 0.;

        for (int i = 1; i<=nf; i++)
          result += 2 * dist->xfxQ(i, x, mu) * quarkCharges2.at(i);

        return result;
      };


      // file output:
      // Print down PDF at "mu"
      for (double x : xlha2)
        file2 << x << ", " << F2totalEvolved.EvaluatexQ(x, mu) << ", " << F2total.EvaluatexQ(x, mu) << ", " << F2totalEvolved.EvaluatexQ(x, mu)/F2total.EvaluatexQ(x, mu) << std::endl;
        // manual calculation of F2 for GRV (only LO):
        // file2 << x << ", " << F2totalEvolved.EvaluatexQ(x, mu) << ", " << F2GRVManual(x, mu) << ", " << F2totalEvolved.EvaluatexQ(x, mu)/F2GRVManual(x, mu) << std::endl; //debug
    }

    // terminal output:
    // print perturbative order and name of LHAPDF-Set
    std::cout << "____________________________________________________________" << std::endl;
    std::cout << "____________________________________________________________\n" << std::endl;
    std::cout << std::defaultfloat;
    std::cout << "Used Perturbative Order   : " << std::to_string(pto) << std::endl;
    std::cout << "Used ZBoson Mass          : " << Qref << " GeV" << std::endl;
    std::cout << "Used alphas @ ZBoson Mass : " << asref << "\n" << std::endl;
    std::cout << "Used Charm Quark Mass     : " << mc << " GeV" << std::endl;
    std::cout << "Used Charm Quark Mass     : " << mc << " GeV" << std::endl;
    std::cout << "Used Bottom Quark Mass    : " << mb << " GeV" << std::endl;
    std::cout << "Used Top Quark Mass       : " << mt << " GeV\n" << std::endl;
    std::cout << "Used LHAPDF Set           : " << NameLHAPDFSet << std::endl;
    
    // close output file
    file2.close();

  return 0;
}