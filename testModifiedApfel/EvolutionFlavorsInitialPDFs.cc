


////////////////////////////////////////////////////////////
// imports
////////////////////////////////////////////////////////////

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include "apfel/apfelxx.h" 

// used to write to files
#include<fstream>
#include<iostream>



////////////////////////////////////////////////////////////
// definitions, can be changed by user
////////////////////////////////////////////////////////////

// Name of the input file with the parameters etc
// const std::string InputFileName = "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/data_InitialPDFs.txt";
const std::string InputFileName = "/home/alexander/Uni/apfelxx_photon_mod/results/Bestandsaufnahme_2022_06_23/data_SAL3HOInitialPDFs.txt";
// const std::string InputFileName = "/home/alexander/Uni/apfelxx_photon_mod/results/Bestandsaufnahme_2022_06_23/data_SAL5HOInitialPDFs.txt";

// Name of the output file for the Evolved PDFs
const std::string OutputFileEvolvedPDFs = "/home/alexander/Uni/apfelxx_photon_mod/results/Bestandsaufnahme_2022_06_23/data_SAL3HOEvolved.txt";
// const std::string OutputFileEvolvedPDFs = "/home/alexander/Uni/apfelxx_photon_mod/results/Bestandsaufnahme_2022_06_23/data_SAL5HOEvolved.txt";

// Decide whether LO or HO should be used
//#define LO
#define HO



////////////////////////////////////////////////////////////
// defintions
////////////////////////////////////////////////////////////

//----Array of final scale values for which data should be output----
double arr_mu[] = {std::sqrt(2.)}; 

//----Vector of test values of x (xlha used for console output, xlha2 used for file output and plotting)----
std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
std::vector<double> xlha2{0.0001, 0.0002, 0.00030000000000000003, 0.0004, 0.0005, 0.0006000000000000001, 0.0007, 0.0008, 0.0009000000000000001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009000000000000001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6000000000000001, 0.7000000000000001, 0.8, 0.9};

// Parameters used by SAL
  // K_S, B_G_HAD, C_G_HAD, A_Q_HAD, B_Q_HAD, C_Q_HAD, A_Q_PL, B_Q_PL, A_G_HAD
std::vector<double> SALParameters{0.3, -0.57, 3, 0.065, -0.16, 1, 4.45, 1.9, 0.027}; //from Vadim, same as SAL Table1 Zeus-TR (b+1)
// std::vector<double> SALParameters{0.3, -0.64, 3, 0.0645, -0.17, 1, 4.04, 1.22, 0.0173}; //from Vadim, same as SAL Table1 Zeus-TR (b+1)




////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////

double SinglePDF(double const& x, std::vector<double> Parameters, double const& eQ2) {
  return Parameters[0] * pow(x, Parameters[1]) * pow((1. - x), Parameters[2]) 
       + eQ2 * Parameters[3] * x * (x * x + (1. - x) * (1. - x)) / ( 1. - Parameters[4] * log(1. - x));
}

std::map<int, double> InitialPDFs(double const& x, double const& Q, std::vector<double> Parameters){ 
  // K_S, B_G_HAD, C_G_HAD, A_Q_HAD, B_Q_HAD, C_Q_HAD, A_Q_PL, B_Q_PL, A_G_HAD
  const std::map<int, double> quarkCharges2 = {{1, 1./9.}, {2, 4./9.}, 
                                                {3, 1./9.}, {4, 4./9.}, 
                                                {5, 1./9.}, {6, 4./9.}};
  return {{-6, 0.*Q}, // Q just included to avoid errors
          {-5, 0.}, 
          {-4, 0.}, 
          {-3, SinglePDF(x, {Parameters[0]*Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]}, quarkCharges2.at(3))}, 
          {-2, SinglePDF(x, {Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]}, quarkCharges2.at(2))}, 
          {-1, SinglePDF(x, {Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]}, quarkCharges2.at(1))}, 
          { 0, SinglePDF(x, {Parameters[8], Parameters[1], Parameters[2], 0., 0.}, 0.)}, 
          { 1, SinglePDF(x, {Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]}, quarkCharges2.at(1))}, 
          { 2, SinglePDF(x, {Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]}, quarkCharges2.at(2))}, 
          { 3, SinglePDF(x, {Parameters[0]*Parameters[3], Parameters[4], Parameters[5], Parameters[6], Parameters[7]}, quarkCharges2.at(3))}, 
          { 4, 0.}, 
          { 5, 0.}, 
          { 6, 0.}};
}

std::vector<double> GetData(std::string InputFileName){

  std::vector<double> result;

  std::fstream InputFile;
  InputFile.open(InputFileName, std::ios::in);

  if (InputFile.is_open()){
    std::string temp;
    std::string delimiter = ", ";

    // unimportant
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);

    //final parameters
    std::getline(InputFile, temp); 
    while(temp.find(delimiter) < std::string::npos)
    {
      result.push_back(std::atof(temp.substr(0, temp.find(delimiter)).c_str())); // make the first number to double and attach to vector
      temp = temp.substr(temp.find(delimiter)+2);
    }
    result.push_back(std::atof(temp.c_str())); // make the remaining number to double and attach to vector

  };

  return result;

}



////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////
int main()
{ 
  ///////////////////////////////////////
  // opening and preparing output file
  ///////////////////////////////////////
std::cout << "test" << std::endl;

  // opening output file
  std::ofstream file;
  file.open(OutputFile);

  // print basic information on the following data
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

  
  // get final Parameters
  std::vector<double> FinalParameters = GetData(InputFileName);
  

#ifdef LO
  const int    pto          = 0;
  const int    pto          = 1;
  const double Qref         = 91.188; // nCTEQ15 input
  const double asref        = 0.1179973; // nCTEQ15 input
  const double mc           = 1.3;
  const double mb           = 4.5;
  const double mt           = 174;
  const double Qin          = 1.3;
#endif //LO

#ifdef HO
  const int    pto          = 1;
  const double Qref         = 91.188; // nCTEQ15 input
  const double asref        = 0.1179973; // nCTEQ15 input
  const double mc           = 1.3;
  const double mb           = 4.5;
  const double mt           = 174;
  const double Qin          = 1.3;
#endif //HO

 
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
  


  /////////////////////////////////////
  // Calculate alpha_s
  /////////////////////////////////////

  // Construct alpha_s object
  apfel::AlphaQCD a{asref, Qref, Masses, Thresholds, pto};
  //AlphaQCD a{asref, Qref, Thresholds, pto}; // This constructor can also be used if thresholds and masses are equal

  // Tabulate alpha_s over an interval between 0.9 and 1001 GeV on a
  // grid with 100 intervals distributed as log(log(Q/Lambda)) and set
  // interpolation degree to 3.
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};

  //----calculate alphas(Q) using the evolution calculated by apfel----
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };



  /////////////////////////////////////
  // Calculate PDFs
  /////////////////////////////////////

  // Define input PDFs as a lambda function using the LHAPDF
  // set. Notice that APFEL++ takes PDFs in the QCD evolution basis
  // (g, Sigma, V, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35). The
  // function "PhysToQCDEv" rotates PDFs from the physical basis
  // (tbar, bbar, ..., g, ..., b, t) into the physical basis. One can
  // use this function to define and set of functions to be evolved.
  const auto InPDFs = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(InitialPDFs(x, Q, FinalParameters)); };

  // Initialise evolution setting "Qin" as an initial scale. This
  // means that the function "InPDFs" will be called at "Qin". It uses
  // alpha_s defined above but any other function can be used.
  auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs, Qin, pto, as);
  //auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), InPDFs, Qin, pto, as); // This constructor can also be used if thresholds and masses are equal

  // Finally tabulate PDFs. The constructor is the same as that of
  // "as".
  //const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, pto, as};



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
  std::cout << "Used Perturbative Order   : " << std::to_string(pto) << std::endl;
  std::cout << "Used ZBoson Mass          : " << Qref << " GeV" << std::endl;
  std::cout << "Used alphas @ ZBoson Mass : " << asref << std::endl;
  std::cout << "Used Q_0                  : " << std::to_string(Qin) << " GeV\n" << std::endl;
  std::cout << "Used Charm Quark Mass     : " << mc << " GeV" << std::endl;
  std::cout << "Used Bottom Quark Mass    : " << mb << " GeV" << std::endl;
  std::cout << "Used Top Quark Mass       : " << mt << " GeV\n" << std::endl;
  std::cout << "Used LHAPDF Set           : " << NameLHAPDFSet << std::endl;
  
  // close output file
  file.close();

  return 0;
}

