


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
const std::string InputFileName = "/home/alexander/Uni/apfelxx_photon_mod/results/Bestandsaufnahme_2022_09_15/data_SAL5HOInitialPDFsErrors.txt";

// Name of the output file for the Structure Functions of the Evolved PDFs
const std::string OutputFileStructureFunctions = "/home/alexander/Uni/apfelxx_photon_mod/results/Bestandsaufnahme_2022_09_15/data_SAL5HOStructureFunctionsErrors.txt";

// Decide whether LO or HO should be used
//#define LO
#define HO



////////////////////////////////////////////////////////////
// defintions
////////////////////////////////////////////////////////////

// array of final scale values for which data should be output (mu^2)
//                 |  0,1,2:ALEPH  |  3,4,5:AMY
//const double arr_mu2[] = {9.9, 20.7, 284, 6.8, 73, 390};   
const double arr_mu2[] = {1.86, 1.9, 2.4, 2.8, 3.7, 3.76, 4.3, 5.0, 5.1, 5.3, 5.9, 6.8, 7.5, 8.9, 9.0, 9.2, 9.9, 10.7, 10.8, 12.0, 14.5, 14.7, 15.3, 16.0, 17.3, 17.5, 17.8, 20.7, 23.0, 23.1, 24.0, 30.0, 45.0, 59.0, 67.2, 73.0, 80.0, 100.0, 135.0, 284.0, 390.0, 780.0};

// xlha is used for console output of the structure function data
std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

// xlha2 is used for file output and plotting of the structure function data
std::vector<double> xlha2{0.0001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110, 0.115, 0.120, 0.125, 0.130, 0.135, 0.140, 0.145, 0.150, 0.155, 0.160, 0.165, 0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.200, 0.205, 0.210, 0.215, 0.220, 0.225, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.300, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355, 0.360, 0.365, 0.370, 0.375, 0.380, 0.385, 0.390, 0.395, 0.400, 0.405, 0.410, 0.415, 0.420, 0.425, 0.430, 0.435, 0.440, 0.445, 0.450, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505, 0.510, 0.515, 0.520, 0.525, 0.530, 0.535, 0.540, 0.545, 0.550, 0.555, 0.560, 0.565, 0.570, 0.575, 0.580, 0.585, 0.590, 0.595, 0.600, 0.605, 0.610, 0.615, 0.620, 0.625, 0.630, 0.635, 0.640, 0.645, 0.650, 0.655, 0.660, 0.665, 0.670, 0.675, 0.680, 0.685, 0.690, 0.695, 0.700, 0.705, 0.710, 0.715, 0.720, 0.725, 0.730, 0.735, 0.740, 0.745, 0.750, 0.755, 0.760, 0.765, 0.770, 0.775, 0.780, 0.785, 0.790, 0.795, 0.800, 0.805, 0.810, 0.815, 0.820, 0.825, 0.830, 0.835, 0.840, 0.845, 0.850, 0.855, 0.860, 0.865, 0.870, 0.875, 0.880, 0.885, 0.890, 0.895, 0.900};//, 0.905, 0.910, 0.915, 0.920, 0.925, 0.930, 0.935, 0.940, 0.945, 0.950, 0.955, 0.960, 0.965, 0.970, 0.975, 0.980, 0.985, 0.990, 0.995}; 

// Parameters used by SAL
  // K_S, B_G_HAD, C_G_HAD, A_Q_HAD, B_Q_HAD, C_Q_HAD, A_Q_PL, B_Q_PL, A_G_HAD
std::vector<double> SALParameters{0.3, -0.57, 3, 0.065, -0.16, 1, 4.45, 1.9, 0.027}; //from Vadim, same as SAL Table1 Zeus-TR (b+1)
// std::vector<double> SALParameters{0.3, -0.64, 3, 0.0645, -0.17, 1, 4.04, 1.22, 0.0173}; //from Vadim, same as SAL Table1 Zeus-TR (b+1)

// struct to store data from the InputFile
struct InputFileDataStruct {   
          std::string FileName;
          std::string NameInputPDF;

          int NumberOfFreeParams;

          std::vector<double> FinalParameters;
          std::vector<std::vector<double>> FinalErrorParametersPlus;
          std::vector<std::vector<double>> FinalErrorParametersMinus;

          double Chi2;
          double Chi2PerDataPoint;
          double DeltaChi2; 
};



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

InputFileDataStruct GetData(std::string InputFileName){

  InputFileDataStruct InputFileData;

  InputFileData.FileName = InputFileName;

  std::fstream InputFile;
  InputFile.open(InputFileData.FileName, std::ios::in);

  if (InputFile.is_open()){
    std::string temp;
    std::string delimiter = ", ";

    // unimportant
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);

    // NameInputPDF
    std::getline(InputFile, InputFileData.NameInputPDF);

    // NumberOfFreeParams
    if (InputFileData.NameInputPDF.substr(0, InputFileData.NameInputPDF.find("SAL")+3) == "INITIALPDFS_SAL")
      InputFileData.NumberOfFreeParams = std::stoi(InputFileData.NameInputPDF.substr(InputFileData.NameInputPDF.find("SAL")+3));
    else
    {
      std::cerr << "This type of Initial PDF is not implemented yet!\n";
      exit(1);
    }

    // unimportant
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);
    std::getline(InputFile, temp);

    //final parameters
    std::getline(InputFile, temp); 
    while(temp.find(delimiter) < std::string::npos)
    {
      InputFileData.FinalParameters.push_back(std::atof(temp.substr(0, temp.find(delimiter)).c_str())); // make the first number to double and attach to vector
      temp = temp.substr(temp.find(delimiter)+2);
    }
    InputFileData.FinalParameters.push_back(std::atof(temp.c_str())); // make the remaining number to double and attach to vector

    // unimportant
    std::getline(InputFile, temp);

    // FinalErrorParametersPlus
    for (int i=0; i<InputFileData.NumberOfFreeParams; i++)
    {
      std::getline(InputFile, temp);
      std::vector<double> tempVector;
      while(temp.find(delimiter) < std::string::npos)
      {
        tempVector.push_back(std::atof(temp.substr(0, temp.find(delimiter)).c_str())); // make the first number to double and attach to vector
        temp = temp.substr(temp.find(delimiter)+2);
      }
      tempVector.push_back(std::atof(temp.c_str())); // make the remaining number to double and attach to vector
      InputFileData.FinalErrorParametersPlus.push_back(tempVector);
    }

    // unimportant
    std::getline(InputFile, temp);

    // FinalErrorParametersMinus
    for (int i=0; i<InputFileData.NumberOfFreeParams; i++)
    {
      std::getline(InputFile, temp);
      std::vector<double> tempVector;
      while(temp.find(delimiter) < std::string::npos)
      {
        tempVector.push_back(std::atof(temp.substr(0, temp.find(delimiter)).c_str())); // make the first number to double and attach to vector
        temp = temp.substr(temp.find(delimiter)+2);
      }
      tempVector.push_back(std::atof(temp.c_str())); // make the remaining number to double and attach to vector
      InputFileData.FinalErrorParametersMinus.push_back(tempVector);
    }

    // unimportant
    std::getline(InputFile, temp);

    // Chi2
    std::getline(InputFile, temp); 
    InputFileData.Chi2 = std::atof(temp.c_str());

    //unimportant
    std::getline(InputFile, temp); 

    // Chi2PerDataPoint
    std::getline(InputFile, temp); 
    InputFileData.Chi2PerDataPoint = std::atof(temp.c_str());

    //unimportant
    std::getline(InputFile, temp); 

    // DeltaChi2
    std::getline(InputFile, temp); 
    InputFileData.DeltaChi2 = std::atof(temp.c_str());
  }

  return InputFileData;

}



////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////
int main()
{     
  std::vector<std::vector<double>> vectorDeltaX(sizeof(arr_mu2)/sizeof(arr_mu2[0]), {0});
  for (int i=0; i<vectorDeltaX.size(); i++)
    vectorDeltaX[i].resize(xlha2.size());

  std::vector<std::vector<double>> vectorX(sizeof(arr_mu2)/sizeof(arr_mu2[0]), {0});
  for (int i=0; i<vectorX.size(); i++)
    vectorX[i].resize(xlha2.size());

  std::vector<std::vector<double>> vectorXSAL(sizeof(arr_mu2)/sizeof(arr_mu2[0]), {0});
  for (int i=0; i<vectorXSAL.size(); i++)
    vectorXSAL[i].resize(xlha2.size());

  std::vector<std::vector<double>> vectorDeltaXTerm(sizeof(arr_mu2)/sizeof(arr_mu2[0]), {0});
  for (int i=0; i<vectorDeltaXTerm.size(); i++)
    vectorDeltaXTerm[i].resize(xlha.size());

  std::vector<std::vector<double>> vectorXTerm(sizeof(arr_mu2)/sizeof(arr_mu2[0]), {0});
  for (int i=0; i<vectorXTerm.size(); i++)
    vectorXTerm[i].resize(xlha.size());

  std::vector<std::vector<double>> vectorXSALTerm(sizeof(arr_mu2)/sizeof(arr_mu2[0]), {0});
  for (int i=0; i<vectorXSALTerm.size(); i++)
    vectorXSALTerm[i].resize(xlha.size());

  InputFileDataStruct InputFileData = GetData(InputFileName);

#ifdef LO
  const int    pto          = 0;
  const double Qref         = 91.188; // nCTEQ15 input
  // const double Qref         = 91.1876;
  const double asref        = 0.1179973; // nCTEQ15 input
  // const double asref        = 0.128;
  const double mc           = 1.3; // nCTEQ15 input
  // const double mc           = 1.5;
  const double mb           = 4.5;
  const double mt           = 174; // nCTEQ15 input
  // const double mt           = 100;
  const double Qin          = 1.3; // nCTEQ15 input
  // const double Qin          = std::sqrt(2.);
  // const double Qin          = 1.295000e+00;
#endif //LO

#ifdef HO
  const int    pto          = 1;
  const double Qref         = 91.188; // nCTEQ15 input
  // const double Qref         = 91.1876;
  const double asref        = 0.1179973; // nCTEQ15 input
  // const double asref        = 0.11087771034313237;
  const double mc           = 1.3; // nCTEQ15 input
  // const double mc           = 1.5;
  const double mb           = 4.5;
  const double mt           = 174; // nCTEQ15 input
  // const double mt           = 100;
  const double Qin          = 1.3; // nCTEQ15 input
  // const double Qin          = std::sqrt(2.);
  // const double Qin          = 1.295000e+00;
#endif //HO



  //////////////////////////////////////////////////////////
  // Things that are the same for all Parameters
  //////////////////////////////////////////////////////////

    // Define Grid
    const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 6e-1, 3}, apfel::SubGrid{50, 8e-1, 5}}};

    // Define Thresholds and Masses
    const std::vector<double> Thresholds = {0, 0, 0, mc, mb};
    const std::vector<double> Masses = Thresholds;

    // Calculate AlphaS(Q)
    apfel::AlphaQCD a{asref, Qref, Masses, Thresholds, pto};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // Effective EW charges for a space-like process (i.e. DIS) (for Structure Functions)
    const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };


  ///////////////////////////////////////////////////////////////////////////
  // Evolving PDFs and Calculating Structure Functions for ErrorParameters
  ///////////////////////////////////////////////////////////////////////////
    for (int i=0; i<InputFileData.NumberOfFreeParams; i++)
    {

      // Define Input PDFs //tobechanged
      const auto InPDFsPlus  = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(InitialPDFs(x, Q, InputFileData.FinalErrorParametersPlus[i])); };
      const auto InPDFsMinus = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(InitialPDFs(x, Q, InputFileData.FinalErrorParametersMinus[i])); };

      // Initialise and Build DglapObjects
      auto EvolvedPDFsPlus  = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFsPlus,  Qin, pto, as);
      auto EvolvedPDFsMinus = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFsMinus, Qin, pto, as);

      // Tabulate PDFs
      //const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};
      const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsPlus{*EvolvedPDFsPlus, 50, 1, 1000, 3, pto, as};
      const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsMinus{*EvolvedPDFsMinus, 50, 1, 1000, 3, pto, as};

      // Define PDFs functions for calculation of Structure Functions
      const auto PDFsEvolvedPlus  = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFsPlus.EvaluateMapxQ(x,Q);  };
      const auto PDFsEvolvedMinus = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFsMinus.EvaluateMapxQ(x,Q); };

      // Initialise and Build Structure Functions
      const auto F2EvolvedPlus  = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolvedPlus,  pto, as, fBq);
      const auto F2EvolvedMinus = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolvedMinus, pto, as, fBq);

      // Tabulate Structure Functions
      const apfel::TabulateObject<apfel::Distribution> F2totalEvolvedPlus  {[=] (double const& Q) -> apfel::Distribution{ return F2EvolvedPlus.at(0).Evaluate(Q);},  50, 1, 1000, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2totalEvolvedMinus {[=] (double const& Q) -> apfel::Distribution{ return F2EvolvedMinus.at(0).Evaluate(Q);}, 50, 1, 1000, 3, Thresholds};

      for (int i_mu=0; i_mu<vectorDeltaX.size(); i_mu++)
        for (int i_x=0; i_x<vectorDeltaX[i_mu].size(); i_x++)
          // calculating deltaX according to nCTEQ15 eq. (2.33)
          vectorDeltaX[i_mu][i_x] += pow((F2totalEvolvedPlus.EvaluatexQ(xlha2[i_x], pow(arr_mu2[i_mu], 0.5)) - F2totalEvolvedMinus.EvaluatexQ(xlha2[i_x], pow(arr_mu2[i_mu], 0.5))), 2);

      for (int i_mu=0; i_mu<vectorDeltaXTerm.size(); i_mu++)
        for (int i_x=0; i_x<vectorDeltaXTerm[i_mu].size(); i_x++)
          // calculating deltaX according to nCTEQ15 eq. (2.33)
          vectorDeltaXTerm[i_mu][i_x] += pow((F2totalEvolvedPlus.EvaluatexQ(xlha[i_x], pow(arr_mu2[i_mu], 0.5)) - F2totalEvolvedMinus.EvaluatexQ(xlha[i_x], pow(arr_mu2[i_mu], 0.5))), 2);
    }

    for (int i_mu=0; i_mu<vectorDeltaX.size(); i_mu++)
      for (int i_x=0; i_x<vectorDeltaX[i_mu].size(); i_x++)
        // calculating deltaX according to nCTEQ15 eq. (2.33)
        vectorDeltaX[i_mu][i_x] = 0.5 * pow(vectorDeltaX[i_mu][i_x], 0.5);

    for (int i_mu=0; i_mu<vectorDeltaXTerm.size(); i_mu++)
      for (int i_x=0; i_x<vectorDeltaXTerm[i_mu].size(); i_x++)
        // calculating deltaX according to nCTEQ15 eq. (2.33)
        vectorDeltaXTerm[i_mu][i_x] = 0.5 * pow(vectorDeltaXTerm[i_mu][i_x], 0.5);
  


  ///////////////////////////////////////////////////////////////////////////
  // Evolving PDFs and Calculating Structure Functions for normal Parameters
  ///////////////////////////////////////////////////////////////////////////
    // Define Input PDFs //tobechanged
    const auto InPDFs    = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(InitialPDFs(x, Q, InputFileData.FinalParameters)); };
    const auto InPDFsSAL = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(InitialPDFs(x, Q, SALParameters)); };

    // Initialise and Build DglapObjects
    auto EvolvedPDFs    = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs,  Qin, pto, as);
    auto EvolvedPDFsSAL = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFsSAL,  Qin, pto, as);

    // Tabulate PDFs
    //const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3, pto, as};
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsSAL{*EvolvedPDFsSAL, 50, 1, 1000, 3, pto, as};

    // Define PDFs functions for calculation of Structure Functions
    const auto PDFsEvolved    = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q);  };
    const auto PDFsEvolvedSAL = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFsSAL.EvaluateMapxQ(x,Q);  };

    // Initialise and Build Structure Functions
    const auto F2Evolved    = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolved,     pto, as, fBq);
    const auto F2EvolvedSAL = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolvedSAL,  pto, as, fBq);

    // Tabulate Structure Functions
    const apfel::TabulateObject<apfel::Distribution> F2totalEvolved    {[=] (double const& Q) -> apfel::Distribution{ return F2Evolved.at(0).Evaluate(Q);},  50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2totalEvolvedSAL {[=] (double const& Q) -> apfel::Distribution{ return F2EvolvedSAL.at(0).Evaluate(Q);},  50, 1, 1000, 3, Thresholds};

    // Save the values at the mu and x values
    for (int i_mu=0; i_mu<vectorX.size(); i_mu++)
      for (int i_x=0; i_x<vectorX[i_mu].size(); i_x++)
        vectorX[i_mu][i_x] = F2totalEvolved.EvaluatexQ(xlha2[i_x], pow(arr_mu2[i_mu], 0.5));

    for (int i_mu=0; i_mu<vectorXTerm.size(); i_mu++)
      for (int i_x=0; i_x<vectorXTerm[i_mu].size(); i_x++)
        vectorXTerm[i_mu][i_x] = F2totalEvolved.EvaluatexQ(xlha[i_x], pow(arr_mu2[i_mu], 0.5));

    // Save the values at the mu and x values for SAL
    for (int i_mu=0; i_mu<vectorXSAL.size(); i_mu++)
      for (int i_x=0; i_x<vectorXSAL[i_mu].size(); i_x++)
        vectorXSAL[i_mu][i_x] = F2totalEvolvedSAL.EvaluatexQ(xlha2[i_x], pow(arr_mu2[i_mu], 0.5));

    for (int i_mu=0; i_mu<vectorXSALTerm.size(); i_mu++)
      for (int i_x=0; i_x<vectorXSALTerm[i_mu].size(); i_x++)
        vectorXSALTerm[i_mu][i_x] = F2totalEvolvedSAL.EvaluatexQ(xlha[i_x], pow(arr_mu2[i_mu], 0.5));


  ///////////////////////////////////////////////////////////
  // print data of Structure Functions to file and terminal
  //////////////////////////////////////////////////////////

    // opening output file
    std::ofstream file;
    file.open(OutputFileStructureFunctions);

    // print basic information on the following data
    file << "# mu2 values:" << std::endl;

    for (double mu2 : arr_mu2) 
    {
      if (mu2 == arr_mu2[0])
        file << std::to_string(mu2);
      else
        file << ", " << std::to_string(mu2);
    }
    
    file << "\n# num_x_vals\n" << xlha2.size() << "\n# x, APFEL++, LHAPDF, ratio" << std::endl;

    for (int i_mu2=0; i_mu2<sizeof(arr_mu2)/sizeof(arr_mu2[0]); i_mu2++)
    {
      double mu = pow(arr_mu2[i_mu2], 0.5);

      // terminal output
      std::cout << std::scientific;
      std::cout << "____________________________________________________________" << std::endl << std::endl;
      std::cout << "Used Mu                   : " << mu <<  " GeV" << std::endl << std::endl;
      std::cout << "x                  SAL values         StructureFunction  DeltaStructureFunction" << std::endl;

      // terminal output
      for (int i_x=0; i_x<xlha.size(); i_x++)
        std::cout << xlha[i_x] << "       " << vectorXSALTerm[i_mu2][i_x] << "       " << vectorXTerm[i_mu2][i_x] << "       " << vectorDeltaXTerm[i_mu2][i_x] << std::endl;

        // file output
      for (int i_x=0; i_x<xlha2.size(); i_x++)
        file      << xlha2[i_x] << ", " << vectorXSAL[i_mu2][i_x] << ", " << vectorX[i_mu2][i_x] << ", " << vectorDeltaX[i_mu2][i_x] << std::endl;
     
    }

    // terminal output
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

    std::cout << "Used InitialPDF           : " << InputFileData.NameInputPDF << std::endl;
    std::cout << "Used DeltaChi2            : " << InputFileData.DeltaChi2 << std::endl;
    std::cout << "Resulting Chi2            : " << InputFileData.Chi2 << std::endl;
    std::cout << "Resulting Chi2/DataPoint  : " << InputFileData.Chi2PerDataPoint << std::endl;
    
    // close output file
    file.close();

  return 0;
}