/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

///////////////////////////////////////
// imports
///////////////////////////////////////

#include "/usr/local/include/minuit-cpp/FunctionMinimum.hh"
#include "/usr/local/include/minuit-cpp/MnUserParameters.hh"
#include "/usr/local/include/minuit-cpp/MnPrint.hh"
#include "/usr/local/include/minuit-cpp/MnMigrad.hh"

#include "minimizationMinuit.h"

#include <string>
#include <map>
#include <vector>
#include <functional>
#include <iostream>


///////////////////////////////////////
// main program
///////////////////////////////////////
int main()
{
    // create FCN function
    StructureFunctionsFcn StructureFunctions(experimentalData, NameLHAPDFSet);

    double chi2 = StructureFunctions(initialParams);
    std::cout << "§ chi2=" << chi2 << std::endl;
    // create initial starting values for parameters with the 
    MinuitCpp::MnUserParameters userParameters;
    for (int i=0; i<initialParams.size(); i++)
        userParameters.Add(ParamsNames[i], initialParams[i], initialErrorParams[i]);

    // create Migrad minimizer
    MinuitCpp::MnMigrad migrad(StructureFunctions, userParameters);

    // minimize
    MinuitCpp::FunctionMinimum min = migrad();

    // output
    std::cout << "§ minimum: " << min << std::endl;

    return 0;
}