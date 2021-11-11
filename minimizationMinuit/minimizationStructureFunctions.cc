///////////////////////////////////////
// imports
///////////////////////////////////////
#include "minimizationMinuit.h"

#include "/usr/local/include/minuit-cpp/FunctionMinimum.hh"
#include "/usr/local/include/minuit-cpp/MnUserParameterState.hh"
#include "/usr/local/include/minuit-cpp/MnPrint.hh"
#include "/usr/local/include/minuit-cpp/MnMigrad.hh"
//#include "/usr/local/include/minuit-cpp/"

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
    StructureFunctionsFcn StructureFunctions(experimentalData);

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