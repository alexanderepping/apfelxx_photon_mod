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
#include <fstream>



///////////////////////////////////////
// main program
///////////////////////////////////////
int main()
{

/**
 * Initialization & Minimization
 */
    // create FCN function
    StructureFunctionsFcn StructureFunctions(experimentalData, NameLHAPDFSet);

    // create initial starting values for parameters with the 
    MinuitCpp::MnUserParameters userParameters;

    // check if lower and/or upper bounds for params are available
    const bool lowerBounds = initialParamsLBounds.find(usedInitialPDFs) != initialParamsLBounds.end();
    const bool upperBounds = initialParamsUBounds.find(usedInitialPDFs) != initialParamsUBounds.end();

    // set initial parameter values
    for (int i=0; i<initialParams.at(usedInitialPDFs).size(); i++)
    {
        userParameters.Add(initialParamsNames.at(usedInitialPDFs)[i], initialParams.at(usedInitialPDFs)[i], initialParamsErrors.at(usedInitialPDFs)[i]);

        // if bounds are defined, set them
        if (lowerBounds & lowerBounds)
            userParameters.SetLimits(initialParamsNames.at(usedInitialPDFs)[i], initialParamsLBounds.at(usedInitialPDFs)[i], initialParamsUBounds.at(usedInitialPDFs)[i]);
        else if (lowerBounds)
            userParameters.SetLowerLimit(initialParamsNames.at(usedInitialPDFs)[i], initialParamsLBounds.at(usedInitialPDFs)[i]);
        else if (upperBounds)
            userParameters.SetUpperLimit(initialParamsNames.at(usedInitialPDFs)[i], initialParamsUBounds.at(usedInitialPDFs)[i]); 
    };

    // create Migrad minimizer
    MinuitCpp::MnMigrad migrad(StructureFunctions, userParameters);

    // minimize
    MinuitCpp::FunctionMinimum min = migrad();

    
/**
 * prepare Data Output
 */
    // get the finalParams; ...Params are vectors etc. with not defined number of parameters
    std::vector<double> finalParams;
    for (int i=0; i<initialParams.at(usedInitialPDFs).size(); i++)
    {
        finalParams.push_back(min.UserParameters().Parameters()[i].Value());
    }

    // get chi2 of finalParams
    const double chi2 = StructureFunctions(finalParams);

    // get the finalParametersMap; ...Parameters are the vectors etc. with all 9(/10 after adding AN_g1) possible parameters
    std::map<int, double> finalParametersMap = StructureFunctions.InitialPDFs(0., 0., finalParams, LHAPDF::mkPDF(NameLHAPDFSet), true);

    // make finalParametersMap into vector
    std::vector<double> finalParameters;
    for (int i=0; i<finalParametersMap.size(); i++)
    {
        finalParameters.push_back(finalParametersMap.at(i));
    }

    // add AN_g1 to vector
    finalParameters.push_back(StructureFunctions.MomentumSumRule(finalParameters));

    

/**
 * File Output
 */
    // opening output file
    std::ofstream file;
    file.open(outputFile);

    file << "# Used InitialPDFs:" << std::endl;
    file << initialPDFsNames.at(usedInitialPDFs) << std::endl;

    file << "# Used experimentalData:" << std::endl;
    for (std::string data : IncludedExperimentalData) 
    {
        if (data == IncludedExperimentalData[IncludedExperimentalData.size()-1])
        file << data << std::endl;
        else
        file << data << ", ";
    }

    file << "# finalParameters names:" << std::endl;
    for (std::string data : initialParamsNames.at(INITIALPDFS_9GDUS)) 
        file << data << ", ";
    file << "AN_g1" << std::endl;

    file << "# finalParameters:" << std::endl;
    for (double data : finalParameters) 
    {
        if (data == finalParameters[finalParameters.size()-1])
        file << data << std::endl;
        else
        file << data << ", ";
    }

    file << "# chi2:" << std::endl;
    file << chi2 << std::endl;

    // close output file
    file.close();

    

/**
 * Terminal Output
 */
    std::cout << "$ Used InitialPDFs:" << std::endl;
    std::cout << "§ " << initialPDFsNames.at(usedInitialPDFs) << std::endl << std::endl;

    std::cout << "§ Used experimentalData:" << std::endl;
    for (std::string data : IncludedExperimentalData) 
    {
        if (data == IncludedExperimentalData[IncludedExperimentalData.size()-1])
        std::cout << "§ " << data << std::endl << std::endl;
        else
        std::cout << "§ " << data << ", ";
    }

    std::cout << "§ finalParameters:" << std::endl;
    for (int i=0; i<finalParameters.size(); i++) 
    {

        if (i == finalParameters.size()-1)
        std::cout << "§ AN_g1:\t" << finalParameters[i] << std::endl << std::endl;
        else
        std::cout << "§ " << initialParamsNames.at(INITIALPDFS_9GDUS)[i] << ":\t" << finalParameters[i] << std::endl;
    }

    std::cout << "§ chi2:" << std::endl;
    std::cout << "§ " << chi2 << std::endl << std::endl << std::endl;
    

    // output
    std::cout << "minimum: " << min << std::endl;

    return 0;
}