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

#include <Eigen/Eigenvalues>
#include <string>
#include <map>
#include <vector>
#include <functional>
#include <cmath> //oder <math.h> //für sqrt
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
    StructureFunctionsFcn StructureFunctions(experimentalData);

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

    // get chi2 for every experiment for finalParams
    std::map<std::string, double> chi2PerExperiment;
    for (std::string DataSet : StructureFunctions.IncludedExpData())
    {
        StructureFunctionsFcn StructureFunctionsTemp(StructureFunctions.ExperimentalData(), {DataSet});
        chi2PerExperiment[DataSet] = StructureFunctionsTemp(finalParams);
    }

    // get the finalParametersMap; ...Parameters are the vectors etc. with all possible parameters of the InitialPDFsMain function
    std::map<int, double> finalParametersMap = StructureFunctions.InitialPDFs(0.5, Qin, finalParams, true);

    // make finalParametersMap into vector
    std::vector<double> finalParameters;
    for (int i=0; i<finalParametersMap.size(); i++)
        finalParameters.push_back(finalParametersMap.at(i));

    if (INITIALPDFS_9GDUS <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_5GQ) // if the InitialPDFs use Main0
        // add AN_g1 to vector
        finalParameters.push_back(StructureFunctions.MomentumSumRule0(finalParameters));
    else if (INITIALPDFS_SAL8 <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_SAL3) // if the InitialPDFs use MainSAL
        finalParameters.push_back(StructureFunctions.MomentumSumRuleSAL(finalParameters, Qin));    


    
/**
 * save results in resultsDataStruct
 */
    resultsDataStruct results;
    results.finalParameters = finalParameters;
    results.chi2 = chi2;
    results.chi2PerExperiment = chi2PerExperiment;



#ifdef ErrorPDFs
/**
 * just small, first output of minimization results (repeated later again, but bad behaviours or sth else can be already seen here)
 */
    TermOutputMinimization(StructureFunctions, results, false);
    std::cout << "minimum: " << min << std::endl; 



/**
 * Calculation Hessian and DeltaChi2
 */
    const Eigen::MatrixXd Hessian = CalculateHessian(StructureFunctions, finalParams);
    // std::cout << "debug: calculated Hessian" << std::endl; //debug

    Eigen::EigenSolver<Eigen::MatrixXd> EigenSolverHessian(Hessian);
    // std::cout << "debug: made the Eigensolver" << std::endl; //debug


/**
 * Calculation Error PDFs
 */
    std::vector<std::vector<double>> errorParamsPlus;
    std::vector<std::vector<double>> errorParamsMinus;

    // std::cout << "debug: before the for loop" << std::endl; //debug
    for (int k=0; k<finalParams.size(); k++)
    {
        std::vector<double> tempPlus;
        std::vector<double> tempMinus;
        // std::cout << "debug: before the second for loop" << std::endl; //debug
        for (int i=0; i<finalParams.size(); i++)
        {
            tempPlus.push_back(finalParams[i] + std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,i)));
            tempMinus.push_back(finalParams[i] - std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,i)));
        }
        errorParamsPlus.push_back(tempPlus);
        errorParamsMinus.push_back(tempMinus);
    }

    // std::cout << "debug: before errorParams" << std::endl; //debug
    std::map<std::string, std::vector<std::vector<double>>> errorParams = {{"+", errorParamsPlus},
                                                                           {"-", errorParamsMinus}};



/**
 * prepare Data Output for Error PDFs
 */
    // make finalParametersMap into vector
    std::vector<std::vector<double>> finalErrorParametersPlus;
    std::vector<std::vector<double>> finalErrorParametersMinus;

    // std::cout << "debug: before for-loop in outputprep" << std::endl; //debug
    for (int k=0; k<finalParams.size(); k++)
    {
        // get the finalErrorParametersMap; ...Parameters are the vectors etc. with all possible parameters of the InitialPDFsMain function
        std::map<int, double> finalErrorParametersMapPlus  = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsPlus[k], true);
        std::map<int, double> finalErrorParametersMapMinus = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsMinus[k], true);

        // std::cout << "debug: before 2nd for-loop" << std::endl; //debug

        std::vector<double> finalErrorParametersPlusTemp;
        std::vector<double> finalErrorParametersMinusTemp;

        // copy data finalErrorParameters from Map into vector
        for (int i=0; i<finalErrorParametersMapPlus.size(); i++)
        {
            // std::cout << "debug: in 2nd for-loop 01" << std::endl; //debug
            finalErrorParametersPlusTemp.push_back(finalErrorParametersMapPlus.at(i));
            finalErrorParametersMinusTemp.push_back(finalErrorParametersMapMinus.at(i));
            // std::cout << "debug: in 2nd for-loop 02" << std::endl; //debug
        }

        // std::cout << "debug: after 2nd for-loop" << std::endl; //debug

        // save finalErrorParameters...Temp in the "big" vector
        finalErrorParametersPlus.push_back(finalErrorParametersPlusTemp);
        finalErrorParametersMinus.push_back(finalErrorParametersMinusTemp);

        // std::cout << "debug: before if" << std::endl; //debug
        if (INITIALPDFS_9GDUS <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_5GQ) // if the InitialPDFs use Main0
            {
                // add AN_g1 to vector
                finalErrorParametersPlus[k].push_back(StructureFunctions.MomentumSumRule0(finalErrorParametersPlus[k]));
                finalErrorParametersMinus[k].push_back(StructureFunctions.MomentumSumRule0(finalErrorParametersMinus[k]));
            }
        else if (INITIALPDFS_SAL8 <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_SAL3) // if the InitialPDFs use MainSAL
            {
                // add A_G_Had to vector
                finalErrorParametersPlus[k].push_back(StructureFunctions.MomentumSumRuleSAL(finalErrorParametersPlus[k], Qin));
                finalErrorParametersMinus[k].push_back(StructureFunctions.MomentumSumRuleSAL(finalErrorParametersMinus[k], Qin));
            }
    }


    
/**
 * save results in resultsDataStruct
 */
    results.finalErrorParametersPlus = finalErrorParametersPlus;
    results.finalErrorParametersMinus = finalErrorParametersMinus;
    results.DeltaChi2 = DeltaChi2;
    results.IncludeErrorPDFs = true;
#endif //ErrorPDFs
    

/**
 * Output
 */
    FileOutputMinimization(StructureFunctions, results, results.IncludeErrorPDFs, outputFile);
    TermOutputMinimization(StructureFunctions, results, results.IncludeErrorPDFs);

#ifdef ErrorPDFs
    std::cout << "Hessian Matrix: " << std::endl << Hessian << std::endl << std::endl;
    std::cout << "Eigenvalues: "    << std::endl << EigenSolverHessian.eigenvalues() << std::endl << std::endl;
    std::cout << "Eigenvectors: "   << std::endl << EigenSolverHessian.eigenvectors() << std::endl << std::endl;
#endif //ErrorPDFs

    // output
    std::cout << "minimum: " << min << std::endl;


    return 0;
}