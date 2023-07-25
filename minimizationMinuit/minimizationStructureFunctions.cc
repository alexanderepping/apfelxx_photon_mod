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
    for (int i=0; i<NumberOfFreeParams; i++)
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

#ifdef CalculateFinalParams
    // create Migrad minimizer
    MinuitCpp::MnMigrad migrad(StructureFunctions, userParameters);

    // minimize
    MinuitCpp::FunctionMinimum min = migrad();
#endif //CalculateFinalParams


    
/**
 * prepare Data Output
 */
    // get the finalParams; ...Params are vectors etc. with not defined number of parameters
    std::vector<double> finalParams;

#ifdef CalculateFinalParams
    for (int i=0; i<NumberOfFreeParams; i++)
    {
        finalParams.push_back(min.UserParameters().Parameters()[i].Value());
    }
#endif //CalculateFinalParams

#ifndef CalculateFinalParams
    finalParams = PrecalculatedFinalParams();
#endif //CalculateFinalParams

    // get chi2 of finalParams and chi2 for every experiment for finalParams
    const std::map<std::string, double> chi2PerExperiment = StructureFunctions.Chi2PerExperiment(finalParams);
    const double chi2 = chi2PerExperiment.at("All");

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



#ifdef CalculateErrorPDFs
/**
 * just small, first output of minimization results (repeated later again, but bad behaviours or sth else can be already seen here)
 */
    TermOutputMinimization(StructureFunctions, results, false);

#ifdef CalculateFinalParams
    std::cout << "minimum: " << min << std::endl; 
#endif //CalculateFinalParams



/**
 * Calculation Hessian
 */
    if (DebugVerbosity > 1)
    {
        DebugString("Changed Apfel-Verbosity from " + std::to_string(apfel::GetVerbosityLevel()), false); //debug
        apfel::SetVerbosityLevel(1);
        DebugString(" to " + std::to_string(apfel::GetVerbosityLevel()), true, false); //debug
    }

    DebugString("before calculating Hessian"); //debug
#ifdef CalculateHessian
    const Eigen::MatrixXd Hessian = CalculateHessianMatrix(StructureFunctions, finalParams);
#endif //CalculateHessian
#ifndef CalculateHessian
    const Eigen::MatrixXd Hessian = PrecalculatedHessian();
#endif //CalculateHessian
    DebugString("before calculating Hessian"); //debug

    TermOutputHessian(Hessian); //debug

    // V_i^(k) is the ith element of the kth eigenvector (V^(k) corresponding to lambda_k)
    // EigenSolverHessian.eigenvectors().col(k) = V^(k)
    // EigenSolverHessian.eigenvectors().col(k)[i] = V_i^(k)
    // EigenSolverHessian.eigenvectors()(i,k) = V_i^(k)
    Eigen::EigenSolver<Eigen::MatrixXd> EigenSolverHessian(Hessian);



/**
 * Calculating the Tolerance, DeltaChi2
 */
    if (DebugVerbosity >= 1 && apfel::GetVerbosityLevel() < 1)
    {
        DebugString("Changed Apfel-Verbosity from " + std::to_string(apfel::GetVerbosityLevel()), false); //debug
        apfel::SetVerbosityLevel(1);
        DebugString(" to " + std::to_string(apfel::GetVerbosityLevel()), true, false); //debug
    }

#ifdef CalculateDeltaChi2
    // calculating the Xi90Rescaled for each experiments
    std::map<std::string, double> Xi90RescaledMap;
    for (std::string DataSet : StructureFunctions.IncludedExpData())
        Xi90RescaledMap[DataSet] = Xi90Rescaled(StructureFunctions.ExperimentalData().at(DataSet).at("Q2Data").size(), chi2PerExperiment.at(DataSet));

    // output the Xi90Rescaled, Xi90 and Xi50
    for (std::string DataSet : StructureFunctions.IncludedExpData())
        DebugString("DataSet: " + DataSet + ", Xi90Rescaled: " + std::to_string(Xi90RescaledMap[DataSet]) + ", Xi90: " + std::to_string(XiP(90,StructureFunctions.ExperimentalData().at(DataSet).at("Q2Data").size())) + ", Xi50: " + std::to_string(XiP(50,StructureFunctions.ExperimentalData().at(DataSet).at("Q2Data").size()))); //debug

    std::vector<std::vector<double>> ziPlusMinus;
    for (int i=0; i<NumberOfFreeParams; i++)
    {
        double ziPlus  = CalculateZiPlusMinus(StructureFunctions, EigenSolverHessian, i, finalParams, Xi90RescaledMap, "plus", 100);
        double ziMinus = CalculateZiPlusMinus(StructureFunctions, EigenSolverHessian, i, finalParams, Xi90RescaledMap, "minus", 100);
        
        ziPlusMinus.push_back({ziPlus, ziMinus});

        DebugString("z_" + std::to_string(i) + "^+ = " + std::to_string(ziPlusMinus[i][0]) + ", z_" + std::to_string(i) + "^- = " + std::to_string(ziPlusMinus[i][1])); //debug
    }
    
    DebugString("", true, false); //debug
    for (int i=0; i<NumberOfFreeParams; i++) //debug
        DebugString("z_" + std::to_string(i) + "^+ = " + std::to_string(ziPlusMinus[i][0]) + ", z_" + std::to_string(i) + "^- = " + std::to_string(ziPlusMinus[i][1])); //debug
    
    double DeltaChi2 = 0.;
    for (std::vector<double> params : ziPlusMinus)
        DeltaChi2 += ((params[0] * params[0]) + (params[1] * params[1])) / (2. * NumberOfFreeParams);

    DebugString("", true, false); //debug
    DebugString("DeltaChi2 = " + std::to_string(DeltaChi2), true, true, 1); //debug
#endif //CalculateDeltaChi2



/**
 * Calculation Error PDFs
 */
    std::vector<std::vector<double>> errorParamsPlus;
    std::vector<std::vector<double>> errorParamsMinus;

    DebugString("before the for loop"); //debug
    for (int k=0; k<NumberOfFreeParams; k++)
    {
        std::vector<double> tempPlus;
        std::vector<double> tempMinus;
        DebugString("before the second for loop"); //debug
        for (int i=0; i<NumberOfFreeParams; i++)
        {
            tempPlus.push_back(finalParams[i] + std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(i,k)));
            tempMinus.push_back(finalParams[i] - std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(i,k)));
        }
        errorParamsPlus.push_back(tempPlus);
        errorParamsMinus.push_back(tempMinus);
    }

    DebugString("before the errorParams"); //debug
    std::map<std::string, std::vector<std::vector<double>>> errorParams = {{"+", errorParamsPlus},
                                                                           {"-", errorParamsMinus}};

    DebugString("", true, false); //debug
    DebugString("errorParamsPlus: "); //debug
    for (int k=0; k<NumberOfFreeParams; k++)
        DebugVectorDoubles(errorParamsPlus[k]); //debug

    DebugString("", true, false); //debug
    DebugString("errorParamsMinus: "); //debug
    for (int k=0; k<NumberOfFreeParams; k++)
        DebugVectorDoubles(errorParamsMinus[k]); //debug




/**
 * prepare Data Output for Error PDFs
 */
    // make finalParametersMap into vector
    std::vector<std::vector<double>> finalErrorParametersPlus;
    std::vector<std::vector<double>> finalErrorParametersMinus;

    DebugString("before for-loop in outputprep"); //debug
    for (int k=0; k<NumberOfFreeParams; k++)
    {
        // get the finalErrorParametersMap; ...Parameters are the vectors etc. with all possible parameters of the InitialPDFsMain function
        std::map<int, double> finalErrorParametersMapPlus  = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsPlus[k], true);
        std::map<int, double> finalErrorParametersMapMinus = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsMinus[k], true);

        DebugString("before 2nd for-loop"); //debug

        std::vector<double> finalErrorParametersPlusTemp;
        std::vector<double> finalErrorParametersMinusTemp;

        // copy data finalErrorParameters from Map into vector
        for (int i=0; i<finalErrorParametersMapPlus.size(); i++)
        {
            DebugString("in 2nd for-loop 01"); //debug
            finalErrorParametersPlusTemp.push_back(finalErrorParametersMapPlus.at(i));
            finalErrorParametersMinusTemp.push_back(finalErrorParametersMapMinus.at(i));
            DebugString("in 2nd for-loop 02"); //debug
        }

        DebugString("after 2nd for-loop"); //debug

        // save finalErrorParameters...Temp in the "big" vector
        finalErrorParametersPlus.push_back(finalErrorParametersPlusTemp);
        finalErrorParametersMinus.push_back(finalErrorParametersMinusTemp);

        DebugString("before if"); //debug
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
#endif //CalculateErrorPDFs
    

/**
 * Output
 */
    FileOutputMinimization(StructureFunctions, results, results.IncludeErrorPDFs, outputFile);
    TermOutputMinimization(StructureFunctions, results, results.IncludeErrorPDFs);

#ifdef CalculateErrorPDFs
    std::cout << "Hessian Matrix: " << std::endl << Hessian << std::endl << std::endl;
    std::cout << "Eigenvalues: "    << std::endl << EigenSolverHessian.eigenvalues() << std::endl << std::endl;
    std::cout << "Eigenvectors: "   << std::endl << EigenSolverHessian.eigenvectors() << std::endl << std::endl;
#endif //CalculateErrorPDFs

    TermOutputHessian(Hessian); //debug

    return 0;
}