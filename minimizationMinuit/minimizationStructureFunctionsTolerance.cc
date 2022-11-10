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
    std::cout << "get verbosity: " << std::to_string(apfel::GetVerbosityLevel()) << std::endl; //debug
    apfel::SetVerbosityLevel(1);
    std::cout << "get verbosity: " << std::to_string(apfel::GetVerbosityLevel()) << std::endl; //debug

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

    
/**
 * prepare Data Output
 */
    // get the finalParams; ...Params are vectors etc. with not defined number of parameters
    std::vector<double> finalParams;

#ifdef LO
    if (usedInitialPDFs == INITIALPDFS_SAL3)
        finalParams = {-0.384581, 0.439775, 0.256241};
    else if (usedInitialPDFs == INITIALPDFS_SAL4VADIM)
        finalParams = {-0.224452, 0.436453, 0.570871, 0.320468};
    else if (usedInitialPDFs == INITIALPDFS_SAL5)
        finalParams = {-0.202210, 0.524953, 0.665115, 0.348421, 0.648468};
#endif //LO
#ifdef HO
    if (usedInitialPDFs == INITIALPDFS_SAL3)
        finalParams = {-0.182473, 0.483074, 0.285813};
    else if (usedInitialPDFs == INITIALPDFS_SAL4VADIM)
        finalParams = {-0.215028, 0.491192, 0.505920, 0.247174};
    else if (usedInitialPDFs == INITIALPDFS_SAL5)
        finalParams = {-0.272855, 0.591875, 0.610880, 0.297112, 1.116165};
#endif //HO


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



#ifdef ErrorPDFs
/**
 * just small, first output of minimization results (repeated later again, but bad behaviours or sth else can be already seen here)
 */
    TermOutputMinimization(StructureFunctions, results, false);



/**
 * Calculation Hessian
 */
    // const std::map<std::string, Eigen::MatrixXd> Hessian = CalculateHessianMap(StructureFunctions, finalParams);
    // const Eigen::MatrixXd HessianAll = Hessian.at("All");
    //std::cout << "debug: calculated Hessian" << std::endl; //debug
    /*
    for (int i=0; i<finalParams.size(); i++)
        for (int j=0; j<finalParams.size(); j++)
            std::cout << "HessianAll(" << std::to_string(i) << ", " << std::to_string(j) << ") = " << std::to_string(HessianAll(i,j)) << ";" << std::endl; //debug
    */
#ifdef LO
    Eigen::MatrixXd HessianAll(finalParams.size(), finalParams.size());
    if (usedInitialPDFs == INITIALPDFS_SAL3)
    {
        HessianAll(0, 0) = 12.3034;
        HessianAll(0, 1) = 59.0161;
        HessianAll(0, 2) = -52.4119;
        HessianAll(1, 0) = 59.0161;
        HessianAll(1, 1) = 11570.5;
        HessianAll(1, 2) = -14398.0;
        HessianAll(2, 0) = -52.4119;
        HessianAll(2, 1) = -14398.0;
        HessianAll(2, 2) = 20965.8;
    }
    else if (usedInitialPDFs == INITIALPDFS_SAL4VADIM)
    {
        HessianAll(0, 0) = 12.918952;
        HessianAll(0, 1) = 33.834566;
        HessianAll(0, 2) = -17.225163;
        HessianAll(0, 3) = 5.593674;
        HessianAll(1, 0) = 33.834566;
        HessianAll(1, 1) = 2393.713372;
        HessianAll(1, 2) = -2160.927266;
        HessianAll(1, 3) = 3339.290884;
        HessianAll(2, 0) = -17.225163;
        HessianAll(2, 1) = -2160.927266;
        HessianAll(2, 2) = 2777.751032;
        HessianAll(2, 3) = -4563.768410;
        HessianAll(3, 0) = 5.593674;
        HessianAll(3, 1) = 3339.290884;
        HessianAll(3, 2) = -4563.768410;
        HessianAll(3, 3) = 7979.411564;
    }
    else if (usedInitialPDFs == INITIALPDFS_SAL5)
    {

    }
#endif //LO

    // V_i^(k) is the ith element of the kth eigenvector (V^(k) corresponding to lambda_k)
    // EigenSolverHessian.eigenvectors().col(k) = V^(k)
    // EigenSolverHessian.eigenvectors().col(k)[i] = V_i^(k)
    // EigenSolverHessian.eigenvectors()(i,k) = V_i^(k)
    Eigen::EigenSolver<Eigen::MatrixXd> EigenSolverHessian(HessianAll);

    //std::cout << "debug: made the Eigensolver" << std::endl; //debug

    //debug begin
    int trash = FindDeltaChi2Limit(finalParams, EigenSolverHessian); //debug
    //debug end



/**
 * Calculating the Tolerance, DeltaChi2
 */
#ifdef CalculateDeltaChi2
    // calculating the Xi90Rescaled for all experiments
    std::map<std::string, double> Xi90RescaledMap;
    for (std::string DataSet : StructureFunctions.IncludedExpData())
        Xi90RescaledMap[DataSet] = Xi90Rescaled(StructureFunctions.ExperimentalData().at(DataSet).at("Q2Data").size(), chi2PerExperiment.at(DataSet));

    //std::map<int, std::vector<std::vector<double>>> zikPlusMinusMap;
    std::map<int, std::vector<double>> zikPlusMap;
    std::map<int, std::vector<double>> zikMinusMap;
    std::vector<std::vector<double>> ziPlusMinus;

    for (int i=0; i<finalParams.size(); i++)
    {
        /* 
        // calculation using CalculateZikPlusMinus
        zikPlusMap[i] = CalculateZikPlusMinus(StructureFunctions, EigenSolverHessian, i, finalParams, Xi90RescaledMap, "plus", 200.0, 1);
        zikMinusMap[i] = CalculateZikPlusMinus(StructureFunctions, EigenSolverHessian, i, finalParams, Xi90RescaledMap, "minus", 2.119172, 0.1);

        // calculate, what is given in nCTEQ15, (A5)
        double ziPlus = zikPlusMap.at(i)[0];
        for (double zi : zikPlusMap.at(i))
        {
            std::cout << "ziPlus = " << std::to_string(ziPlus) << ", zi = " << std::to_string(zi) << std::endl; //debug
            if (zi < ziPlus)
                ziPlus = zi;
        }

        // calculate, what is given in nCTEQ15, (A5)
        double ziMinus = zikMinusMap.at(i)[0];
        for (double zi : zikMinusMap.at(i))
        {
            std::cout << "ziMinus = " << std::to_string(ziMinus) << ", zi = " << std::to_string(zi) << std::endl; //debug
            if (zi > ziMinus)
                ziMinus = zi;
        }
        */

        double ziPlus  = CalculateZiPlusMinus(StructureFunctions, EigenSolverHessian, i, finalParams, Xi90RescaledMap, "plus", 100);
        double ziMinus = CalculateZiPlusMinus(StructureFunctions, EigenSolverHessian, i, finalParams, Xi90RescaledMap, "minus", 100);

        
        ziPlusMinus.push_back({ziPlus, ziMinus});
//debug
        std::cout << "§§ z_" << std::to_string(i) << "^+ = " << std::to_string(ziPlusMinus[i][0]) << ", z_" << std::to_string(i) << "^- = " << std::to_string(ziPlusMinus[i][1]) << std::endl;

        /*
        for (std::string DataSet : StructureFunctions.IncludedExpData())
            std::cout << DataSet << " ";

        std::cout << std::endl << "z_i^(k+): ";
        for (double params : zikPlusMap.at(i))
            std::cout << std::to_string(params) << "\t";

        std::cout << std::endl << "z_i^(k-): ";
        for (double params : zikMinusMap.at(i))
            std::cout << std::to_string(params) << "\t";
        
        std::cout << std::endl << std::endl;
        */
//debug end
    }
    
    std::cout << std::endl; //debug
    for (int i=0; i<finalParams.size(); i++) //debug
        std::cout << "§§ z_" << std::to_string(i) << "^+ = " << std::to_string(ziPlusMinus[i][0]) << ", z_" << std::to_string(i) << "^- = " << std::to_string(ziPlusMinus[i][1]) << std::endl; //debug
    
    double DeltaChi2 = 0.;
    for (std::vector<double> params : ziPlusMinus)
        DeltaChi2 += ((params[0] * params[0]) + (params[1] * params[1])) / (2. * finalParams.size());

    std::cout << std::endl << "DeltaChi2 = " << std::to_string(DeltaChi2) << std::endl;
#endif //CalculateDeltaChi2


/**
 * Calculation Error PDFs
 */
    std::cout << std::endl << "DeltaChi2 = " << std::to_string(DeltaChi2) << std::endl; //debug
    std::vector<std::vector<double>> errorParamsPlus;
    std::vector<std::vector<double>> errorParamsMinus;

    std::cout << "debug: before the for loop" << std::endl; //debug
    for (int k=0; k<finalParams.size(); k++)
    {
        std::vector<double> tempPlus;
        std::vector<double> tempMinus;
        std::cout << "debug: before the second for loop" << std::endl; //debug
        for (int i=0; i<finalParams.size(); i++)
        {
            // here was an error with the indizes of the eigenvectors
            tempPlus.push_back(finalParams[i] + std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(i,k)));
            tempMinus.push_back(finalParams[i] - std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(i,k)));
        }
        errorParamsPlus.push_back(tempPlus);
        errorParamsMinus.push_back(tempMinus);
    }

    std::cout << "debug: before errorParams" << std::endl; //debug
    std::map<std::string, std::vector<std::vector<double>>> errorParams = {{"+", errorParamsPlus},
                                                                           {"-", errorParamsMinus}};

    std::cout << "\nerrorParamsPlus: " << std::endl;; //debug
    for (int k=0; k<finalParams.size(); k++)
    {
    std::cout << "§ "; //debug
        for (int i=0; i<finalParams.size(); i++)
            std::cout << std::to_string(errorParamsPlus[k][i]) << ", ";
    std::cout << std::endl; //debug
    }

    std::cout << "\nerrorParamsMinus: " << std::endl;; //debug
    for (int k=0; k<finalParams.size(); k++)
    {
    std::cout << "§ "; //debug
        for (int i=0; i<finalParams.size(); i++)
            std::cout << std::to_string(errorParamsMinus[k][i]) << ", ";
    std::cout << std::endl; //debug
    }




/**
 * prepare Data Output for Error PDFs
 */
    // make finalParametersMap into vector
    std::vector<std::vector<double>> finalErrorParametersPlus;
    std::vector<std::vector<double>> finalErrorParametersMinus;

    std::cout << "debug: before for-loop in outputprep" << std::endl; //debug
    for (int k=0; k<finalParams.size(); k++)
    {
        // get the finalErrorParametersMap; ...Parameters are the vectors etc. with all possible parameters of the InitialPDFsMain function
        std::map<int, double> finalErrorParametersMapPlus  = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsPlus[k], true);
        std::map<int, double> finalErrorParametersMapMinus = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsMinus[k], true);

        std::cout << "debug: before 2nd for-loop" << std::endl; //debug

        std::vector<double> finalErrorParametersPlusTemp;
        std::vector<double> finalErrorParametersMinusTemp;

        // copy data finalErrorParameters from Map into vector
        for (int i=0; i<finalErrorParametersMapPlus.size(); i++)
        {
            std::cout << "debug: in 2nd for-loop 01" << std::endl; //debug
            finalErrorParametersPlusTemp.push_back(finalErrorParametersMapPlus.at(i));
            finalErrorParametersMinusTemp.push_back(finalErrorParametersMapMinus.at(i));
            std::cout << "debug: in 2nd for-loop 02" << std::endl; //debug
        }

        std::cout << "debug: after 2nd for-loop" << std::endl; //debug

        // save finalErrorParameters...Temp in the "big" vector
        finalErrorParametersPlus.push_back(finalErrorParametersPlusTemp);
        finalErrorParametersMinus.push_back(finalErrorParametersMinusTemp);

        std::cout << "debug: before if" << std::endl; //debug
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
    std::cout << "Hessian Matrix: " << std::endl << HessianAll << std::endl << std::endl;
    std::cout << "Eigenvalues: "    << std::endl << EigenSolverHessian.eigenvalues() << std::endl << std::endl;
    std::cout << "Eigenvectors: "   << std::endl << EigenSolverHessian.eigenvectors() << std::endl << std::endl;
#endif //ErrorPDFs

    return 0;
}