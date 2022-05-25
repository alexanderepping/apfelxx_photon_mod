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
 * prepare Data Output and Calculation of Error PDFs
 */
    // get the finalParams; ...Params are vectors etc. with not defined number of parameters
    std::vector<double> finalParams;
    for (int i=0; i<initialParams.at(usedInitialPDFs).size(); i++)
    {
        finalParams.push_back(min.UserParameters().Parameters()[i].Value());
    }

    // get chi2 of finalParams
    const double chi2 = StructureFunctions(finalParams);

    // get the finalParametersMap; ...Parameters are the vectors etc. with all possible parameters of the InitialPDFsMain function
    std::map<int, double> finalParametersMap = StructureFunctions.InitialPDFs(0.5, Qin, finalParams, LHAPDF::mkPDF(NameLHAPDFSet), true);

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
 * just small, first output of minimization results (repeated later again, but bad behaviours or sth else can be already seen here)
 */
#ifdef ErrorPDFs
    std::cout << "§ chi2:" << std::endl;
    std::cout << "§ " << chi2 << std::endl;

    std::cout << "§ chi2/NumberOfDataPoints:" << std::endl;
    std::cout << "§ " << chi2 / StructureFunctions.F2Gamma().size() << std::endl << std::endl << std::endl;
    

    // output
    std::cout << "minimum: " << min << std::endl; 
#endif //ErrorPDFs


/**
 * Calculation Error PDFs
 */
#ifdef ErrorPDFs
    const Eigen::MatrixXd Hessian = CalculateHessian(StructureFunctions, finalParams);

    Eigen::EigenSolver<Eigen::MatrixXd> EigenSolverHessian(Hessian);


    std::vector<std::vector<double>> errorParamsPlus;
    std::vector<std::vector<double>> errorParamsMinus;

    for (int k=0; k<finalParams.size(); k++)
    {
        std::vector<double> tempPlus;
        std::vector<double> tempMinus;
        for (int i=0; i<finalParams.size(); i++)
        {
            tempPlus.push_back(finalParams[i] + std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,i)));
            tempMinus.push_back(finalParams[i] - std::sqrt(DeltaChi2 / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,i)));
            // tempPlus.push_back( finalParams[i] + std::sqrt(DeltaChi2 / EigenSolverHessian.eigenvalues()[k].real()) * EigenSolverHessian.eigenvectors()(k,i)).real());
            // tempMinus.push_back(finalParams[i] - std::sqrt(DeltaChi2 / EigenSolverHessian.eigenvalues()[k].real()) * EigenSolverHessian.eigenvectors()(k,i)).real());
        }
        errorParamsPlus.push_back(tempPlus);
        errorParamsMinus.push_back(tempMinus);
    }

    std::map<std::string, std::vector<std::vector<double>>> errorParams = {{"+", errorParamsPlus},
                                                                           {"-", errorParamsMinus}};
#endif //ErrorPDFs


/**
 * prepare Data Output for Error PDFs
 */
#ifdef ErrorPDFs
    // make finalParametersMap into vector
    std::vector<std::vector<double>> finalErrorParametersPlus;
    std::vector<std::vector<double>> finalErrorParametersMinus;

    for (int k=0; k<finalParams.size(); k++)
    {
        // get the finalErrorParametersMap; ...Parameters are the vectors etc. with all possible parameters of the InitialPDFsMain function
        std::map<int, double> finalErrorParametersMapPlus  = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsPlus[k], LHAPDF::mkPDF(NameLHAPDFSet), true);
        std::map<int, double> finalErrorParametersMapMinus = StructureFunctions.InitialPDFs(0.5, Qin, errorParamsMinus[k], LHAPDF::mkPDF(NameLHAPDFSet), true);

        for (int i=0; i<finalErrorParametersMapPlus.size(); i++)
        {
            finalErrorParametersPlus[k].push_back(finalErrorParametersMapPlus.at(i));
            finalErrorParametersMinus[k].push_back(finalErrorParametersMapMinus.at(i));
        }

        if (INITIALPDFS_9GDUS <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_5GQ) // if the InitialPDFs use Main0
            {
                // add AN_g1 to vector
                finalErrorParametersPlus[k].push_back(finalParameters[finalParameters.size()-1]);
                finalErrorParametersMinus[k].push_back(finalParameters[finalParameters.size()-1]);
            }
        else if (INITIALPDFS_SAL8 <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_SAL3) // if the InitialPDFs use MainSAL
            {
                // add A_G_Had to vector
                finalErrorParametersPlus[k].push_back(finalParameters[finalParameters.size()-1]);
                finalErrorParametersMinus[k].push_back(finalParameters[finalParameters.size()-1]);
            }
    }
#endif //ErrorPDFs
    

/**
 * File Output Minimization and ErrorPDFs
 */
    // opening output file
    std::ofstream file;
    file.open(outputFile);

    file << "# " << initialPDFsNames.at(usedInitialPDFs) << std::endl;

    file << "## Used InitialPDFs:" << std::endl;
    file << initialPDFsNames.at(usedInitialPDFs) << std::endl;

    file << "## Used experimentalData:" << std::endl;
    for (std::string data : IncludedExperimentalData) 
    {
        if (data == IncludedExperimentalData[IncludedExperimentalData.size()-1])
            file << data << std::endl;
        else
            file << data << ", ";
    }

    file << "## finalParameters names:" << std::endl;

    if (INITIALPDFS_9GDUS <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_5GQ) // if the InitialPDFs use Main0
    {
        for (std::string data : initialParamsNames.at(INITIALPDFS_9GDUS)) 
            file << data << ", ";
        file << "AN_g1" << std::endl;
    }
    else if (INITIALPDFS_SAL8 <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_SAL3) // if the InitialPDFs use MainSAL
    {
        for (std::string data : initialParamsNames.at(INITIALPDFS_SAL8)) 
            file << data << ", ";
        file << "A_G_HAD" << std::endl;
    }

    file << "## finalParameters:" << std::endl;
    for (double data : finalParameters) 
    {
        if (data == finalParameters[finalParameters.size()-1])
            file << data << std::endl;
        else
            file << data << ", ";
    }

#ifdef ErrorPDFs
    file << "## finalErrorParametersPlus:" << std::endl;
    for (std::vector<double> vectorData : finalErrorParametersPlus) 
    {
        for (double data : vectorData) 
        {
            if (data == vectorData[vectorData.size()-1])
                file << data << std::endl;
            else
                file << data << ", ";
        }
    }

    file << "## finalErrorParametersMinus:" << std::endl;
    for (std::vector<double> vectorData : finalErrorParametersMinus) 
    {
        for (double data : vectorData) 
        {
            if (data == vectorData[vectorData.size()-1])
                file << data << std::endl;
            else
                file << data << ", ";
        }
    }
#endif //ErrorPDFs

    file << "## chi2:" << std::endl;
    file << chi2 << std::endl;

    file << "## chi2/NumberOfDataPoints:" << std::endl;
    file << chi2 / StructureFunctions.F2Gamma().size() << std::endl;

#ifdef ErrorPDFs
    file << "## delta chi2:" << std::endl;
    file << DeltaChi2 << std::endl;
#endif //ErrorPDFs

    // close output file
    file.close();
    


/**
 * Terminal Output Minimization
 */
    std::cout << "§ Used InitialPDFs:" << std::endl;
    std::cout << "§ " << initialPDFsNames.at(usedInitialPDFs) << std::endl << std::endl;

    std::cout << "§ Used experimentalData:" << std::endl;
    for (std::string data : IncludedExperimentalData) 
    {
        if (data == IncludedExperimentalData[IncludedExperimentalData.size()-1])
        std::cout << "§ " << data << std::endl << std::endl;
        else
        std::cout << "§ " << data << ", ";
    }

    if (INITIALPDFS_9GDUS <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_5GQ) // if the InitialPDFs use Main0
    {
        std::cout << "§ finalParameters:" << std::endl;
        for (int i=0; i<finalParameters.size(); i++) 
        {

            if (i == finalParameters.size()-1)
                std::cout << "§ AN_g1:\t" << finalParameters[i] << std::endl << std::endl;
            else
            std::cout << "§ " << initialParamsNames.at(INITIALPDFS_9GDUS)[i] << ":\t" << finalParameters[i] << std::endl;
        }
    }
    else if (INITIALPDFS_SAL8 <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_SAL3) // if the InitialPDFs use MainSAL
    {
        std::cout << "§ finalParameters:" << std::endl;
        for (int i=0; i<finalParameters.size(); i++) 
        {

            if (i == finalParameters.size()-1)
                std::cout << "§ A_G_HAD:\t" << finalParameters[i] << std::endl << std::endl;
            else
            std::cout << "§ " << initialParamsNames.at(INITIALPDFS_SAL8)[i] << ":\t" << finalParameters[i] << std::endl;
        }
    }

    

    std::cout << "§ chi2:" << std::endl;
    std::cout << "§ " << chi2 << std::endl;

    std::cout << "§ chi2/NumberOfDataPoints:" << std::endl;
    std::cout << "§ " << chi2 / StructureFunctions.F2Gamma().size() << std::endl << std::endl << std::endl;
    

    // output
    std::cout << "minimum: " << min << std::endl;

    

/**
 * Terminal Output Error PDFs
 */
#ifdef ErrorPDFs
    std::cout << "Delta Chi^2: "    << std::endl << DeltaChi2 << std::endl << std::endl;
    std::cout << "Hessian Matrix: " << std::endl << Hessian << std::endl << std::endl;
    std::cout << "Eigenvalues: "    << std::endl << EigenSolverHessian.eigenvalues() << std::endl << std::endl;
    std::cout << "Eigenvectors: "   << std::endl << EigenSolverHessian.eigenvectors() << std::endl << std::endl;

    std::cout << "Error Parameters:" << std::endl;
    std::cout << "          \t";
    for (int i=0; i<initialParamsNames.at(usedInitialPDFs).size(); i++)
    {
        std::string name = initialParamsNames.at(usedInitialPDFs)[i];
        name.append(8-name.length(), ' ');
        std::cout << name << "\t ";
    }

    for (int k=0; k<finalParams.size(); k++)
    {
        std::cout << std::endl << "k=" << std::to_string(k) << "\tak+ = \t";
        for (int i=0; i<finalParams.size(); i++)
        {
            std::string val = std::to_string(errorParams.at("+")[k][i]);
            val.insert(val.begin(), 9-val.length(), ' ');
            std::cout << val << "\t ";
        }

        std::cout << std::endl << "   " << "\tak- = \t";
        for (int i=0; i<finalParams.size(); i++)
        {
            std::string val = std::to_string(errorParams.at("-")[k][i]);
            val.insert(val.begin(), 9-val.length(), ' ');
            std::cout << val << "\t ";
        }
    }
    std::cout << std::endl;
#endif //ErrorPDFs


    return 0;
}