/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "HelperFunctions.h"
#include "StructureFunctionsFcn.h"
#include "configMinuit.h"

#include <Eigen/Eigenvalues>
#include <vector>
#include <iostream>
#include <fstream>


void DebugString(std::string const& message, bool const& endl, bool const& prefix, int const& requiredVerbosity)
{
    if (DebugVerbosity >= requiredVerbosity)
    {
        if (prefix)
            std::cout << DebugPrefix;
        std::cout << message;
        if (endl)
            std::cout << std::endl;
    }
}

void DebugVectorDoubles(std::vector<double> const& vector, bool const& endl, bool const& prefix, int const& requiredVerbosity)
{
    if (DebugVerbosity >= requiredVerbosity)
    {
        if (prefix)
            std::cout << DebugPrefix;
        for (int i=0; i<vector.size(); i++) 
        {
            if (i != vector.size()-1)
                std::cout << std::to_string(vector[i]) << ", ";
            else
                std::cout << std::to_string(vector[i]);
        }
        if (endl)
            std::cout << std::endl;
    }
}

void DebugVectorString(std::vector<std::string> const& vector, bool const& endl, bool const& prefix, int const& requiredVerbosity)
{
    if (DebugVerbosity >= requiredVerbosity)
    {
        if (prefix)
            std::cout << DebugPrefix;
        for (int i=0; i<vector.size(); i++) 
        {
            if (i != vector.size()-1)
                std::cout << vector[i] << ", ";
            else
                std::cout << vector[i];
        }
        if (endl)
            std::cout << std::endl;
    }
}

int FileOutputMinimization(StructureFunctionsFcn const& StructureFunctions,
                           resultsDataStruct     const& results,
                           bool                  const& PrintErrorPDFs,
                           std::string           const& outputFile)
{
        // opening output file
        std::ofstream file;
        file.open(outputFile);

        std::string temp="";

    if (PrintErrorPDFs)
        file << "# " << nameUsedInitialPDFs << " ERRORS" << std::endl;
    else
        file << "# " << nameUsedInitialPDFs << std::endl;

        file << "## Used InitialPDFs:" << std::endl;
        file << nameUsedInitialPDFs << std::endl;

        file << "## Used experimentalData:" << std::endl;
        for (std::string data : StructureFunctions.IncludedExpData()) 
            temp += data + ", ";
        temp.resize(temp.size()-2);
        file << temp << std::endl;

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

        // std::cout << "## finalParameters:" << std::endl; //debug
        file << "## finalParameters:" << std::endl;
        temp = "";
        for (double data : results.finalParameters) 
            temp += std::to_string(data) + ", ";
        temp.resize(temp.size()-2);
        file << temp << std::endl;

    if (PrintErrorPDFs)
    {
        file << "## finalErrorParametersPlus:" << std::endl;
        for (std::vector<double> vectorData : results.finalErrorParametersPlus) 
        {
            temp = "";
            for (double data : vectorData) 
                temp += std::to_string(data) + ", ";
            temp.resize(temp.size()-2);
            file << temp << std::endl;
        }

        file << "## finalErrorParametersMinus:" << std::endl;
        for (std::vector<double> vectorData : results.finalErrorParametersMinus) 
        {
            temp = "";
            for (double data : vectorData) 
                temp += std::to_string(data) + ", ";
            temp.resize(temp.size()-2);
            file << temp << std::endl;
        }
    }

        file << "## chi2:" << std::endl;
        file << results.chi2 << std::endl;

        file << "## chi2 per experimentalDataSet:" << std::endl;
        temp = "";
        for (std::string data : StructureFunctions.IncludedExpData()) 
            temp += std::to_string(results.chi2PerExperiment.at(data)) + ", ";
        temp.resize(temp.size()-2);
        file << temp << std::endl;

        file << "## chi2/NumberOfDataPoints:" << std::endl;
        file << results.chi2 / StructureFunctions.F2Gamma().size() << std::endl;

        file << "## chi2/NumberOfDataPoints per experimentalDataSet:" << std::endl;
        temp = "";
        for (std::string data : StructureFunctions.IncludedExpData()) 
            temp += std::to_string(results.chi2PerExperiment.at(data)/StructureFunctions.ExperimentalData().at(data).at("F2Gamma").size()) + ", ";
        temp.resize(temp.size()-2);
        file << temp << std::endl;

    if (PrintErrorPDFs)
    {
        file << "## delta chi2:" << std::endl;
        file << results.DeltaChi2 << std::endl;
    }

        // close output file
        file.close();

        return 0;
}


int TermOutputMinimization(StructureFunctionsFcn const& StructureFunctions,
                           resultsDataStruct     const& results,
                           bool                  const& PrintErrorPDFs)
{
        std::string temp="";

        std::cout << std::endl << std::endl << "Used InitialPDFs:" << std::endl;
        std::cout << "\t" << nameUsedInitialPDFs << std::endl;

        std::cout << "Used experimentalData:" << std::endl;
        for (std::string data : StructureFunctions.IncludedExpData()) 
            temp += data + ", ";
        temp.resize(temp.size()-2);
        std::cout << "\t" << temp << std::endl;

            std::cout << "finalParameters:" << std::endl;
        if (INITIALPDFS_9GDUS <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_5GQ) // if the InitialPDFs use Main0
        {
            for (int i=0; i<results.finalParameters.size(); i++) 
            {

                if (i == results.finalParameters.size()-1)
                    std::cout << "\tAN_g1:\t" << results.finalParameters[i] << std::endl << std::endl;
                else
                {
                    if (initialParamsNames.at(INITIALPDFS_9GDUS)[i].size() < 7)
                        std::cout << "\t" << initialParamsNames.at(INITIALPDFS_9GDUS)[i] << ":\t\t" << results.finalParameters[i] << std::endl;
                    else 
                        std::cout << "\t" << initialParamsNames.at(INITIALPDFS_9GDUS)[i] << ":\t" << results.finalParameters[i] << std::endl;
                }
            }
        }
        else if (INITIALPDFS_SAL8 <= usedInitialPDFs && usedInitialPDFs <= INITIALPDFS_SAL3) // if the InitialPDFs use MainSAL
        {
            for (int i=0; i<results.finalParameters.size(); i++) 
            {

                if (i == results.finalParameters.size()-1)
                    std::cout << "\tA_G_HAD:\t" << results.finalParameters[i] << std::endl;
                else
                {
                    if (initialParamsNames.at(INITIALPDFS_SAL8)[i].size() < 7)
                        std::cout << "\t" << initialParamsNames.at(INITIALPDFS_SAL8)[i] << ":\t\t" << results.finalParameters[i] << std::endl;
                    else
                        std::cout << "\t" << initialParamsNames.at(INITIALPDFS_SAL8)[i] << ":\t" << results.finalParameters[i] << std::endl;
                }
            }
        }

    if (PrintErrorPDFs)
    {
        std::cout << "finalErrorParametersPlus:" << std::endl;
        for (std::vector<double> vectorData : results.finalErrorParametersPlus) 
        {
            temp = "";
            for (double data : vectorData) 
                temp += std::to_string(data) + ", ";
            temp.resize(temp.size()-2);
            std::cout << "\t" << temp << std::endl;
        }

        std::cout << "finalErrorParametersMinus:" << std::endl;
        for (std::vector<double> vectorData : results.finalErrorParametersMinus) 
        {
            temp = "";
            for (double data : vectorData) 
                temp += std::to_string(data) + ", ";
            temp.resize(temp.size()-2);
            std::cout << "\t" << temp << std::endl;
        }
    }

        std::cout << "chi2:" << std::endl;
        std::cout << "\t" << results.chi2 << std::endl;

        std::cout << "chi2 per experimentalDataSet:" << std::endl;
        temp = "";
        for (std::string data : StructureFunctions.IncludedExpData()) 
            temp += data + ", ";
        temp.resize(temp.size()-2);
        std::cout << "\t" << temp << std::endl << "\t";
        
        for (std::string data : StructureFunctions.IncludedExpData()) 
        {
            if (data == StructureFunctions.IncludedExpData()[StructureFunctions.IncludedExpData().size()-1])
                std::cout << results.chi2PerExperiment.at(data) << std::endl;
            else
                std::cout << std::to_string(results.chi2PerExperiment.at(data)) << ", ";
        }
        
        std::cout << "chi2/NumberOfDataPoints:" << std::endl;
        std::cout << "\t" << results.chi2 / StructureFunctions.F2Gamma().size() << std::endl;

        std::cout << "chi2/NumberOfDataPoints per experimentalDataSet:" << std::endl << "\t";
        for (std::string data : StructureFunctions.IncludedExpData()) 
        {
            if (data == StructureFunctions.IncludedExpData()[StructureFunctions.IncludedExpData().size()-1])
                std::cout << results.chi2PerExperiment.at(data)/StructureFunctions.ExperimentalData().at(data).at("F2Gamma").size() << std::endl;
            else
                std::cout << std::to_string(results.chi2PerExperiment.at(data)/StructureFunctions.ExperimentalData().at(data).at("F2Gamma").size()) << ", ";
        }

    if (PrintErrorPDFs)
    {
        std::cout << "Delta Chi^2: "    << std::endl << "\t" << results.DeltaChi2 << std::endl;
    }

        std::cout << std::endl;

        return 0;
}

void TermOutputHessian(Eigen::MatrixXd const& Hessian)
{
    for (int i=0; i<Hessian.innerSize(); i++)
        for (int j=0; j<Hessian.innerSize(); j++)
            DebugString("Hessian("+std::to_string(i) + ", "+std::to_string(j) + ") = " + std::to_string(Hessian(i,j)) + ";", true, false, 1);
}

std::vector<double> PrecalculatedFinalParams()
{
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
    return finalParams;
}


Eigen::MatrixXd PrecalculatedHessian()
{
    Eigen::MatrixXd Hessian(NumberOfFreeParams, NumberOfFreeParams);

#ifdef LO
    if (usedInitialPDFs == INITIALPDFS_SAL3)
    {
        Hessian(0, 0) =     12.3034;
        Hessian(0, 1) =     59.0161;
        Hessian(0, 2) =    -52.4119;
        Hessian(1, 0) =     59.0161;
        Hessian(1, 1) =  11570.5000;
        Hessian(1, 2) = -14398.0000;
        Hessian(2, 0) =    -52.4119;
        Hessian(2, 1) = -14398.0000;
        Hessian(2, 2) =  20965.8000;
    }
    else if (usedInitialPDFs == INITIALPDFS_SAL4VADIM)
    {
        Hessian(0, 0) =    12.918952;
        Hessian(0, 1) =    33.834566;
        Hessian(0, 2) =   -17.225163;
        Hessian(0, 3) =     5.593674;
        Hessian(1, 0) =    33.834566;
        Hessian(1, 1) =  2393.713372;
        Hessian(1, 2) = -2160.927266;
        Hessian(1, 3) =  3339.290884;
        Hessian(2, 0) =   -17.225163;
        Hessian(2, 1) = -2160.927266;
        Hessian(2, 2) =  2777.751032;
        Hessian(2, 3) = -4563.768410;
        Hessian(3, 0) =     5.593674;
        Hessian(3, 1) =  3339.290884;
        Hessian(3, 2) = -4563.768410;
        Hessian(3, 3) =  7979.411564;
    }
    else if (usedInitialPDFs == INITIALPDFS_SAL5)
    {
        Hessian(0, 0) =    13.133108;
        Hessian(0, 1) =    30.608777;
        Hessian(0, 2) =   -18.199265;
        Hessian(0, 3) =    -4.647307;
        Hessian(0, 4) =    -1.904935;
        Hessian(1, 0) =    30.608777;
        Hessian(1, 1) =  1651.131172;
        Hessian(1, 2) = -1586.578415;
        Hessian(1, 3) =  2262.220228;
        Hessian(1, 4) =   -99.486074;
        Hessian(2, 0) =   -18.199265;
        Hessian(2, 1) = -1586.578415;
        Hessian(2, 2) =  2315.570915;
        Hessian(2, 3) = -3941.330014;
        Hessian(2, 4) =    82.481197;
        Hessian(3, 0) =    -4.647307;
        Hessian(3, 1) =  2262.220228;
        Hessian(3, 2) = -3941.330014;
        Hessian(3, 3) =  7646.089002;
        Hessian(3, 4) =  -109.864891;
        Hessian(4, 0) =    -1.904935;
        Hessian(4, 1) =   -99.486074;
        Hessian(4, 2) =    82.481197;
        Hessian(4, 3) =  -109.864891;
        Hessian(4, 4) =     9.429489;
    }
#endif //LO

#ifdef HO
    if (usedInitialPDFs == INITIALPDFS_SAL3)
    {
        Hessian(0, 0) = 11.275728;
        Hessian(0, 1) = 127.361011;
        Hessian(0, 2) = -207.095668;
        Hessian(1, 0) = 127.361011;
        Hessian(1, 1) = 9892.509396;
        Hessian(1, 2) = -12871.692718;
        Hessian(2, 0) = -207.095668;
        Hessian(2, 1) = -12871.692718;
        Hessian(2, 2) = 19336.950953;
    }
    else if (usedInitialPDFs == INITIALPDFS_SAL4VADIM)
    {
        Hessian(0, 0) = 18.776399;
        Hessian(0, 1) = 115.308802;
        Hessian(0, 2) = -182.573521;
        Hessian(0, 3) = 239.261653;
        Hessian(1, 0) = 115.308802;
        Hessian(1, 1) = 3453.409756;
        Hessian(1, 2) = -3720.603211;
        Hessian(1, 3) = 4381.288628;
        Hessian(2, 0) = -182.573521;
        Hessian(2, 1) = -3720.603211;
        Hessian(2, 2) = 5007.902371;
        Hessian(2, 3) = -6123.447014;
        Hessian(3, 0) = 239.261653;
        Hessian(3, 1) = 4381.288628;
        Hessian(3, 2) = -6123.447014;
        Hessian(3, 3) = 7737.402750;
    }
    else if (usedInitialPDFs == INITIALPDFS_SAL5)
    {
        Hessian(0, 0) = 21.671858;
        Hessian(0, 1) = 95.200237;
        Hessian(0, 2) = -175.013106;
        Hessian(0, 3) = 260.917683;
        Hessian(0, 4) = -2.996099;
        Hessian(1, 0) = 95.200237;
        Hessian(1, 1) = 2333.970536;
        Hessian(1, 2) = -2718.078311;
        Hessian(1, 3) = 2993.359887;
        Hessian(1, 4) = -83.276462;
        Hessian(2, 0) = -175.013106;
        Hessian(2, 1) = -2718.078311;
        Hessian(2, 2) = 4037.734106;
        Hessian(2, 3) = -5021.838038;
        Hessian(2, 4) = 94.693121;
        Hessian(3, 0) = 260.917683;
        Hessian(3, 1) = 2993.359887;
        Hessian(3, 2) = -5021.838038;
        Hessian(3, 3) = 6910.487481;
        Hessian(3, 4) = -102.556293;
        Hessian(4, 0) = -2.996099;
        Hessian(4, 1) = -83.276462;
        Hessian(4, 2) = 94.693121;
        Hessian(4, 3) = -102.556293;
        Hessian(4, 4) = 3.967599;
    }
#endif //HO
    return Hessian;
}