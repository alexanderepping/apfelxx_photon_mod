/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "OutputFunctions.h"
#include "StructureFunctionsFcn.h"
#include "ResultsFunctions.h"
#include "configMinuit.h"

#include <vector>
#include <iostream>
#include <fstream>


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
