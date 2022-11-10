/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include "/usr/local/include/minuit-cpp/FCNBase.hh"

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"

#include <Eigen/Eigenvalues>
#include <string>
#include <vector>
#include <map>

/**
 * @brief struct, containing all the important results like final Parameters (with or wothout ErrorPDFs), chi2, etc.
 */
struct resultsDataStruct {
    std::vector<double>                 finalParameters;
    std::vector<std::vector<double>>    finalErrorParametersPlus;
    std::vector<std::vector<double>>    finalErrorParametersMinus;
    double                              chi2;
    std::map<std::string, double>       chi2PerExperiment;
    double                              DeltaChi2;
    bool                                IncludeErrorPDFs = false;
};



/** @brief print debug message to terminal */
void DebugString(std::string const& message, bool const& endl=true, bool const& prefix=true, int const& requiredVerbosity=2);

/** @brief print debug vector of doubles to terminal */
void DebugVectorDoubles(std::vector<double> const& vector, bool const& endl=true, bool const& prefix=true, int const& requiredVerbosity=2);

/** @brief print debug vector of strings to terminal */
void DebugVectorString(std::vector<std::string> const& vector, bool const& endl=true, bool const& prefix=true, int const& requiredVerbosity=2);



/**
 * @brief takes the result-data and writes them to a file. Can output with or without Error PDFs.
 * 
 * @param StructureFunctions: StructureFunctionsFcn Object, for which the minimization was done
 * @param results: struct, containing all the important results like final Parameters etc
 * @param PrintErrorPDFs: bool to decide, if ErrorPDFs data should be output or not
 * @param outputFile: string of filename, where to write to. Default is outputFile from configMinuit.h.
 * 
 * @return returns 0 if output succeeded
 */
int FileOutputMinimization(StructureFunctionsFcn const& StructureFunctions,
                           resultsDataStruct     const& results,
                           bool                  const& PrintErrorPDFs,
                           std::string           const& outputFile = outputFile);



/**
 * @brief takes the result-data and writes them to the terminal. Can output with or without Error PDFs.
 * 
 * @param StructureFunctions: StructureFunctionsFcn Object, for which the minimization was done
 * @param results: struct, containing all the important results like final Parameters etc
 * @param PrintErrorPDFs: bool to decide, if ErrorPDFs data should be output or not
 * 
 * @return returns 0 if output succeeded
 */
int TermOutputMinimization(StructureFunctionsFcn const& StructureFunctions,
                           resultsDataStruct     const& results,
                           bool                  const& PrintErrorPDFs);



/**
 * @brief outputs the elements of the Hessian in such a way, that it is easy to just copy them into a file so that they don't have to be calculated every time.
 * note, that the Hessian changes, if the parameters are changed!
 * 
 * @param Hessian: the Hessian Matrix
 */
void TermOutputHessian(Eigen::MatrixXd const& Hessian);



/**
 * @brief returns a precalculated Hessian, which reduces the runtime, because you don't have to calculate it.
 * note, that the Hessian changes, if the parameters are changed and the precalculated Hessian might not be correct anymore.
 */
Eigen::MatrixXd PrecalculatedHessian();