/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include "/usr/local/include/minuit-cpp/FCNBase.hh"

#include "StructureFunctionsFcn.h"
#include "ResultsFunctions.h"
#include "configMinuit.h"

#include <vector>

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
