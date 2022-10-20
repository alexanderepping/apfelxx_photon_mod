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

int FileOutputMinimization(StructureFunctionsFcn const& StructureFunctions,
                           resultsDataStruct     const& results,
                           bool                  const& PrintErrorPDFs,
                           std::string           const& outputFile = outputFile);


int TermOutputMinimization(StructureFunctionsFcn const& StructureFunctions,
                           resultsDataStruct     const& results,
                           bool                  const& PrintErrorPDFs);
