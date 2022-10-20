
/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include <string>
#include <vector>
#include <map>

struct resultsDataStruct {
    std::vector<double>                 finalParameters;
    std::vector<std::vector<double>>    finalErrorParametersPlus;
    std::vector<std::vector<double>>    finalErrorParametersMinus;
    double                              chi2;
    std::map<std::string, double>       chi2PerExperiment;
    double                              DeltaChi2;
    bool                                IncludeErrorPDFs = false;
};