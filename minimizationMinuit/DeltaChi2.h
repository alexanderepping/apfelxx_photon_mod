
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

#include <boost/math/special_functions/gamma.hpp>
#include <vector>

/**
 * @brief function to calculate Xi_p from nCTEQ15 (A2); inverse, incomplete gamma function of N/2 and p%
 * 
 * @param p: percentage p%, that fit with Chi2 > Xi_p genuinely describes given set of data
 * @param N: number of data points
 * 
 * @return Xi_p
 */
double XiP(int const& p, 
           int const& N);

/**
 * @brief function to calculate Xi90, rescaled, from nCTEQ15 (A3)
 * 
 * @param Nk: number of data points for experiment k
 * @param Chi2k0: minimum Chi2 for experiment k
 * 
 * @return rescaled Xi90
 */
double Xi90Rescale(int    const& Nk,
                   double const& Chi2k0);