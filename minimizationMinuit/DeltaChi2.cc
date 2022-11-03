
/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "DeltaChi2.h"
#include "configMinuit.h"

#include <boost/math/special_functions/gamma.hpp>
#include <vector>

double XiP(int const& p, 
           int const& N)
{
    return 2. * boost::math::gamma_p_inv(N/2., p/100.);
}

double Xi90Rescaled(int    const& Nk,
                   double const& Chi2k0)
{
    return XiP(90, Nk) / XiP(50, Nk) * Chi2k0;
}