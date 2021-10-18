//
// APFEL++ 2017
//
// Author: Alexander Epping: a_eppi01@uni-muenster.de
//

#pragma once

#include "apfel/constants.h"

#include <functional>
#include <map>

namespace apfel
{
    /**
     * @brief Perturbative order of the pointlike contributions.
     * To be manually changed in "apfel/pointlikecontributions.h".
     */
    const int ptoPL = 0;
    
    /**
     * @brief Constant version of the structure constant: 1/137.
     */
    const double alphaQED = 1./137.;
    
    /**
     * @brief QED coefficient for the pointlike contribution:
     * alphaQED / (2*Pi).
     */
    const double coeffQED = alphaQED / (2. * M_PI) ;
    
    /**
     * @brief QCD coefficient for the pointlike contribution:
     * alphaQCD(Q) / (2*Pi). 
     * @param alphasAtQ: value for the alphaQCD at specific Q
     */
    double coeffQCD(double const& alphasAtQ);
    
    /**
     * @brief general coefficient
     * @param x: x-value
     */
    double coeffGeneral(double const& x);
    
    /**
     * @brief Expectation value of the squared quark charge, eq^2. nf = 4 assumed.
     */
    const double eExp2    = 10. / 9. / 4. ;
    
    /**
     * @brief Expectation value of the quartic quark charge eq^4. nf = 4 assumed.
     */
    const double eExp4    = 34. / 81. / 4. ; 


    /**
     * @brief Pointlike contribution for specific particle at
     * specific perturbative order.
     * @param particle: particle for which pointlike contribution is to be calculated
     * @param PerturbativeOrderPointlike: Perturbative order of the pointlike contributions. To be manually changed in "apfel/pointlikecontributions.h".
     * @param nf: number of active flavours
     * @param Alphas: the function returning the strong coupling
     */
    std::function<double(double const&)> PointlikeContribution (int    const& particle,
                                                                int    const& ptoPL, 
                                                                int    const& nf,
                                                                double const& alphasAtQ);
}