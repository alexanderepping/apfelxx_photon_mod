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
     * @brief Expectation value of the squared quark charge, eq^2. nf = 4 assumed.
     */
    const double eExp2    = 10. / 9. / 4. ;
    
    /**
     * @brief Expectation value of the quartic quark charge eq^4. nf = 4 assumed.
     */
    const double eExp4    = 34. / 81. / 4. ; 
    
    /**
     * @brief QCD coefficient for the pointlike contribution:
     * alphaQCD(x) / (2*Pi). 
     * @param Alphas: the function returning the strong coupling
     * @param x: wrong, probably should be some energy
     */
    double coeffQCD(std::function<double(double const&)> const& Alphas, double const& x) const {return 0 ;};
    //double coeffQCD(std::function<double(double const&)> const& Alphas, double const& x) const {return Alphas(x) / (2. * M_PI) ;};


    /**
     * @brief Perturbative order of the pointlike contributions.
     * To be manually changed in "apfel/pointlikecontributions.h".
     * @param particle: particle for which pointlike contribution is to be calculated
     * @param PerturbativeOrderPointlike: Perturbative order of the pointlike contributions. To be manually changed in "apfel/pointlikecontributions.h".
     * @param nf: number of active flavours
     * @param Alphas: the function returning the strong coupling
     */
    std::function<double(double const&)> PointlikeContribution (int                                  const& particle,
                                                                int                                  const& ptoPL, 
                                                                int                                  const& nf,
                                                                std::function<double(double const&)> const& Alphas) const;

    
    const std::map<int, std::function<double(double const&, double const&)>> GluonMap;
    const std::map<int, std::function<double(double const&, double const&)>> SigmaMap;
    const std::map<int, std::function<double(double const&, double const&)>> ValenceMap;
    const std::map<int, std::function<double(double const&, double const&)>> NSMap;

    const std::map<int, std::map<int, std::function<double(double const&, double const&)>>> PointlikeMap;

    
}