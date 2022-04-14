//
// APFEL++ 2017
//
// Author: Alexander Epping: a_eppi01@uni-muenster.de
//

#include "apfel/constants.h"
#include "apfel/specialfunctions.h"
#include "apfel/pointlikecontributions.h"

#include <cmath>
#include <algorithm>
#include <iostream>



namespace apfel 
{
// coefficients for the whole contribution
    double coeffGeneral(double const& x)
    {
        // this coefficient is applied because the PDFs are also multiplied by it in the LHAPDF format
        return x / alphaQED; 
    };


    double coeffQCD(double const& alphasAtQ)
    {
        return alphasAtQ / (2. * M_PI);
    };




// Expectation Values of e^2 and e^4
    double ExpE2(int const& nf)
    {
        double result = 0.; 
        for (int i_nf = 1; i_nf<= nf; i_nf++)
            {
                // loop through all active quark flavors and add their squared charges together
                result += quarkCharges2.at(i_nf);
            };   
        return result / (nf);
    };


    double ExpE4(int const& nf)
    {
        double result = 0.; 
        for (int i_nf = 1; i_nf<= nf; i_nf++)
            {
                // loop through all active quark flavors and add their squared squared charges together
                result += quarkCharges2.at(i_nf) * quarkCharges2.at(i_nf);
            };   
        return result / (nf);
    };

    


// coefficients for the photon-parton splitting functions
    double coeffSinglet(int const& nf)
    { 
        return 3. * nf * ExpE2(nf); 
    };


    double coeffNonSinglet(int const& nf)
    { 
        return 3. * nf * (ExpE4(nf) - ExpE2(nf) * ExpE2(nf)); 
    };


    double coeffNonSingletContribution(int const& nf)
    { 
        return 1. / (quarkCharges2.at(nf) - ExpE2(nf)); 
    };




// functions used in the calculation of the pointlike contributions
    double k0q(double const& x)
    {
        return 2. * (x * x + (1 -x) * (1 -x));
    };


    double k1q(double const& x)
    {
        return 4 / 3. * (
                +                                         (20 * x * x -29 * x +14.)
                - Pi2 / 3.                              * ( 4 * x * x - 4 * x + 2.)
                - log(x)                                * (           - 4 * x + 1.)
                + log(1 - x)                            * (                     4.)
                - log(x) * log(x)                       * (           - 2 * x + 1.)
                - log((1 - x) / x)                      * ( 8 * x * x - 8 * x + 4.)
                + log((1 - x) / x) * log((1 - x) / x)   * ( 4 * x * x - 4 * x + 2.));
    };


    double k1g(double const& x)
    {
        return 4 / (3. * x) * (
                + 1. / 3.           * (20 * x * x * x +24 * x * x -48 * x + 4)
                - log(x)            * (                10 * x * x + 6 * x    )
                - log(x) * log(x)   * (                 2 * x * x + 2 * x    ));
    };


    double ConvolutionP0qqBgamma(double const& x)
    {
        return - 4 * CF * ( 
                + 5 * x - 3.5
                + Pi2 / 6.                      * ( 4 * x * x -  4 * x + 2. ) 
                - log(x)                        * (16 * x * x -  8 * x + 0.5)
                + log(1 - x)                    * (16 * x * x - 18 * x + 2.5)
                + log((1 - x) / x) * log(x)     * ( 4 * x * x -  2 * x + 1. ) 
                - log((1 - x) / x) * log(1 - x) * ( 4 * x * x -  4 * x + 2. )
                - dilog((-1 + x) / x)           * ( 4 * x * x -  2 * x + 1. ));
    };


    double ConvolutionP0gqBgamma(double const& x)
    {
        return 4. / (3. * x) * CF * (
                +                     (2 - 20 * x +  2 * x * x + 16 * x * x * x)
                + Pi2               * (         x +      x * x) 
                + log(1 - x)        * (4 +  3 * x -  3 * x * x -  4 * x * x * x) 
                - log(x)            * (     3 * x + 15 * x * x -  4 * x * x * x)
                - log(x) * log(x)   * (     3 * x +  3 * x * x)
                - dilog(x)          * (     6 * x +  6 * x * x));
    };




// photon-parton splitting functions for gluon
    double k0Gluon(double const& x, int const& nf)
    { 
        return 0*x*nf; 
    };


    double k1GluonMS(double const& x, int const& nf)
    { 
        return coeffSinglet(nf) * k1g(x); 
    };


    double k1GluonDIS(double const& x, int const& nf)
    { 
        return k1GluonMS(x, nf) - coeffSinglet(nf) / 2. * ConvolutionP0gqBgamma(x); 
    };




// photon-parton splitting functions for quark singlet
    double k0Singlet(double const& x, int const& nf)
    { 
        return coeffSinglet(nf) * k0q(x); 
    };


    double k1SingletMS(double const& x, int const& nf)
    { 
        return coeffSinglet(nf) * k1q(x); 
    };


    double k1SingletDIS(double const& x, int const& nf)
    { 
        return k1SingletMS(x, nf) - coeffSinglet(nf) / 2. * ConvolutionP0qqBgamma(x); 
    };




// photon-parton splitting functions for quark non-singlet
    double k0NonSinglet(double const& x, int const& nf)
    { 
        return coeffNonSinglet(nf) * k0q(x); 
    };


    double k1NonSingletMS(double const& x, int const& nf)
    { 
        return coeffNonSinglet(nf) * k1q(x); 
    };


    double k1NonSingletDIS(double const& x, int const& nf)
    { 
        return k1NonSingletMS(x, nf) - coeffNonSinglet(nf) / 2. * ConvolutionP0qqBgamma(x); 
    };




// calculation of the pointlike contribution
    std::map<int, double> PointlikeMap(int const& contribution, int const& nf, double const& x)
    {
        if (contribution == GLUONS)
            return {{0, k0Gluon(x, nf)},
                    {1, k1GluonDIS(x, nf)}};
        
        else if (contribution == SINGLET)
            return {{0, k0Singlet(x, nf)},
                    {1, k1SingletDIS(x, nf)}};

        else if (contribution == NONSINGLET2)
            return {{0, coeffNonSingletContribution(2) * k0NonSinglet(x, 2)},
                    {1, coeffNonSingletContribution(2) * k1NonSingletDIS(x, 2)}};
        
        else if (contribution == NONSINGLET3)
            return {{0, coeffNonSingletContribution(3) * k0NonSinglet(x, 3)},
                    {1, coeffNonSingletContribution(3) * k1NonSingletDIS(x, 3)}};
        
        else if (contribution == NONSINGLET4)
            return {{0, coeffNonSingletContribution(4) * k0NonSinglet(x, 4)},
                    {1, coeffNonSingletContribution(4) * k1NonSingletDIS(x, 4)}};
        
        else if (contribution == NONSINGLET5)
            return {{0, coeffNonSingletContribution(5) * k0NonSinglet(x, 5)},
                    {1, coeffNonSingletContribution(5) * k1NonSingletDIS(x, 5)}};
        
        else if (contribution == NONSINGLET6)
            return {{0, coeffNonSingletContribution(6) * k0NonSinglet(x, 6)},
                    {1, coeffNonSingletContribution(6) * k1NonSingletDIS(x, 6)}};
        
        else 
            return {{0, 0},
                    {1, 0}};
        
    };


    std::function<double(double const&)> PointlikeContribution (int    const& particleComb,
                                                                int    const& nf,
                                                                int    const& pto, 
                                                                double const& alphasAtQ)
    {
        return [&] (double const& x) -> double 
            {
                int particleCombTemp = particleComb;
                int testInt = 0;

                // test, if quarks, that aren't produced at that nf, would contribute to that particleComb
                // can only happen for the Ts
                for(int j = nf+1; j <=6; j++)
                    testInt += std::abs(contributingKs.at(particleComb)[j]);
                
                // if quarks, that aren't produced at that nf, would contribute, 
                // the pointlike contribution is the same as the one of the Singlet
                if (testInt > 0)
                    particleCombTemp = SIGMA;

                double result = 0;

                // loop through all perturbation orders <= to the current 
                for (int i = 0; i<= pto; i++)
                {
                    for (int contribution = ContributionsEnumerator::GLUONS; contribution <= ContributionsEnumerator::NONSINGLET6; contribution++)
                    {
                        double tempResult = contributingKs.at(particleCombTemp)[contribution] * PointlikeMap(contribution, nf, x).at(i);
                        result += coeffGeneral(x) * coeffQED * pow(coeffQCD(alphasAtQ), i) * tempResult;
                    };
                };

                return result;
            };
    };
}