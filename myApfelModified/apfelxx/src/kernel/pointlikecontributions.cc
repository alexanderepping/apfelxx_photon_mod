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
// fixed alphas and coefficients
    double coeffGeneral(double const& x)
    {
        // this coefficient is applied because the PDFs are also multiplied by it in the LHAPDF format
        return x / alphaQED; 
    };


    double coeffQCD(double const& alphasAtQ)
    {
        return alphasAtQ / (2. * M_PI);
    };




// stuff used to calculate the expectation value of the squared quark charge
    double ExpE2(int const& particleComb, int const& nf)
    {
        double result = 0.; 
        for (int i_nf = 1; i_nf<= nf; i_nf++)
            {
                // loop through all active quark flavors and add their squared charges together, following the QCD evolution basis
                result += quarkCharges2.at(i_nf) * (1 + rulesMultiplication.at(particleComb).at(0)) * rulesMultiplication.at(particleComb).at(i_nf);
            };   
        return result / (nf * 2.);
    };




// functions used to actually calculate the pointlike contributions
    double k0q(double const& x, int const& particleComb, int const& nf)
    {
        return 3 * nf * ExpE2(particleComb, nf) * 2 * (x * x + (1 -x) * (1 -x));
    };


    double k1qMS(double const& x, int const& particleComb, int const& nf)
    {
        double k_x = 4 / 3. * (4 - 9 * x -(1 - 4 * x) * log(x) -(1 - 2 * x) * log(x) * log(x) + 4 * log(1 - x) 
                               + (4 * log(x) - 4 * log(x) * log(1 - x) + 2 * log(x) * log(x) - 4 * log(1 - x) + 2 * log(1 - x) * log(1 - x) - 2 / 3. * Pi2 + 10) 
                               * (x * x + (1 - x) * (1 - x)));
        return 3 * nf * ExpE2(particleComb, nf) * k_x;
    };


    double k1qDIS(double const& x, int const& particleComb, int const& nf)
    {
        double ConvolutionP0qqBgamma = - 4 * CF * ( 
                        + 5 * x - 3.5
                        + Pi2 / 6.                      * ( 4 * x * x -  4 * x + 2. )
                        - log(x)                        * (16 * x * x -  8 * x + 0.5)
                        + log(1 - x)                    * (16 * x * x - 18 * x + 2.5)
                        + log((1 - x) / x) * log(x)     * ( 4 * x * x -  2 * x + 1. ) 
                        - log((1 - x) / x) * log(1 - x) * ( 4 * x * x -  4 * x + 2. )
                        - dilog((-1 + x) / x)           * ( 4 * x * x -  2 * x + 1. ));

        return k1qMS(x, particleComb, nf) - 3. * nf * ExpE2(particleComb, nf) / 2. * ConvolutionP0qqBgamma;
        // return k1qMS(x, particleComb, nf) - 3. * nf * ExpE2(particleComb, nf) / 4. * ConvolutionP0qqBgamma;
        // return k1qMS(x, particleComb, nf) - 3. * nf * ExpE2(particleComb, nf) * ConvolutionP0qqBgamma;
        // return - 3. * nf * ExpE2(particleComb, nf) / 2. * ConvolutionP0qqBgamma;
    };


    double k1gMS(double const& x, int const& particleComb, int const& nf)
    {
        return 3 * nf * ExpE2(particleComb, nf) * 4 / 3. * (-16 + 8 * x + 20 / 3. * x * x + 4 / (3. * x) -(6 + 10 * x) * log(x) - 2 * (1 + x) * log(x) * log(x));
    };


    double k1gDIS(double const& x, int const& particleComb, int const& nf)
    {
        double ConvolutionP0gqBgamma = 4. / (3. * x) * CF * (
                        +                     (2 - 20 * x +  2 * x * x + 16 * x * x * x)
                        + Pi2               * (         x +      x * x) 
                        + log(1 - x)        * (4 +  3 * x -  3 * x * x -  4 * x * x * x) 
                        - log(x)            * (     3 * x + 15 * x * x -  4 * x * x * x)
                        - log(x) * log(x)   * (     3 * x +  3 * x * x)
                        - dilog(x)          * (     6 * x +  6 * x * x));

        return k1gMS(x, particleComb, nf) - 3. * nf * ExpE2(particleComb, nf) / 2. * ConvolutionP0gqBgamma;
        // return k1gMS(x, particleComb, nf) - 3. * nf * ExpE2(particleComb, nf) / 4. * ConvolutionP0gqBgamma;
        // return - 3. * nf * ExpE2(particleComb, nf) / 2. * ConvolutionP0gqBgamma;
    };


    std::map<int, double> PointlikeMap(int const& particleComb, int const& nf, double const& x)
    {
        if (particleComb == GLUON)
        {
            return {{0, 0*x*nf*particleComb},       // k0g
                    {1, k1gDIS(x, particleComb, nf)}}; // k1gDIS
        }
        else
        {           
            return {{0, k0q(x, particleComb, nf)},  // k0q
                    {1, k1qDIS(x, particleComb, nf)}}; // k1qDIS
        }
    };



    std::function<double(double const&)> PointlikeContribution (int    const& particleComb,
                                                                int    const& nf,
                                                                int    const& pto, 
                                                                double const& alphasAtQ)
    {
        return [&] (double const& x) -> double 
            {
                double result = 0;
                for (int i = 0; i<= pto; i++)
                {
                    result += coeffGeneral(x) * PointlikeMap(particleComb, nf, x).at(i) * coeffQED * pow(coeffQCD(alphasAtQ), i);
                };
                return result;
            };
    };
}