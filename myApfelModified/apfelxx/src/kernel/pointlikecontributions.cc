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
// manually changeable functions
    double coeffGeneral(double const& x)
    {
        // this coefficient is defined because the PDFs are also multiplied by it in the LHAPDF format
        return x / alphaQED; 
    };




// fixed alphas and coefficients
    double coeffQCD(double const& alphasAtQ)
    {
        return alphasAtQ / (2. * M_PI);
    };




// stuff used to calculate the expectation value of the squared quark charge
    double ExpE2(int const& particle, int const& nf)
    {
        double result = 0.; 
        for (int i_nf = 1; i_nf<= nf; i_nf++)
            {
                // loop through all active quark flavors and add their squared charges together, following the QCD evolution basis
                result += quarkCharges2.at(i_nf) * (1 + rulesMultiplication.at(particle).at(0)) * rulesMultiplication.at(particle).at(i_nf);
            };   
        return result / (nf * 1.);
    };




// functions used to actually calculate the pointlike contributions
    double k0q(double const& x, int const& particle, int const& nf)
    {
        return nf * ExpE2(particle, nf) * 3 * (x * x + (1 -x) * (1 -x));
    };


    double k1qMS(double const& x, int const& particle, int const& nf)
    {
        double k_x = 4 / 3. * (4 - 9 * x -(1 - 4 * x) * log(x) -(1 - 2 * x) * log(x) * log(x) + 4 * log(1 - x) 
                               + (4 * log(x) - 4 * log(x) * log(1 - x) + 2 * log(x) * log(x) - 4 * log(1 - x) + 2 * log(1 - x) * log(1 - x) - 2 / 3. * M_PI * M_PI + 10) 
                               * (x * x + (1 - x) * (1 - x)));
        return nf * ExpE2(particle, nf) * 3 * k_x;
    };


    double k1qDIS(double const& x, int const& particle, int const& nf)
    {
        double ConvolutionP0qqBgamma = - 4 * CF * (
                                                + 5 * x - 3.5
                                                + Pi2 * (4 * x * x - 4 * x + 2)
                                                - log(x)                        * (16 * x * x -  8 * x + 0.5)
                                                + log(1 - x)                    * (16 * x * x - 18 * x + 2.5)
                                                + log((1 - x) / x) * log(x)     * ( 4 * x * x -  2 * x + 1) 
                                                - log((1 - x) / x) * log(1 - x) * ( 4 * x * x -  4 * x + 2)
                                                - dilog((-1 + x) / x)           * ( 4 * x * x -  2 * x + 1));
        return k1qMS(x, particle, nf) - 3. * nf * ExpE2(particle, nf) / 2. * ConvolutionP0qqBgamma;
    };


    double k1gMS(double const& x, int const& particle, int const& nf)
    {
        return nf * ExpE2(particle, nf) * 3 * 4 / 3. * (-16 + 8 * x + 20 / 3. * x * x + 4 / (3. * x) -(6 + 10 * x) * log(x) - 2 * (1 + x) * log(x) * log(x));
    };


    double k1gDIS(double const& x, int const& particle, int const& nf)
    {
        double ConvolutionP0gqBgamma = 4. / 3. * CF * (
                                                    + 16 * x *  + 2 * x - 20 + 2. / x
                                                    - log(x)                    * 24 * x
                                                    + log(1 - x)                * (6 * x + 4. / x)
                                                    + log((1 - x) / x)          * (3 - 9 * x - 4 * x * x) 
                                                    + log((1 - x) / x) * log(x) * (6 * x + 6) 
                                                    - dilog((-1 + x) / x)       * (6 * x + 6));
        return k1gMS(x, particle, nf) - 3. * nf * ExpE2(particle, nf) / 2. * ConvolutionP0gqBgamma;
    };


    std::map<int, double> PointlikeMap(int const& particle, int const& nf, double const& x)
    {
        if (particle == GLUON)
        {
            return {{0, 0*x*nf*particle},       // k0g
                    {1, k1gDIS(x, particle, nf)}}; // k1gDIS
        }
        else
        {           
            return {{0, k0q(x, particle, nf)},  // k0q
                    {1, k1qDIS(x, particle, nf)}}; // k1qDIS
        }
    };



    std::function<double(double const&)> PointlikeContribution (int    const& particle,
                                                                int    const& ptoPL, 
                                                                int    const& nf,
                                                                double const& alphasAtQ)
    {
        return [&] (double const& x) -> double 
            {
                double result = 0;
                for (int i = 0; i<= ptoPL; i++)
                {
                    result += coeffGeneral(x) * PointlikeMap(particle, nf, x).at(i) * coeffQED * pow(coeffQCD(alphasAtQ), i);
                };
                return result;
            };
    };
}