//
// APFEL++ 2017
//
// Author: Alexander Epping: a_eppi01@uni-muenster.de
//

#include "apfel/constants.h"
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
    double nfExpE2(int const& particle, int const& nf)
    {
        double result = 0.; 
        for (int i_nf = 1; i_nf<= nf; i_nf++)
            {
                // loop through all active quark flavors and add their squared charges together, following the QCD evolution basis
                result += quarkCharges2.at(i_nf) * (1 + rulesMultiplication.at(particle).at(0)) * rulesMultiplication.at(particle).at(i_nf);
            };
        /*//debug:
        std::string a_Particles[13] = {"GLUON  ", "SIGMA  ", "VALENCE", "T3     ", "V3     ", "T8     ", "V8     ", "T15    ", "V15    ", "T24    ", "V24    ", "T35    ", "V35    "};
        std::cout << "<e²> of " << a_Particles[particle] << result << std::endl;
        //:debug*/
        return result;
    };




// functions used to actually calculate the pointlike contributions
    double k0(double const& x)
    {
        return 3 * ( x * x + ( 1 - x ) * ( 1 - x ) );
    };


    std::map<int, double> PointlikeMap(int const& particle, int const& nf, double const& x)
    {
        if (particle == GLUON)
        {
            return {{0, 0*x*nf*particle}, // k0g
                    {1, 0*x*nf*particle}};// k1g
        }
        else
        {           
            return {{0, k0(x) * nfExpE2(particle,nf)},  // k0q
                    {1, 0 * nf * x * particle       }}; // k1q
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