//
// APFEL++ 2017
//
// Author: Alexander Epping: a_eppi01@uni-muenster.de
//

#include "apfel/constants.h"
#include "apfel/pointlikecontributions.h"

#include <cmath>
#include <algorithm>



namespace apfel 
{
    const std::map<int, std::function<double(double const&, double const&)>> GluonMap   = {{0, [] (double const& x, double const& nf) -> double {
                                                                                                return 0*x*nf;}} //k0g
                                                                                          };

    const std::map<int, std::function<double(double const&, double const&)>> SigmaMap   = {{0, [] (double const& x, double const& nf) -> double {
                                                                                                return 3 * nf * eExp2 * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );}} //k0q
                                                                                          };

    const std::map<int, std::function<double(double const&, double const&)>> ValenceMap = {{0, [] (double const& x, double const& nf) -> double {
                                                                                                return 3 * nf * ( eExp4 - eExp2 * eExp2 ) * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );}} //k0ns
                                                                                          };

    const std::map<int, std::function<double(double const&, double const&)>> NSMap      = {{0, [] (double const& x, double const& nf) -> double {
                                                                                                return 3 * nf * ( eExp4 - eExp2 * eExp2 ) * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );}} //k0ns
                                                                                          };
                                                    
    const std::map<int, std::map<int, std::function<double(double const&, double const&)>>> PointlikeMap = {{ 0, GluonMap}, { 1, SigmaMap}, { 2, ValenceMap},
                                                                                                            { 3, NSMap}, { 4, NSMap}, { 5, NSMap}, { 6, NSMap}, { 7, NSMap}, 
                                                                                                            { 8, NSMap}, { 9, NSMap}, {10, NSMap}, {11, NSMap}, {12, NSMap}};


    double coeffQCD(double const& alphasAtQ){ return alphasAtQ / (2. * M_PI); };



    std::function<double(double const&)> PointlikeContribution (int    const& particle,
                                                                int    const& ptoPL, 
                                                                int    const& nf,
                                                                double const& alphasAtQ)
    {
        return [&] (double const& x) -> double 
            {
                double result = 0;
                for (int i = 0; i<= ptoPL; i++){result = result + PointlikeMap.at(particle).at(0)(x, nf) * coeffQED * pow(coeffQCD(alphasAtQ), i);};
                return result;
            };
    };




}