//
// APFEL++ 2017
//
// Author: Alexander Epping: a_eppi01@uni-muenster.de
//

#pragma once

#include "apfel/constants.h"

#include <functional>
#include <map>
#include <vector>


namespace apfel
{
// manually changeable variables & functions
    /**
     * @brief Perturbative order of the pointlike contributions.
     * To be manually changed in "apfel/pointlikecontributions.h".
     */
    const int ptoPL = 0;
    
    
    /**
     * @brief general coefficient
     * @param x: x-value
     */
    double coeffGeneral(double const& x);




// fixed alphas and coefficients
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




// stuff used to calculate the expectation value of the squared quark charge
    /**
     * @brief The map enumerators, describing the QCD evolution basis
     */
    enum ParticleEnumerator: int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};


    /**
     * @brief Quark charges squared. To be changed manually.
     *///{1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}
    const std::map<int, double> quarkCharges2 = {{1, 1./9.}, {2, 4./9.}, 
                                           {3, 1./9.}, {4, 4./9.}, 
                                           {5, 1./9.}, {6, 4./9.}};


    /**
     * @brief Multiplication rules of the QCD evolution basis.
     * rulesMultiplication.at(particle).at(0) is the sign of the q; +1 for q+qbar, -1 for q-qbar.
     * The following numbers are the coefficients of the quarks:             d,  u,  s,  c,  b,  t. 
     */                                             
    const std::map<int, std::vector<int>> rulesMultiplication = {{GLUON,   {+1,  1,  1,  1,  1,  1,  1}},
                                                                 {SIGMA,   {+1,  1,  1,  1,  1,  1,  1}},
                                                                 {VALENCE, {-1,  1,  1,  1,  1,  1,  1}},
                                                                 {T3,      {+1, -1,  1,  0,  0,  0,  0}},
                                                                 {V3,      {-1, -1,  1,  0,  0,  0,  0}},
                                                                 {T8,      {+1,  1,  1, -2,  0,  0,  0}},
                                                                 {V8,      {-1,  1,  1, -2,  0,  0,  0}},
                                                                 {T15,     {+1,  1,  1,  1, -3,  0,  0}},
                                                                 {V15,     {-1,  1,  1,  1, -3,  0,  0}},
                                                                 {T24,     {+1,  1,  1,  1,  1, -4,  0}},
                                                                 {V24,     {-1,  1,  1,  1,  1, -4,  0}},
                                                                 {T35,     {+1,  1,  1,  1,  1,  1, -5}},
                                                                 {V35,     {-1,  1,  1,  1,  1,  1, -5}}};
    

    /**
     * @brief calculation of *<e²>
     * @param particle: "particle" for which pointlike contribution is to be calculated. Particle meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double ExpE2(int const& particle, int const& nf);




// functions used to actually calculate the pointlike contributions
    /**
     * - k0(x) = n_f * <e²> * N_c * (x² +(1-x)²)
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01)
     * @brief k0-function, used with all 0th order pointlike contributions.
     * @param x: x-value
     * @param particle: "particle" for which pointlike contribution is to be calculated. Particle meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k0(double const& x, int const& particle, int const& nf);


    /**
     * - k1(x) = n_f * <e²> * N_c * 4/3 * {4 -9x -(1 -4x)log(x) -(1-2x)log²(x) +4log(1-x) +[4log(x) -4log(x)log(1-x) +2log²(x) -4log(1-x) +2log²(1-x) -2/3*pi² +10][x² +(1-x)²]}
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01)
     * - might be wrong (bc of wrong scheme)
     * @brief k1-function, used with all st order quark pointlike contributions.
     * @param x: x-value
     * @param particle: "particle" for which pointlike contribution is to be calculated. Particle meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k1(double const& x, int const& particle, int const& nf);


    /**
     * - k1g(x) = n_f * <e²> * N_c * 4/3 * [-16 +8x +12/3*x² +4/(3x) -(6 +10x)log(x) -2(1+x)log²(x)]
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01)
     * - might be wrong (bc of wrong scheme)
     * @brief k1g-function, used by the 1st order pointlike contribution of the gluon. Taken from
     * @param x: x-value
     * @param particle: "particle" for which pointlike contribution is to be calculated. Particle meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k1g(double const& x, int const& particle, int const& nf);


    /**
     * @brief Calculate the pointlike contribution (w/out any coefficients), e.g. k0q.
     * @param particle: "particle" for which pointlike contribution is to be calculated. Particle meaning in the QCD evolution basis
     * @param nf: number of active flavours
     * @param x: x-value
     */
    std::map<int, double> PointlikeMap(int const& particle, int const& nf, double const& x);


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