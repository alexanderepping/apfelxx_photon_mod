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
    const int ptoPL = 1;
    

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
    enum ParticleCombinationEnumerator: int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};


    /**
     * @brief Quark charges squared. To be changed manually.
     *///{1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}
    const std::map<int, double> quarkCharges2 = {{1, 1./9.}, {2, 4./9.}, 
                                           {3, 1./9.}, {4, 4./9.}, 
                                           {5, 1./9.}, {6, 4./9.}};


    /**
     * @brief Multiplication rules of the QCD evolution basis.
     * rulesMultiplication.at(particleComb).at(0) shows whether q+qbar or q-qbar is used; +1 for q+qbar, -1 for q-qbar.
     * The following numbers are the coefficients of the quarks:                 d,  u,  s,  c,  b,  t. 
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
     * @brief calculation of <e²>.
     * <e²> becomes zero, if the particleCombination consists of only q-qbar.
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double ExpE2(int const& particleComb, int const& nf);




// functions used to actually calculate the pointlike contributions
    /**
     * - k0q(x) = n_f * <e²> * N_c * 2 * (x² +(1-x)²)
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     * - also see my notes (15.2.3 fourth attempt correction) for calculations
     * @brief k0q-function, used with all 0th order pointlike contributions.
     * @param x: x-value
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k0q(double const& x, int const& particleComb, int const& nf);


    /**
     * - k1q_MS(x) = n_f * <e²> * N_c * 4/3 * {4 -9x -(1 -4x)log(x) -(1-2x)log²(x) +4log(1-x) +[4log(x) -4log(x)log(1-x) +2log²(x) -4log(1-x) +2log²(1-x) -2/3*pi² +10][x² +(1-x)²]}
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01)
     * @brief MSbar part of the k1q-function, used with all 1st order quark pointlike contributions.
     * @param x: x-value
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k1qMS(double const& x, int const& particleComb, int const& nf);


    /**
     * - k1q_DIS(x) = k1_MS(x) - 3 * n_f / 2 * <e²> * Convolution(P^0_{qq},B_gamma)(x)
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. (3.1)
     * - calculated using Calculations_etc/CalculationDISAdditionalPartsForK_Pqq.nb and Calculations_etc/CalculationDISPqqBgamma.txt
     * @brief k1q-function, used with all 1st order quark pointlike contributions, including the additional DIS-part
     * @param x: x-value
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k1qDIS(double const& x, int const& particleComb, int const& nf);


    /**
     * - k1g_MS(x) = n_f * <e²> * N_c * 4/3 * [-16 +8x +20/3*x² +4/(3x) -(6 +10x)log(x) -2(1+x)log²(x)]
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01)
     *      - it was wrong in Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01)
     * @brief MSbar part of the k1g-function, used by the 1st order pointlike contribution of the gluon.
     * @param x: x-value
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k1gMS(double const& x, int const& particleComb, int const& nf);


    /**
     * - k1g_DIS(x) = k1g_MS(x) - 3 * n_f / 2 * <e²> * Convolution(P^0_{gq},B_gamma)(x)
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. (3.1)
     * - calculated using Calculations_etc/CalculationDISAdditionalPartsForK_Pgq.nb
     * @brief k1g-function, used by the 1st order pointlike contribution of the gluon, including the additional DIS-part
     * @param x: x-value
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     */
    double k1gDIS(double const& x, int const& particleComb, int const& nf);


    /**
     * @brief Calculate the pointlike contribution (w/out any coefficients), e.g. k0q.
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     * @param x: x-value
     */
    std::map<int, double> PointlikeMap(int const& particleComb, int const& nf, double const& x);


    /**
     * @brief Pointlike contribution for specific particleCombination at
     * specific perturbative order.
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param PerturbativeOrderPointlike: Perturbative order of the pointlike contributions. To be manually changed in "apfel/pointlikecontributions.h".
     * @param nf: number of active flavours
     * @param Alphas: the function returning the strong coupling
     */
    std::function<double(double const&)> PointlikeContribution (int    const& particleComb,
                                                                int    const& ptoPL, 
                                                                int    const& nf,
                                                                double const& alphasAtQ);
}