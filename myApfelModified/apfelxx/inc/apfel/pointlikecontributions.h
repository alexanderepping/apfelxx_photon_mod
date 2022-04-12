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
// fixed alphas and coefficients for the whole contribution
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
     * @brief general coefficient 
     * @param x: x-value
     */
    double coeffGeneral(double const& x);


    /**
     * @brief QCD coefficient for the pointlike contribution:
     * alphaQCD(Q) / (2*Pi). 
     * @param alphasAtQ: value for the alphaQCD at specific Q
     */
    double coeffQCD(double const& alphasAtQ);




// some general definitions
    /**
     * @brief The map enumerators, describing the QCD evolution basis
     */
    enum ParticleCombinationEnumerator: int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The map enumerators, describing the different kinds of pointlike contributions. 
     * Used in
     */
    enum ContributionsEnumerator: int {GLUONS, SINGLET, NONSINGLET2, NONSINGLET3, NONSINGLET4, NONSINGLET5, NONSINGLET6};


    /**
     * @brief Quark charges squared.
     * {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {6, "t"}
     */
    const std::map<int, double> quarkCharges2 = {{1, 1./9.}, {2, 4./9.}, 
                                                 {3, 1./9.}, {4, 4./9.}, 
                                                 {5, 1./9.}, {6, 4./9.}};


    /**
     * @brief Map to represent which photon-parton splitting functions contribute
     * with what factor to which particle combination (in QCD evolution basis).
     * The order of the photon-parton splitting functions is given in ContributionsEnumerator.
     */
    const std::map<int, std::vector<int>> contributingKs = {{GLUON,   { 1,  0,  0,  0,  0,  0,  0}},
                                                            {SIGMA,   { 0,  1,  0,  0,  0,  0,  0}},
                                                            {VALENCE, { 0,  0,  0,  0,  0,  0,  0}},
                                                            {T3,      { 0,  0,  1,  0,  0,  0,  0}},
                                                            {V3,      { 0,  0,  0,  0,  0,  0,  0}},
                                                            {T8,      { 0,  0, -3, -2,  0,  0,  0}},
                                                            {V8,      { 0,  0,  0,  0,  0,  0,  0}},
                                                            {T15,     { 0,  0,  0, -2, -3,  0,  0}},
                                                            {V15,     { 0,  0,  0,  0,  0,  0,  0}},
                                                            {T24,     { 0,  0,  0,  0, -5, -4,  0}},
                                                            {V24,     { 0,  0,  0,  0,  0,  0,  0}},
                                                            {T35,     { 0,  0,  0,  0,  0, -4, -5}},
                                                            {V35,     { 0,  0,  0,  0,  0,  0,  0}}};




// Expectation Values of e^2 and e^4
    /**
     * @brief calculation of <e^2>.
     * @param nf: number of active flavours
     */
    double ExpE2(int const& nf);


    /**
     * @brief calculation of <e^4>.
     * @param nf: number of active flavours
     */
    double ExpE4(int const& nf);

    


// coefficients for the photon-parton splitting functions (k-functions)
    /**
     * @brief Coefficient used in the singlet and gluon k-functions. 
     * coeffSinglet = N_C * nf * <e^2>
     * @param nf: number of active flavours
     */
    double coeffSinglet(int const& nf);


    /**
     * @brief Coefficient used in the non singlet photon-parton splitting function (k-function). 
     * coeffNonSinglet = N_C * nf * (<e^4> - <e^2>^2)
     * @param nf: number of active flavours
     */
    double coeffNonSinglet(int const& nf);


    /**
     * @brief Coefficient used for calculating the T-functions in the QCD evolution basis.
     * It is multiplied with the non singlet photon-parton splitting functions.
     * coeffNonSingletContribution = 1 / (e_q^2 - <e^2>)
     * @param nf: number of active flavours
     */
    double coeffNonSingletContribution(int const& nf);




// functions used in the calculation of the pointlike contributions
    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     *      - note, that this one will be multiplied with coeffSinglet / coeffNonSinglet later to receive the functions from the paper
     * @brief k0q-function, used in the calculation of the Singlet and NonSinglet photon-parton splitting function at 0th perturbative order
     * @param x: x-value
     */
    double k0q(double const& x);


    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01) (function below (3.6), named k(x))
     * @brief k1q-function, used in the calculation of the Singlet and NonSinglet photon-parton splitting functions at 1st perturbative order
     * @param x: x-value
     */
    double k1q(double const& x);


    /**
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01)
     *      - it was wrong in Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01)
     *      - note, that this one will be multiplied with coeffSinglet later to receive the function from the paper
     * @brief k1g-function, used in the calculation of the gluon photon-parton splitting function at 1st perturbative order
     * @param x: x-value
     */
    double k1g(double const& x);


    /**
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. A11
     * @brief convolution of P0qq and Bgamma, 
     * used to calculate the additional DISgamma parts to the Singlet and NonSinglet photon-parton splitting functions at 1st perturbative order
     * @param x: x-value
     */
    double ConvolutionP0qqBgamma(double const& x);


    /**
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. A12
     * @brief convolution of P0gq and Bgamma, 
     * used to calculate the additional DISgamma parts to the Gluon photon-parton splitting function at 1st perturbative order
     * @param x: x-value
     */
    double ConvolutionP0gqBgamma(double const& x);




// photon-parton splitting functions for gluon
    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     * @brief photon-parton splitting function for the gluon at 0th perturbative order
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k0Gluon(double const& x, int const& nf);


    /**
     * - taken from Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. (2.7)
     * @brief photon-parton splitting function for the gluon at 1st perturbative order in the MSbar scheme
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k1GluonMS(double const& x, int const& nf);


    /**
     * - calculated with Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. 3.1
     * @brief photon-parton splitting function for the gluon at 1st perturbative order in the DISgamma scheme
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k1GluonDIS(double const& x, int const& nf);




// photon-parton splitting functions for quark singlet
    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     * @brief photon-parton splitting function for the singlet at 0th perturbative order
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k0Singlet(double const& x, int const& nf);


    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     * @brief photon-parton splitting function for the singlet at 1st perturbative order in the MSbar scheme
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k1SingletMS(double const& x, int const& nf);


    /**
     * - calculated with Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. 3.1
     * @brief photon-parton splitting function for the singlet at 1st perturbative order in the DISgamma scheme
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k1SingletDIS(double const& x, int const& nf);




// photon-parton splitting functions for quark non-singlet
    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     * @brief photon-parton splitting function for the non singlet at 0th perturbative order
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k0NonSinglet(double const& x, int const& nf);


    /**
     * - taken from Glück & Reya - Physical Review D, Volume 28, Number 11 (1983.12.01), eq. (3.6)
     * @brief photon-parton splitting function for the non singlet at 1st perturbative order in the MSbar scheme
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k1NonSingletMS(double const& x, int const& nf);


    /**
     * - calculated with Glück, Reya & Vogt - Physical Review D, Volume 45, Number 11 (1992.06.01), eq. 3.1
     * @brief photon-parton splitting function for the non singlet at 1st perturbative order in the DISgamma scheme
     * @param x: x-value
     * @param nf: number of active flavours
     */
    double k1NonSingletDIS(double const& x, int const& nf);

    


// calculation of the pointlike contribution
    /**
     * @brief Maps perturbative order and contribution to the correct photon-parton splitting function
     * @param contribution: which photon-parton splitting function should be returned
     * @param nf: number of active flavours
     * @param x: x-value
     */
    std::map<int, double> PointlikeMap(int const& contribution, int const& nf, double const& x);


    /**
     * @brief Pointlike contribution for specific particleCombination at
     * specific perturbative order.
     * @param particleComb: particleCombination for which pointlike contribution is to be calculated. particleCombination meaning in the QCD evolution basis
     * @param nf: number of active flavours
     * @param pto: Perturbative order 
     * @param alphasAtQ: the function returning the strong coupling
     */
    std::function<double(double const&)> PointlikeContribution (int    const& particleComb,
                                                                int    const& nf,
                                                                int    const& pto, 
                                                                double const& alphasAtQ);
}