#pragma once

#include "apfel/apfelxx.h"
#include "LHAPDF/LHAPDF.h"

#include <map>
#include <vector>
#include <functional>
#include <cmath>


#ifdef vadimInitialPDFs
/**
 * @brief just some function to return the InitialPDFs from LHAPDF file. 
 * Uses 1 / alpha_em * an * x**a * (1-x)**b for gluon, up and down; the GRV PDFs for the rest.
 * 
 * @param x
 * @param Q
 * @param params: vector with 9 parameters: 
 *        AN_glu1 (0), A_glu1 (1), B_glu1 (2), AN_dbar1 (3), A_dbar1 (4), B_dbar1 (5), AN_ubar1 (6), A_ubar1 (7), B_ubar1 (8)
 */
const std::function<std::map<int, double>(double const&, double const&, std::vector<double> const&)> InitialPDFs = 
    [] (double const& x, double const& Q, std::vector<double> const& params) -> std::map<int, double>
    {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t

        LHAPDF::PDF* dist = LHAPDF::mkPDF("GRVCustomSetLO");  

        std::map<int, double> result;

        //result.insert(std::pair<int, double>(-6, dist->xfxQ(x, Q).at(-6)));
        result.insert(std::pair<int, double>(-5, dist->xfxQ(x, Q).at(-5)));
        result.insert(std::pair<int, double>(-4, dist->xfxQ(x, Q).at(-4)));
        result.insert(std::pair<int, double>(-3, dist->xfxQ(x, Q).at(-3)));
        result.insert(std::pair<int, double>(-2, 1./137. * params[6] * pow(x, params[7]) * pow( (1.0 - x) , params[8])));
        result.insert(std::pair<int, double>(-1, 1./137. * params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5])));
        result.insert(std::pair<int, double>( 0, 1./137. * params[0] * pow(x, params[1]) * pow( (1.0 - x) , params[2])));
        result.insert(std::pair<int, double>( 1, 1./137. * params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5])));
        result.insert(std::pair<int, double>( 2, 1./137. * params[6] * pow(x, params[7]) * pow( (1.0 - x) , params[8])));
        result.insert(std::pair<int, double>( 3, dist->xfxQ(x, Q).at(3)));
        result.insert(std::pair<int, double>( 4, dist->xfxQ(x, Q).at(4)));
        result.insert(std::pair<int, double>( 5, dist->xfxQ(x, Q).at(5)));
        //result.insert(std::pair<int, double>( 6, dist->xfxQ(x, Q).at(6)));

        return apfel::PhysToQCDEv(result);
    };
#endif //exampleInitialPDFs


#ifdef exampleInitialPDFs
/**
 * @brief just some function to return the InitialPDFs from LHAPDF file
 * 
 * @param x
 * @param Q
 * @param params: not used
 */
const std::function<std::map<int, double>(double const&, double const&, std::vector<double> const&)> InitialPDFs = 
    [] (double const& x, double const& Q, std::vector<double> const& params) -> std::map<int, double>
    {
        LHAPDF::PDF* dist = LHAPDF::mkPDF("GRVCustomSetLO");  
        return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); 
    };
#endif //exampleInitialPDFs