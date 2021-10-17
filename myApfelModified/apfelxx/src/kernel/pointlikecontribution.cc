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
    const std::function<double(double const&)> PointlikeContribution(int const& particle,
                                                                     int const& PerturbativeOrderPointlike, 
                                                                     int const& nf)
    {
        switch (particle)
        {
        case EvolutionBasisQCD::Object::GLUON:
            return GluonPointlike(PerturbativeOrderPointlike, nf);
            break;
        case EvolutionBasisQCD::Object::SIGMA:
            return SigmaPointlike(PerturbativeOrderPointlike, nf);
            break
        default:
            return NonSingletPointlike(PerturbativeOrderPointlike, nf);
            break;
        }
    }

    const std::function<double(double const&)> SigmaPointlike(int const& PerturbativeOrderPointlike, 
                                                              int const& nf) // viellelicht einfach für alle gleich machen? also mit map. die vorfaktoren sind gleich
    {
        if (PerturbativeOrderPointlike == 0)
            return [&] (double const& x) -> double {return coeffQED * k0q(x, nf);};

        else if (PerturbativeOrderPointlike == 1)
            return [&] (double const& x) -> double {return coeffQED * ( k0q(x, nf) + coeffQCD * k1q(x, nf) );};
    }


    double k0q(double const& x, int const& nf)
    {
        return x * 3 * nf * eExp2 * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );
    }

    double k1q(double const& x, int const& nf)
    {
        return x * 3 * nf * eExp2 * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );
    }




}