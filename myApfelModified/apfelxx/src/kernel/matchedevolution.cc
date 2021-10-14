//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchedevolution.h"
#include "apfel/constants.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/ode.h"
#include "apfel/doubleobject.h"
#include "apfel/evolutionbasisqcd.h"

#include <cmath>
#include <algorithm>



namespace apfel 
  //addition, everything in this namespace
  // This class is defined because using T would result in errors. 
  // The reason for that is that T could also be e.g. double which doesn't have the .GetObjects() function.
  // The EvolveObject function adds the pointlike term to the rhs of the dglap evolution equqation. 
{
  //_________________________________________________________________________
  template<>
  MatchedEvolution<Set<Distribution>>::MatchedEvolution(Set<Distribution>   const& ObjRef,// is of type T = Set<Distribution>
                                                        double              const& MuRef,
                                                        std::vector<double> const& Thresholds,
                                                        int                 const& nsteps):
    _ObjRef(ObjRef),// is of type T = Set<Distribution>
    _MuRef(MuRef),
    _Thresholds(Thresholds),
    _nsteps(nsteps)
  {
    // Compute squared reference scale.
    _MuRef2 = pow(MuRef,2);

    // Compute log of the squared final scale.
    _LogMuRef2 = log(_MuRef2);

    // Compute squared thresholds.
    for (auto const& th : Thresholds)
      {
        const double th2 = pow(th, 2);
        _Thresholds2.push_back(th2);
        _LogThresholds2.push_back(( th2 > 0 ? log(th2) : -100));
      }

    // Sort the quark thresholds and logs.
    if (_Thresholds2.size() > 1)
      std::sort(_Thresholds2.begin(), _Thresholds2.end());
  }

  //_________________________________________________________________________________
  template<>
  Set<Distribution> MatchedEvolution<Set<Distribution>>::EvolveObject(int const& nf, double const& t0, double const& t1, Set<Distribution> const& Obj0) const
  // Obj0 is of type T = Set<Distribution>, therefore the return is also of T = Set<Distribution>
  // t0 and t1 are the initial and final scale (ln(Q^2)) of the evolution
  {
    // Return immediately "Obj0" if "t0" and "t1" are equal.
    if (t0 == t1)
      return Obj0;

    //std::cout << "\nCalled Evolve Object";//debug


    const std::function<Set<Distribution>(double const&, Set<Distribution> const&)> f = [&] (double const& t, Set<Distribution> const& Obj) -> Set<Distribution>
      {
        // rright hand side of homogeneous DGLAP evolution equation
        Set<Distribution> rhsGeneral = Derivative(nf, t, Obj);
        
        // add inhomogeneous part to the homogeneous 
        std::map<int, Distribution> tempMap;

        // loop through map of rhsGeneral and only add the inhomogeneous part (Distributions) at some keys
        for (int it=0; it!=static_cast<int>(rhsGeneral.GetObjects().size()); ++it)
          {
            Grid const& rhsGrid = rhsGeneral.GetObjects().at(it).GetGrid();

            switch (it)
            {
            case EvolutionBasisQCD::Object::GLUON:
              {
                tempMap.insert({it, rhsGeneral.GetObjects().at(it)});
                break;
              }
            case EvolutionBasisQCD::Object::SIGMA:
              {
                Distribution rhsParicular(rhsGrid, [&] (double const& x) -> double  { const double eExp2    = 10. / 9. / 4. ;
                                                                                      //const double alphaQED = 1./137.;
                                                                                      const double k0q      = x * 3 * nf * eExp2 * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );
                                                                                      //return alphaQED / (2 * M_PI) * k0q; });
                                                                                      //return alphaQED*alphaQED / (2 * M_PI) * k0q; });
                                                                                      return 1 / (2 * M_PI) * k0q; });

                // insert general + particular solution into tempMap
                tempMap.insert({it, rhsGeneral.GetObjects().at(it) + rhsParicular});
                break;
              }
            default:
              {
                Distribution rhsParicular(rhsGrid, [&] (double const& x) -> double {  const double eExp2    = 10. / 9. / 4. ;
                                                                                      const double eExp4    = 34. / 81. / 4. ; 
                                                                                      //const double alphaQED = 1./137.;
                                                                                      const double k0ns     = x * 3 * nf * ( eExp4 - eExp2 * eExp2 ) * 2 * ( x * x + ( 1 - x ) * ( 1 - x ) );
                                                                                      //return alphaQED / (2 * M_PI) * k0ns; });
                                                                                      //return alphaQED*alphaQED / (2 * M_PI) * k0ns; });
                                                                                      return 1 / (2 * M_PI) * k0ns; });

                // insert general + particular solution into tempMap
                tempMap.insert({it, rhsGeneral.GetObjects().at(it) + rhsParicular});
                break;
              }
            }
          }
                    
        rhsGeneral.SetObjects(tempMap);

        return rhsGeneral;
      };

    const auto dObj = rk4<Set<Distribution>>(f);
    //const auto dObj = rk4Photon<T>([&] (double const& t, T const& Obj) -> T{ return Derivative(nf, t, Obj); }); //debug


    // T is obviously the same as U. The type of the input object is T = Set<Distribution>,
    // meaning the function Derivative also returns a Set<Distribution>.

    // Use "_nsteps" steps for the evolution.
    double t        = t0;
    Set<Distribution>      Obj      = Obj0;
    const double dt = ( t1 - t0 ) / _nsteps;
    for (int k = 0; k < _nsteps; k++)
      {
        Obj += dObj(t, Obj, dt);
        t   += dt;
      }
    return Obj;
  }

  //_________________________________________________________________________
  template<>
  Set<Distribution> MatchedEvolution<Set<Distribution>>::Evaluate(double const& mu) const // result will be T = Set<Distribution>
  {
    //std::cout << "\nCalled Evaluate";//debug
    const double mu2  = pow(mu,2);
    const double lmu2 = log(mu2);

    // Find initial and final number of flavours.
    const int nfi = NF(_MuRef2, _Thresholds2);
    const int nff = NF(mu2, _Thresholds2);

    // Don't do the matching is initial and final number of flavours
    // are equal.
    if (nfi == nff)
      return EvolveObject(nfi, _LogMuRef2, lmu2, _ObjRef);

    // Direction of the evolution
    const bool sgn = std::signbit(nfi - nff);

    // Create a vector of objects containing the object right above
    // each threshold to make sure that every time a threshold is
    // crossed a new object with a different convolution map is
    // created (effective only when a "Set" object is evolved).
    Set<Distribution>      vobj = _ObjRef;
    double ti   = _LogMuRef2;
    for (int inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
        // Final scale
        const double tf = _LogThresholds2[(sgn ? inf : inf - 1)];

        // Do the matching
        vobj = MatchObject(sgn, inf, EvolveObject(inf, ti, tf, vobj));

        // Update initial scale and displace it by "eps8" to make sure
        // to be above (below) the threshold
        ti = tf * ( 1 + (sgn ? 1 : -1) * eps8 );
      }
    return EvolveObject(nff, ti, lmu2, vobj);
  }
}

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  MatchedEvolution<T>::MatchedEvolution(T                   const& ObjRef,// is of type T = Set<Distribution>
                                        double              const& MuRef,
                                        std::vector<double> const& Thresholds,
                                        int                 const& nsteps):
    _ObjRef(ObjRef),// is of type T = Set<Distribution>
    _MuRef(MuRef),
    _Thresholds(Thresholds),
    _nsteps(nsteps)
  {
    // Compute squared reference scale.
    _MuRef2 = pow(MuRef,2);

    // Compute log of the squared final scale.
    _LogMuRef2 = log(_MuRef2);

    // Compute squared thresholds.
    for (auto const& th : Thresholds)
      {
        const double th2 = pow(th, 2);
        _Thresholds2.push_back(th2);
        _LogThresholds2.push_back(( th2 > 0 ? log(th2) : -100));
      }

    // Sort the quark thresholds and logs.
    if (_Thresholds2.size() > 1)
      std::sort(_Thresholds2.begin(), _Thresholds2.end());
  }

  //_________________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::EvolveObject(int const& nf, double const& t0, double const& t1, T const& Obj0) const
  // Obj0 is of type T = Set<Distribution>, therefore the return is also of T = Set<Distribution>
  // t0 and t1 are the initial and final scale (ln(Q^2)) of the evolution
  {
    // Return immediately "Obj0" if "t0" and "t1" are equal.
    if (t0 == t1)
      return Obj0;

    // Numerical solution of the evolution equation with fourth-order
    // Runge-Kutta.
    // The output of rk4  is std::function<U(double const&, U const&, double const&)>,
    // meaning a function of double const& t, T const& Obj and double const& dt.
    // It takes std::function<U(double const& t, U const& Obj)> const& f as an input,
    // meaning it requires an const& f of the type function<U(double const& t, U const& Obj)>,
    // where the function is the lambda function, the U in front of the brackets means 
    // the output type of the function and t and Obj are the inputs.
    const auto dObj = rk4<T>([&] (double const& t, T const& Obj) -> T{ return Derivative(nf, t, Obj); });
    // T is obviously the same as U. The type of the input object is T = Set<Distribution>,
    // meaning the function Derivative also returns a Set<Distribution>.

    // Use "_nsteps" steps for the evolution.
    double t        = t0;
    T      Obj      = Obj0;
    const double dt = ( t1 - t0 ) / _nsteps;
    for (int k = 0; k < _nsteps; k++)
      {
        Obj += dObj(t, Obj, dt);
        t   += dt;
      }
    return Obj;
  }

  //_________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::Evaluate(double const& mu) const // result will be T = Set<Distribution>
  {
    //std::cout << "\nCalled Evaluate";//debug
    const double mu2  = pow(mu,2);
    const double lmu2 = log(mu2);

    // Find initial and final number of flavours.
    const int nfi = NF(_MuRef2, _Thresholds2);
    const int nff = NF(mu2, _Thresholds2);

    // Don't do the matching is initial and final number of flavours
    // are equal.
    if (nfi == nff)
      return EvolveObject(nfi, _LogMuRef2, lmu2, _ObjRef);

    // Direction of the evolution
    const bool sgn = std::signbit(nfi - nff);

    // Create a vector of objects containing the object right above
    // each threshold to make sure that every time a threshold is
    // crossed a new object with a different convolution map is
    // created (effective only when a "Set" object is evolved).
    T      vobj = _ObjRef;
    double ti   = _LogMuRef2;
    for (int inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
        // Final scale
        const double tf = _LogThresholds2[(sgn ? inf : inf - 1)];

        // Do the matching
        vobj = MatchObject(sgn, inf, EvolveObject(inf, ti, tf, vobj));

        // Update initial scale and displace it by "eps8" to make sure
        // to be above (below) the threshold
        ti = tf * ( 1 + (sgn ? 1 : -1) * eps8 );
      }
    return EvolveObject(nff, ti, lmu2, vobj);
  }

  // template fixed types
  template class MatchedEvolution<double>;                                    //<! Single coupling
  template class MatchedEvolution<Distribution>;                              //<! Single distribution
  //template class MatchedEvolution<Set<Distribution>>;                         //<! Set of distributions //addition
  template class MatchedEvolution<DoubleObject<Distribution>>;                //<! Double object of distributions
  template class MatchedEvolution<Operator>;                                  //<! Single Operator
  template class MatchedEvolution<Set<Operator>>;                             //<! Set of Operators
  template class MatchedEvolution<DoubleObject<Operator>>;                    //<! Double object of operators
  template class MatchedEvolution<DoubleObject<Distribution, Operator>>;      //<! Double object of distributions and operators
  template class MatchedEvolution<Set<DoubleObject<Distribution, Operator>>>; //<! Set of double object of distributions and operators
}
