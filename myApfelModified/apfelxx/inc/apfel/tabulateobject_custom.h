//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/qgrid.h"
#include "apfel/matchedevolution.h"

#include <map>

namespace apfel
{
  /**
   * @brief The template TabulateObject class is a derived of the
   * QGrid class that tabulates on object of type T (it can be a
   * double, a Distribution, an Operator, Set<Distribution>, a
   * Set<Operator>) over a grid in Q, taking into account the possible
   * presence of thresholds, and provides the method to evaluate the
   * tabulated object at any generic value of Q.
   */
  template<class T>
  class TabulateObject: public QGrid<T>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The TabulateObject constructor for an "evolving"
     * object (MatchedEvolution).
     * @param Object: the MatchedEvolution type object to be tabulated in Q
     * @param nQ: the number of on nodes of the grid in Q
     * @param QMin: the lower bound of the grid in Q
     * @param QMax: the upper bound of the grid in Q
     * @param InterDegree: the interpolation degree on the grid in Q
     * @param Alphas: the function returning the strong coupling
     * @param Lambda: the value of the parameter in the function ln(ln(Q<SUP>2</SUP>/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.25)
     */
    TabulateObject(MatchedEvolution<T>                       & Object,
                   int                                  const& nQ,
                   double                               const& QMin,
                   double                               const& QMax,
                   int                                  const& InterDegree,
                   std::function<double(double const&)> const& Alphas,
                   double                               const& Lambda = 0.25);
}
