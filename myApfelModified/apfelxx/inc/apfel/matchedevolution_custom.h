//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <cmath>
#include <vector>

namespace apfel
{
  /**
   * @brief The MatchedEvolution class is a template mother class for
   * the computation of the running of a generic quantity in a
   * VFNS. It provides the basic ingredients for the computation and
   * the heavy-quark threshold matching of the running of a given
   * object.
   */
  template<>
  class MatchedEvolution
  {
  public:
    /**
     * @brief Virtual function for the computation of the evolution.
     * @param nf: the number of active flavours
     * @param Obj0: the starting object
     * @param mu02: the squared starting scale
     * @param mu2: the squared final scale
     * @param Alphas: the function returning the strong coupling
     * @return the object evolved at the scale mu2
     */
    virtual Set<Distribution> EvolveObject(int const& nf, double const& mu02, double const& mu2, T const& Obj0, std::function<double(double const&)> const& Alphas) const;

    /**
     * @brief Function that returns the evolved object.
     * @param mu: the final scale
     * @param Alphas: the function returning the strong coupling
     * @return the evolved object.
     */
    Set<Distribution> Evaluate(double const& mu, std::function<double(double const&)> const& Alphas) const;
}
