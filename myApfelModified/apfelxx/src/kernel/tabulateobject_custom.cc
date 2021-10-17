//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
// edited by: Alexander Epping: a_eppi01@uni-muenster.de
//

#include "apfel/tabulateobject.h"
#include "apfel/tabulateobject_custom.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"
#include "apfel/doubleobject.h"
#include "apfel/set.h"
#include "apfel/timer.h"

#include<algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  template<class T>
  // in Evolution.cc: TabulateObject<apfel::Set<apfel::Distribution>>
  // would mean, that Object is MatchedEvolution<Set<Distribution>>
  TabulateObject<T>::TabulateObject(MatchedEvolution<T>                       & Object,
                                    int                                  const& nQ,
                                    double                               const& QMin,
                                    double                               const& QMax,
                                    int                                  const& InterDegree,
                                    std::function<double(double const&)> const& Alphas,
                                    double                               const& Lambda):
    QGrid<T>
  {
    nQ, QMin, QMax, InterDegree, Object.GetThresholds(), Lambda
  }
  {
    report("Tabulating object... ");
    Timer t;

    // Save initial conditions.
    const int    nsteps = Object.GetNumberOfSteps();
    const T      ObjRef = Object.GetObjectRef();
    const double MuRef  = Object.GetMuRef();

    // Set number of steps of the RK algorith to 1.
    Object.SetNumberOfSteps(1);

    // Find the point on the QGrid right below MuRef.
    const int tQ = std::lower_bound(this->_Qg.begin(), this->_Qg.end(), MuRef) - this->_Qg.begin() - 1;

    // Loop on "_Qg" below "MuRef".
    for (int iQ = tQ; iQ >= 0; iQ--)
      {

        //std::cout << "\nCalling Evaluate with iQ="<<std::to_string(iQ);//debug
        const T o = Object.Evaluate(this->_Qg[iQ], Alphas); // result will be T = Set<Distribution>
        this->_GridValues.push_back(o);
        Object.SetObjectRef(o);
        Object.SetMuRef(this->_Qg[iQ]);
      }

    // Reverse order of the elements.
    std::reverse(this->_GridValues.begin(), this->_GridValues.end());

    // Loop on "_Qg" above "MuRef".
    Object.SetObjectRef(ObjRef);
    Object.SetMuRef(MuRef);
    for (int iQ = tQ + 1; iQ < (int) this->_Qg.size(); iQ++)
      {
        //std::cout << "\nCalling Evaluate with iQ="<<std::to_string(iQ);//debug
        const T o = Object.Evaluate(this->_Qg[iQ], Alphas);
        this->_GridValues.push_back(o);
        Object.SetObjectRef(o);
        Object.SetMuRef(this->_Qg[iQ]);
      }

    // Reset initial conditions.
    Object.SetNumberOfSteps(nsteps);
    Object.SetObjectRef(ObjRef);
    Object.SetMuRef(MuRef);

    t.stop();
  }
}
