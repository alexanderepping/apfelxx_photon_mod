// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_GaussRandomGen_H_
#define MN_GaussRandomGen_H_

#include <cmath>
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

namespace MinuitCpp {

class GaussRandomGen {

public:

  GaussRandomGen() : fMean(0.), fSigma(1.) {}

  GaussRandomGen(double mean, double sigma) : fMean(mean), fSigma(sigma) {}

  ~GaussRandomGen() {}

  double Mean() const {return fMean;}

  double Sigma() const {return fSigma;}

  double operator()() const {
    //need to random variables flat in [0,1)
    double r1 = std::rand()/double(RAND_MAX);
    double r2 = std::rand()/double(RAND_MAX);

    //two possibilities to generate a random gauss variable (m=0,s=1)
    double s = sqrt(-2.*log(1.-r1))*cos(2.*M_PI*r2);
//     double s = sqrt(-2.*log(1.-r1))*sin(2.*M_PI*r2);

    //scale to desired gauss
    return Sigma()*s + Mean();
  }

private:

  double fMean;
  double fSigma;
  
};

}  // namespace Minuit2

#endif //MN_GaussRandomGen_H_