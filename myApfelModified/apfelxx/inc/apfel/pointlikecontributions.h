//
// APFEL++ 2017
//
// Author: Alexander Epping: a_eppi01@uni-muenster.de
//

namespace apfel
{
    
    const int PerturbativeOrderPointlike = 0;
    const double eExp2    = 10. / 9. / 4. ;
    const double eExp4    = 34. / 81. / 4. ; 
    const double alphaQED = 1./137.;
    const double coeffQED = 1. / (2. * M_PI) ;
    const double coeffQCD = 0;


    const std::function<double(double const&)> PointlikeContribution(int    const& particle,
                                                                     int    const& PerturbativeOrderPointlike, 
                                                                     int    const& nf);

    
}