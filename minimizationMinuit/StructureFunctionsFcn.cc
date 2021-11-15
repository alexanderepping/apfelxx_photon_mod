/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "apfel/apfelxx.h" 
#include "LHAPDF/LHAPDF.h"

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"

#include <string>
#include <map>
#include <vector>
#include <functional>
#include <math.h> // for gamma functions
#include <cassert> 



StructureFunctionsFcn::StructureFunctionsFcn(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                             std::string const& NameLHAPDFSet,
                                             double      const& ErrorDef): 
                _Energies(combineData(experimentalData, "Energies")),
                _xData(combineData(experimentalData, "xData")),
                _xError(combineData(experimentalData, "xError")),
                _F2Gamma(combineData(experimentalData, "F2Gamma")),
                _yError(combineData(experimentalData, "yError")),
                _NameLHAPDFSet(NameLHAPDFSet),
                _LHAPDFSet(LHAPDF::mkPDF(NameLHAPDFSet)),
                _ErrorDef(ErrorDef)
    {}

StructureFunctionsFcn::StructureFunctionsFcn(std::vector<double> const& Energies,
                                             std::vector<double> const& xData,
                                             std::vector<double> const& xError,
                                             std::vector<double> const& F2Gamma,
                                             std::vector<double> const& yError,
                                             std::string const& NameLHAPDFSet,
                                             double              const& ErrorDef): 
                _Energies(Energies),
                _xData(xData),
                _xError(xError),
                _F2Gamma(F2Gamma),
                _yError(yError),
                _NameLHAPDFSet(NameLHAPDFSet),
                _LHAPDFSet(LHAPDF::mkPDF(NameLHAPDFSet)),
                _ErrorDef(ErrorDef)
    {}


double StructureFunctionsFcn::operator()(std::vector<double> const& params) const 
{
    /////////////////////////////////////
    // Calculate Structure Function
    /////////////////////////////////////

    // Open LHAPDF set
    //LHAPDF::PDF* dist = LHAPDF::mkPDF("GRVCustomSetLO");  

    // Define Grid
    const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 6e-1, 3}, apfel::SubGrid{50, 8e-1, 5}}};

    // Define Thresholds and Masses
    const std::vector<double> Thresholds = {0, 0, 0, mc, mb};
    const std::vector<double> Masses = Thresholds;

    // Calculate AlphaS(Q)
    apfel::AlphaQCD a{asref, Qref, Masses, Thresholds, pto};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // Define Input PDFs
    auto InPDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { 
        switch (usedInitialPDFs)
        {
        case INITIALPDFS_9GDU_00:
            return InitialPDFs_9gdu_00(x, Q, params, _LHAPDFSet);
            break;

        case INITIALPDFS_9GDU_01:
            return InitialPDFs_9gdu_01(x, Q, params);
            break;

        case INITIALPDFS_3G:
            return InitialPDFs_3g(x, Q, params, _LHAPDFSet);
            break;

        case INITIALPDFS_2G:
            return InitialPDFs_2g(x, Q, params, _LHAPDFSet);
            break;
        
        default:
            break;
        }
         };

    // Initialise and Build DglapObjects
    auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs, Qin, pto, as);

    // Tabulate PDFs
    // const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3, as};

    // Effective EW charges for a space-like process (i.e. DIS)
    const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

    // Define PDFs functions 
    //const auto PDFs        = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); };
    const auto PDFsEvolved = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };

    // Initialise and Build Structure Functions
    //const auto F2        = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, pto, as, fBq);
    const auto F2Evolved = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolved, pto, as, fBq);

    // Tabulate Structure Functions
    //const apfel::TabulateObject<apfel::Distribution> F2total        {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2totalEvolved {[&] (double const& Q) -> apfel::Distribution{ return F2Evolved.at(0).Evaluate(Q);}, 50, 1, 1000, 3, Thresholds};

    const std::function<double (double const&, double const&)> F2EvolvedAtXQ = [&] (double const& x, double const& Q) -> double {return F2totalEvolved.EvaluatexQ(x, Q);};
    


    /////////////////////////////////////
    // Calulate chi2
    /////////////////////////////////////

    double chi2 = 0.;
    
    for (int i=0; i<Energies().size(); i++)
    { //unsure because of error and if i maybe need the matrix variant?!
        chi2 += ((F2EvolvedAtXQ(xData()[i], Energies()[i]) - F2Gamma()[i]) * (F2EvolvedAtXQ(xData()[i], Energies()[i]) - F2Gamma()[i]) / (yError()[i] * yError()[i]));
    }

    std::cout << "§ "; //debug
    for (int i=0; i<initialParams.at(usedInitialPDFs).size(); i++)
        std::cout << params[i] << ", ";
    std::cout << chi2 << std::endl;
    
    return chi2;
}



std::vector<double> StructureFunctionsFcn::combineData(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                                       std::string                                                       const& dataName) const
{
    std::vector<double> returnData{};

    for (std::string DataSet : IncludedExperimentalData)
        returnData.insert(std::end(returnData), std::begin(experimentalData.at(DataSet).at(dataName)), std::end(experimentalData.at(DataSet).at(dataName)));

    return returnData;
}




double StructureFunctionsFcn::betaFunction(double const & x, double const& y) const
{
    assert(x  >= 0.);
    assert(y  >= 0.);

    return (std::tgamma(x) * std::tgamma(y)) / std::tgamma(x + y);
}



std::map<int, double> StructureFunctionsFcn::InitialPDFs_9gdu_00(double                const& x,
                                                                 double                const& Q,
                                                                 std::vector<double>   const& params,
                                                                 LHAPDF::PDF*                 dist) const
    {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t

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

std::map<int, double> StructureFunctionsFcn::InitialPDFs_9gdu_01(double                const& x,
                                                                 double                const& Q,
                                                                 std::vector<double>   const& params) const
    {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t

        std::map<int, double> result;
        // top-quark
        result.insert(std::pair<int, double>( 6, 0.));
        result.insert(std::pair<int, double>(-6, result.at(6)));

        // bottom-quark
        result.insert(std::pair<int, double>( 5, 0.));
        result.insert(std::pair<int, double>(-5, result.at(5)));

        // charm-quark
        result.insert(std::pair<int, double>( 4, 0.));
        result.insert(std::pair<int, double>(-4, result.at(4)));

        // up-quark
        result.insert(std::pair<int, double>( 2, params[6] * pow(x, params[7]) * pow( (1.0 - x) , params[8])));
        result.insert(std::pair<int, double>(-2, result.at(2)));

        // down-quark
        result.insert(std::pair<int, double>( 1, params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5])));
        result.insert(std::pair<int, double>(-1, result.at(1)));

        // strange-quark
        result.insert(std::pair<int, double>( 3, params[0] / 2. * (result.at(2) + result.at(1))));
        result.insert(std::pair<int, double>(-3, result.at(3)));

        // gluon
        const double AN_g1 = (0.5 - (1 + params[0]/2.) * (  params[6] * betaFunction(params[7]+1, params[8]+1) 
                                                          + params[3] * betaFunction(params[4]+1, params[5]+1))) 
                             / betaFunction(params[1]+1, params[2]+1);
        result.insert(std::pair<int, double>( 0, AN_g1 * pow(x, params[1]) * pow( (1.0 - x) , params[2])));

        return apfel::PhysToQCDEv(result);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_3g(double                const& x,
                                                           double                const& Q,
                                                           std::vector<double>   const& params,
                                                           LHAPDF::PDF*                 dist) const
    {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t

        std::map<int, double> result;

        //result.insert(std::pair<int, double>(-6, dist->xfxQ(x, Q).at(-6)));
        result.insert(std::pair<int, double>(-5, dist->xfxQ(x, Q).at(-5)));
        result.insert(std::pair<int, double>(-4, dist->xfxQ(x, Q).at(-4)));
        result.insert(std::pair<int, double>(-3, dist->xfxQ(x, Q).at(-3)));
        result.insert(std::pair<int, double>(-2, dist->xfxQ(x, Q).at(-2)));
        result.insert(std::pair<int, double>(-1, dist->xfxQ(x, Q).at(-1)));
        result.insert(std::pair<int, double>( 0, 1./137. * params[0] * pow(x, params[1]) * pow( (1.0 - x) , params[2])));
        result.insert(std::pair<int, double>( 1, dist->xfxQ(x, Q).at(1)));
        result.insert(std::pair<int, double>( 2, dist->xfxQ(x, Q).at(2)));
        result.insert(std::pair<int, double>( 3, dist->xfxQ(x, Q).at(3)));
        result.insert(std::pair<int, double>( 4, dist->xfxQ(x, Q).at(4)));
        result.insert(std::pair<int, double>( 5, dist->xfxQ(x, Q).at(5)));
        //result.insert(std::pair<int, double>( 6, dist->xfxQ(x, Q).at(6)));

        return apfel::PhysToQCDEv(result);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_2g(double                const& x,
                                                           double                const& Q,
                                                           std::vector<double>   const& params,
                                                           LHAPDF::PDF*                 dist) const
    {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t

        std::map<int, double> result;

        //result.insert(std::pair<int, double>(-6, dist->xfxQ(x, Q).at(-6)));
        result.insert(std::pair<int, double>(-5, dist->xfxQ(x, Q).at(-5)));
        result.insert(std::pair<int, double>(-4, dist->xfxQ(x, Q).at(-4)));
        result.insert(std::pair<int, double>(-3, dist->xfxQ(x, Q).at(-3)));
        result.insert(std::pair<int, double>(-2, dist->xfxQ(x, Q).at(-2)));
        result.insert(std::pair<int, double>(-1, dist->xfxQ(x, Q).at(-1)));
        result.insert(std::pair<int, double>( 0, 1./137. * pow(x, params[0]) * pow( (1.0 - x) , params[1])));
        result.insert(std::pair<int, double>( 1, dist->xfxQ(x, Q).at(1)));
        result.insert(std::pair<int, double>( 2, dist->xfxQ(x, Q).at(2)));
        result.insert(std::pair<int, double>( 3, dist->xfxQ(x, Q).at(3)));
        result.insert(std::pair<int, double>( 4, dist->xfxQ(x, Q).at(4)));
        result.insert(std::pair<int, double>( 5, dist->xfxQ(x, Q).at(5)));
        //result.insert(std::pair<int, double>( 6, dist->xfxQ(x, Q).at(6)));

        return apfel::PhysToQCDEv(result);
    };