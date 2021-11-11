#include "apfel/apfelxx.h" 
#include "LHAPDF/LHAPDF.h"

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"
#include "InitialPDFs.h"

#include <string>
#include <map>
#include <vector>
#include <functional>


StructureFunctionsFcn::StructureFunctionsFcn(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                             double const& ErrorDef): 
                _Energies(combineData(experimentalData, "Energies")),
                _xData(combineData(experimentalData, "xData")),
                _xError(combineData(experimentalData, "xError")),
                _F2Gamma(combineData(experimentalData, "F2Gamma")),
                _yError(combineData(experimentalData, "yError")),
                _ErrorDef(ErrorDef)
    {}

StructureFunctionsFcn::StructureFunctionsFcn(std::vector<double> const& Energies,
                                             std::vector<double> const& xData,
                                             std::vector<double> const& xError,
                                             std::vector<double> const& F2Gamma,
                                             std::vector<double> const& yError,
                                             double              const& ErrorDef): 
                _Energies(Energies),
                _xData(xData),
                _xError(xError),
                _F2Gamma(F2Gamma),
                _yError(yError),
                _ErrorDef(ErrorDef)
    {}


double StructureFunctionsFcn::operator()(std::vector<double> const& params) const 
{
    /////////////////////////////////////
    // Calculate Structure Function
    /////////////////////////////////////

    // Open LHAPDF set
    LHAPDF::PDF* dist = LHAPDF::mkPDF(NameLHAPDFSet);   

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
        /* //debug
        if (Q != 1.29)
        {
            std::cout << "§ Q=" << Q << std::endl;
        }*/
        
        return InitialPDFs(x, Q, params);};

    // Initialise and Build DglapObjects
    auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs, Qin, pto, as);

    // Tabulate PDFs
    // const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3, as};

    // Effective EW charges for a space-like process (i.e. DIS)
    const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

    // Define PDFs functions 
    // Cutoff not used!!
    const auto PDFs        = [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(dist->xfxQ(x, Q)); };
    const auto PDFsEvolved = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };

    // Initialise and Build Structure Functions
    const auto F2        = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, pto, as, fBq);
    const auto F2Evolved = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolved, pto, as, fBq);

    // Tabulate Structure Functions
    const apfel::TabulateObject<apfel::Distribution> F2total        {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
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