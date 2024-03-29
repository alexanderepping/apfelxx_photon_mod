/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "apfel/apfelxx.h"  

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"
#include "HelperFunctions.h"

#include <string>
#include <map>
#include <vector>
#include <functional>
#include <math.h> // for gamma functions
#include <cmath> // for the ExpIntegralEi/expint Well, I started a new project and I use <string> again. Now it works. 
#include <cassert> 



StructureFunctionsFcn::StructureFunctionsFcn(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                             std::vector<std::string>                                          const& IncludedExpData,
                                             double                                                            const& ErrorDef): 
                _ExperimentalData(experimentalData),
                _IncludedExpData(IncludedExpData),
                _ErrorDef(ErrorDef)
    {}



double StructureFunctionsFcn::operator()(std::vector<double> const& params) const 
{
    return Chi2PerExperiment(params).at("All");
}



std::map<std::string, double> StructureFunctionsFcn::Chi2PerExperiment(std::vector<double> const& params) const 
{
    /////////////////////////////////////
    // Calculate Structure Function
    /////////////////////////////////////

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
    auto InPDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return InitialPDFs(x, Q, params); };

    // Initialise and Build DglapObjects
    auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Masses, Thresholds), InPDFs, Qin, pto, as);

    // Tabulate PDFs
    // const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3, as};
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3, pto, as}; //notTestedYet

    // Effective EW charges for a space-like process (i.e. DIS)
    const std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

    // Define PDFs functions 
    const auto PDFsEvolved = [&] (double const& x, double const& Q) -> std::map<int, double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };

    // Initialise and Build Structure Functions
    //const auto F2        = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFs, pto, as, fBq);
    const auto F2Evolved = BuildStructureFunctions(InitializeF2NCObjectsZM(g, Thresholds), PDFsEvolved, pto, as, fBq);

    // Tabulate Structure Functions
    //const apfel::TabulateObject<apfel::Distribution> F2total        {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2totalEvolved {[&] (double const& Q) -> apfel::Distribution{ return F2Evolved.at(0).Evaluate(Q);}, 50, 1, 1000, 3, Thresholds};

    const std::function<double (double const&, double const&)> F2EvolvedAtXQ = [&] (double const& x, double const& Q) -> double {return F2totalEvolved.EvaluatexQ(x, Q);};
    


    /////////////////////////////////////
    // Calulate chi2PerExperiment
    /////////////////////////////////////

    std::map<std::string, double> chi2PerExperiment;
    chi2PerExperiment["All"] = 0.;

    for (std::string DataSet : IncludedExpData())
    {
        chi2PerExperiment[DataSet] = 0.;

        const std::vector<double> Q2Data     = ExperimentalData().at(DataSet).at("Q2Data");
        const std::vector<double> xData      = ExperimentalData().at(DataSet).at("xData");
        const std::vector<double> F2Gamma    = ExperimentalData().at(DataSet).at("F2Gamma");
        const std::vector<double> F2GammaErr = ExperimentalData().at(DataSet).at("F2GammaErr");


        for (int i=0; i<Q2Data.size(); i++)
        {
            const double QData = pow(Q2Data[i], 0.5);
            chi2PerExperiment[DataSet] += ((F2EvolvedAtXQ(xData[i], QData) - F2Gamma[i]) * (F2EvolvedAtXQ(xData[i], QData) - F2Gamma[i]) / (F2GammaErr[i] * F2GammaErr[i]));
        }
        chi2PerExperiment["All"] += chi2PerExperiment[DataSet];
    }
    
    return chi2PerExperiment;
}



std::map<int, double> StructureFunctionsFcn::InitialPDFs(double              const& x,
                                                         double              const& Q,
                                                         std::vector<double> const& params,
                                                         bool                const& returnParameters) const
{   
    switch (usedInitialPDFs)
    {
    case INITIALPDFS_9GDUS:
        return InitialPDFs_9gdus(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_8GDU:
        return InitialPDFs_8gdu(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_6GQS:
        return InitialPDFs_6gqs(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_5GQ:
        return InitialPDFs_5gq(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_SAL8:
        return InitialPDFs_SAL8(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_SAL6:
        return InitialPDFs_SAL6(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_SAL5:
        return InitialPDFs_SAL5(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_SAL4:
        return InitialPDFs_SAL4(x, Q, params, returnParameters);
        break;
    case INITIALPDFS_SAL3:
        return InitialPDFs_SAL3(x, Q, params, returnParameters);
        break;
    default:
        abort();
    } 
}



std::vector<double> StructureFunctionsFcn::combineData(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                                       std::vector<std::string>                                          const& IncludedExpData,
                                                       std::string                                                       const& dataName) const
{
    std::vector<double> returnData{};

    for (std::string DataSet : IncludedExpData)
        returnData.insert(std::end(returnData), std::begin(experimentalData.at(DataSet).at(dataName)), std::end(experimentalData.at(DataSet).at(dataName)));

    return returnData;
}



double StructureFunctionsFcn::betaFunction(double const & x, double const& y) const
{
    //assert(x  >= 0. & y  >= 0.);
    //DebugString("betaFunction - x = "+std::to_string(x)+", y = "+std::to_string(y)); //debug
    if (x <= 0. || y <= 0.)
        DebugString("betaFunction - x = "+std::to_string(x)+", y = "+std::to_string(y)); //debug

    return (std::tgamma(x) * std::tgamma(y)) / std::tgamma(x + y);
}



double StructureFunctionsFcn::MomentumSumRule0(std::vector<double> const& params,
                                              double const& totMom) const
{
    return ( totMom - (2 + params[0]) * (  params[6] * betaFunction(params[7]+1, params[8]+1) 
                                                + params[3] * betaFunction(params[4]+1, params[5]+1))) 
            / betaFunction(params[1]+1, params[2]+1);
}



double StructureFunctionsFcn::MomentumSumRuleSAL(std::vector<double> const& params,
                                                 double const& Q) const 
{//K_S (0), B_G_HAD(1), C_G_HAD(2), A_Q_HAD(3), B_Q_HAD(4), C_Q_HAD(5), A_Q_PL(6), B_Q_PL(7)

    const std::function<double(int const&, double const&)> helperFunction = [&] (int const& n, double const& B_Q_PL) -> double {return std::exp(n / B_Q_PL) * std::expint(-n / B_Q_PL);};



    const double EQ2_u = 4./9.; // electric charge squared of the up quark
    const double EQ2_d = 1./9.; // electric charge squared of the down quark
    const double EQ2_s = 1./9.; // electric charge squared of the strange quark

    // right hand side / result of the Momentum Sum Rule
    const double rhs = 1. + 2. / (3. * M_PI) * log( Q * Q / 4. ); 

    // Integral over x x^B_G_HAD (1-x)^C_G_HAD
    const double lhsHadG = betaFunction(params[1]+1, params[2]+1); 

    // Integral over A_Q_HAD x x^B_Q_HAD (1-x)^C_Q_HAD
    const double lhsHadQ = params[3] * betaFunction(params[4]+1, params[5]+1); 

    // Integral over A_Q_PL x (x^2 + (1-x)^2) / (1 - B_Q_PL ln(1-x)a
    // depending on, if B_Q_PL is zero or not, differnt values for lhsPLQ
    double lhsPLQ; 
    if (params[7] != 0)
        lhsPLQ  = -1 * params[6] / params[7] * ( helperFunction(1,params[7]) - 3 * helperFunction(2,params[7]) + 4 * helperFunction(3,params[7]) - 2 * helperFunction(4,params[7]) );
    else
        lhsPLQ  = params[6] / 3.;

    const double A_G_Had = (rhs - 2. * ((1.+1.+params[0])*lhsHadQ + (EQ2_u+EQ2_d+EQ2_s)*lhsPLQ)) / lhsHadG;

    DebugString("Parameters:", true, true, 3); //debug
    DebugVectorDoubles(params, true, true, 3); //debug

    DebugString("lhsHadG, lhsHadQ, lhsPLQ, rhs: ", true, true, 3);
    DebugString(std::to_string(lhsHadG) + ", " + std::to_string(lhsHadQ) + ", " + std::to_string(lhsPLQ) + ", " + std::to_string(rhs), true, true, 3);

    if (params[7] != 0)
    {
        DebugString("helperFunction(1,params[7]), helperFunction(2,params[7]), helperFunction(3,params[7]), helperFunction(4,params[7]): ", true, true, 3);
        DebugString(std::to_string(helperFunction(1,params[7])) + ", " + std::to_string(helperFunction(2,params[7])) + ", " + std::to_string(helperFunction(3,params[7])) + ", " + std::to_string(helperFunction(4,params[7])), true, true, 3);
    }
    DebugString("", true, false, 3); //debug

    return A_G_Had;
}



std::map<int, double> StructureFunctionsFcn::InitialPDFsMain0(double                const& x,
                                                              double                const& Q,
                                                              std::vector<double>   const& params,
                                                              bool                  const& outputAN_g1) const
    {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t
        // parameters:  K_s1 (0), A_g1 (1), B_g1 (2), AN_d1 (3), A_d1 (4), B_d1 (5), AN_u1 (6), A_u1 (7), B_u1 (8)

        std::map<int, double> result;

        const double AN_g1 = MomentumSumRule0(params);

        result.insert(std::pair<int, double>( 6, 0.));                                                          // top-quark
        result.insert(std::pair<int, double>( 5, 0.));                                                          // bottom-quark
        result.insert(std::pair<int, double>( 4, 0.));                                                          // charm-quark
        result.insert(std::pair<int, double>( 2, params[6] * pow(x, params[7]) * pow( (1.0 - x) , params[8]))); // up-quark
        result.insert(std::pair<int, double>( 1, params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5]))); // down-quark
        result.insert(std::pair<int, double>( 3, params[0] / 2. * (result.at(2) + result.at(1))));              // strange-quark
        result.insert(std::pair<int, double>( 0, AN_g1 * pow(x, params[1]) * pow( (1.0 - x) , params[2])));     // gluon

        // anti-quarks
        for (int i=1; i<7; i++)
            result.insert(std::pair<int, double>(-1*i, result.at(i)));
        
        // output AN_g1
        if (outputAN_g1)
            std::cout << "§ AN_g1\t = "<< AN_g1 << std::endl; // debug

        return apfel::PhysToQCDEv(result);
    };



 std::map<int, double> StructureFunctionsFcn::InitialPDFsMainSAL(double                const& x,
                                                                 double                const& Q,
                                                                 std::vector<double>   const& params,
                                                                 bool                  const& outputA_G_HAD) const
     {   // particles:   0: gluon, 1: d, 2: u, 3: s, 4: c, 5: b, 6: t
        // parameters:  K_S (0), B_G_HAD(1), C_G_HAD(2), A_Q_HAD(3), B_Q_HAD(4), C_Q_HAD(5), A_Q_PL(6), B_Q_PL(7)

        std::map<int, double> result;

        const double EQ2_s = 1./9.; // electric charge squared of the strange quark
        const double EQ2_u = 4./9.; // electric charge squared of the up quark
        const double EQ2_d = 1./9.; // electric charge squared of the down quark

        const double A_G_HAD = MomentumSumRuleSAL(params, Q);

        result.insert(std::pair<int, double>( 6, 0.));                                                                                  // top-quark
        result.insert(std::pair<int, double>( 5, 0.));                                                                                  // bottom-quark
        result.insert(std::pair<int, double>( 4, 0.));                                                                                  // charm-quark                                                                
        result.insert(std::pair<int, double>( 3, params[0] * params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5])
                                                + EQ2_s * params[6] * ( x * x + (1 - x) * (1 - x) ) / ( 1 - params[7] * log(1 - x))));  // strange-quark
        result.insert(std::pair<int, double>( 2, params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5])
                                                + EQ2_u * params[6] * ( x * x + (1 - x) * (1 - x) ) / ( 1 - params[7] * log(1 - x))));  // up-quark
        result.insert(std::pair<int, double>( 1, params[3] * pow(x, params[4]) * pow( (1.0 - x) , params[5])
                                                + EQ2_d * params[6] * ( x * x + (1 - x) * (1 - x) ) / ( 1 - params[7] * log(1 - x))));  // down-quark
        result.insert(std::pair<int, double>( 0, A_G_HAD  * pow(x, params[1]) * pow( (1.0 - x) , params[2])));                          // gluon

        // anti-quarks
        for (int i=1; i<7; i++)
            result.insert(std::pair<int, double>(-1*i, result.at(i)));
        
        // output A_G_HAD
        if (outputA_G_HAD)
            std::cout << "§ A_G_HAD\t = "<< A_G_HAD << std::endl; // debug

        return apfel::PhysToQCDEv(result);
    };



std::map<int, double> StructureFunctionsFcn::InitialPDFs_9gdus(double              const& x,
                                                               double              const& Q,
                                                               std::vector<double> const& params,
                                                               bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(params[0]);            //parameters[0] = K_s1
        parameters.push_back(params[1]);            //parameters[1] = A_g1
        parameters.push_back(params[2]);            //parameters[2] = B_g1
        parameters.push_back(params[3]);            //parameters[3] = AN_d1
        parameters.push_back(params[4]);            //parameters[4] = A_d1
        parameters.push_back(params[5]);            //parameters[5] = B_d1
        parameters.push_back(params[6]);            //parameters[6] = AN_u1
        parameters.push_back(params[7]);            //parameters[7] = A_u1
        parameters.push_back(params[8]);            //parameters[8] = B_u1
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMain0(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_8gdu(double              const& x,
                                                              double              const& Q,
                                                              std::vector<double> const& params,
                                                              bool                const& returnParameters) const
    {   
        std::vector<double> parameters;
        parameters.push_back(0.5);                  //parameters[0] = K_s1 //0.5 because the strange quark should contribute less than 1/2(u+d)
        parameters.push_back(params[0]);            //parameters[1] = A_g1
        parameters.push_back(params[1]);            //parameters[2] = B_g1
        parameters.push_back(params[2]);            //parameters[3] = AN_d1
        parameters.push_back(params[3]);            //parameters[4] = A_d1
        parameters.push_back(params[4]);            //parameters[5] = B_d1
        parameters.push_back(params[5]);            //parameters[6] = AN_u1
        parameters.push_back(params[6]);            //parameters[7] = A_u1
        parameters.push_back(params[7]);            //parameters[8] = B_u1
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMain0(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_6gqs(double              const& x,
                                                              double              const& Q,
                                                              std::vector<double> const& params,
                                                              bool                const& returnParameters) const
    {   
        std::vector<double> parameters;
        parameters.push_back(params[0]);            //parameters[0] = K_s1
        parameters.push_back(params[1]);            //parameters[1] = A_g1
        parameters.push_back(params[2]);            //parameters[2] = B_g1
        parameters.push_back(params[3]);            //parameters[3] = AN_d1
        parameters.push_back(params[4]);            //parameters[4] = A_d1
        parameters.push_back(params[5]);            //parameters[5] = B_d1
        parameters.push_back(params[3]);            //parameters[6] = AN_u1
        parameters.push_back(params[4]);            //parameters[7] = A_u1
        parameters.push_back(params[5]);            //parameters[8] = B_u1
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMain0(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_5gq(double              const& x,
                                                             double              const& Q,
                                                             std::vector<double> const& params,
                                                             bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(0.5);                  //parameters[0] = K_s1 //0.5 because the strange quark should contribute less than 1/2(u+d)
        parameters.push_back(params[0]);            //parameters[1] = A_g1
        parameters.push_back(params[1]);            //parameters[2] = B_g1
        parameters.push_back(params[2]);            //parameters[3] = AN_d1
        parameters.push_back(params[3]);            //parameters[4] = A_d1
        parameters.push_back(params[4]);            //parameters[5] = B_d1
        parameters.push_back(params[2]);            //parameters[6] = AN_u1
        parameters.push_back(params[3]);            //parameters[7] = A_u1
        parameters.push_back(params[4]);            //parameters[8] = B_u1
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMain0(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_SAL8(double             const& x,
                                                             double              const& Q,
                                                             std::vector<double> const& params,
                                                             bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(params[0]);            //parameters[0] = K_S
        parameters.push_back(params[1]);            //parameters[1] = B_G_HAD
        parameters.push_back(params[2]);            //parameters[2] = C_G_HAD
        parameters.push_back(params[3]);            //parameters[3] = A_Q_HAD
        parameters.push_back(params[4]);            //parameters[4] = B_Q_HAD
        parameters.push_back(params[5]);            //parameters[5] = C_Q_HAD
        parameters.push_back(params[6]);            //parameters[6] = A_Q_PL
        parameters.push_back(params[7]);            //parameters[7] = B_Q_PL
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMainSAL(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_SAL6(double             const& x,
                                                             double              const& Q,
                                                             std::vector<double> const& params,
                                                             bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(params[0]);            //parameters[0] = K_S
        parameters.push_back(params[1]);            //parameters[1] = B_G_HAD
        parameters.push_back(3.);                   //parameters[2] = C_G_HAD
        parameters.push_back(params[2]);            //parameters[3] = A_Q_HAD
        parameters.push_back(params[3]);            //parameters[4] = B_Q_HAD
        parameters.push_back(1.);                   //parameters[5] = C_Q_HAD
        parameters.push_back(params[4]);            //parameters[6] = A_Q_PL
        parameters.push_back(params[5]);            //parameters[7] = B_Q_PL
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMainSAL(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_SAL5(double             const& x,
                                                             double              const& Q,
                                                             std::vector<double> const& params,
                                                             bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(0.3);                  //parameters[0] = K_S
        parameters.push_back(params[0]);            //parameters[1] = B_G_HAD
        parameters.push_back(3.);                   //parameters[2] = C_G_HAD
        parameters.push_back(params[1]);            //parameters[3] = A_Q_HAD
        parameters.push_back(params[2]);            //parameters[4] = B_Q_HAD
        parameters.push_back(1.);                   //parameters[5] = C_Q_HAD
        parameters.push_back(params[3]);            //parameters[6] = A_Q_PL
        parameters.push_back(params[4]);            //parameters[7] = B_Q_PL
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMainSAL(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_SAL4(double             const& x,
                                                              double              const& Q,
                                                              std::vector<double> const& params,
                                                              bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(0.3);                  //parameters[0] = K_S
        parameters.push_back(params[0]);            //parameters[1] = B_G_HAD
        parameters.push_back(3.);                   //parameters[2] = C_G_HAD
        parameters.push_back(params[1]);            //parameters[3] = A_Q_HAD
        parameters.push_back(params[2]);            //parameters[4] = B_Q_HAD
        parameters.push_back(1.);                   //parameters[5] = C_Q_HAD
        parameters.push_back(params[3]);            //parameters[6] = A_Q_PL
        parameters.push_back(0.);                   //parameters[7] = B_Q_PL
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMainSAL(x, Q, parameters);
    };

std::map<int, double> StructureFunctionsFcn::InitialPDFs_SAL3(double             const& x,
                                                             double              const& Q,
                                                             std::vector<double> const& params,
                                                             bool                const& returnParameters) const
    {
        std::vector<double> parameters;
        parameters.push_back(0.3);                  //parameters[0] = K_S
        parameters.push_back(params[0]);            //parameters[1] = B_G_HAD
        parameters.push_back(3.);                   //parameters[2] = C_G_HAD
        parameters.push_back(params[1]);            //parameters[3] = A_Q_HAD
        parameters.push_back(params[2]);            //parameters[4] = B_Q_HAD
        parameters.push_back(1.);                   //parameters[5] = C_Q_HAD
        parameters.push_back(0.);                   //parameters[6] = A_Q_PL
        parameters.push_back(0.);                   //parameters[7] = B_Q_PL
        
        if (returnParameters)
        {
            std::map<int, double> parametersMap;

            for (int i=0; i<parameters.size(); i++)
                parametersMap.insert(std::pair<int, double>( i, parameters[i]));

            return parametersMap;
        } 
        else
            return InitialPDFsMainSAL(x, Q, parameters);
    };