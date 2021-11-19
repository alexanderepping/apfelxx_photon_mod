/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include "LHAPDF/LHAPDF.h"

#include "/usr/local/include/minuit-cpp/FCNBase.hh"

#include "configMinuit.h"

#include <string>
#include <map>
#include <vector>

/**
 * @brief Class derived from the minuit class FCNBase.
 * Calculates the StructureFunctions and returns chi2.
 * Multiple different initial functions for the PDFs 
 * are also included; further can be added. 
 * They can be chosen and configured in configMinuit.h.
 */
class StructureFunctionsFcn : public MinuitCpp::FCNBase 
{
public:
    /**
     * @name Constructors and Destructors
     * List of constructors and destructors.
     */
    ///@{
    /**
     * @brief The StructureFunctionsFcn constructor for the experimentalData map
     * found in experimentalData.h
     * 
     * @param experimentalData: the experimental data in the format seen in experimentalData.h
     * @param ErrorDef: default value = 1.
     */
    StructureFunctionsFcn(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                          std::string const& NameLHAPDFSet,
                          double      const& ErrorDef = 1. );

    /**
     * @brief The StructureFunctionsFcn constructor for the seperate input of experimental data vectors.
     * 
     * @param Energies: vector of Energies from different experiments
     * @param xData: vector of x data from different experiments
     * @param xError: vector of x errors for the xData vector
     * @param F2Gamma: vector of F2Gamma values from different experiments; F2Gamma[i] = F2Gamma(xData[i], Energies[i])
     * @param yError: vector of y errors, the errors of the F2Gamma vector
     * @param ErrorDef: default value = 1.
     * 
     */
    StructureFunctionsFcn(std::vector<double> const& Energies,
                          std::vector<double> const& xData,
                          std::vector<double> const& xError,
                          std::vector<double> const& F2Gamma,
                          std::vector<double> const& yError,
                          std::string         const& NameLHAPDFSet,
                          double              const& ErrorDef = 1. );

    /**
     * @brief The StructureFunctionsFcn destructor.
     */
    ~StructureFunctionsFcn() {}
    ///@}


    /**
     * @brief virtual function to implement the operator ()
     * 
     * @param params: parameters for the initial PDFs
     * 
     * @return returns the chi2 of the initial PDF with the given parameters
     */
    virtual double operator()(std::vector<double> const& params) const;

    /**
     * @brief virtual function to return _ErrorDef
     */
    virtual double Up() const {return _ErrorDef;}

    /**
     * @brief function to call the correct InitialPDF 
     * 
     * @param x
     * @param Q
     * @param params: vector with parameters 
     * @param dist: LHAPDF data set
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMain0)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs(double              const& x,
                                      double              const& Q,
                                      std::vector<double> const& params,
                                      LHAPDF::PDF*               dist,
                                      bool                const& returnParameters = false) const;


    /**
     * @brief function to combine the data from several experiments into only one vector
     * 
     * @param experimentalData: the experimental data in the format seen in experimentalData.h
     * @param dataName: "Name" of the data whose vectors should be combined (e.g. "xData" or "Energies")
     * 
     * @return vector with combined data
     */
    std::vector<double> combineData(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                    std::string                                                       const& dataName) const;

    /** 
     * @brief calculate the beta function = int^1_0 t^(x-1) * (1-t)^(y-1) dt
     * 
     * @param x
     * @param y
     */
    double betaFunction(double const & x, double const& y) const;

    /**
     * @brief Calculate the normalization constant AN_g1 for InitialPDFsMain0 using the momentum sum rule
     * 
     * @param params: vector of parameters
     * @param totMom: total Momentum (lhs of the M.S.R.); default is totalMomentum (see configMinuit.h)
     * 
     * @return AN_g1
     * 
     */
double MomentumSumRule(std::vector<double> const& params,
                       double              const& totMom = totalMomentum) const;


    /**
     * @name get functions
     * List of functions to return the different member data
     */
    ///@{
    /** @brief return vector of Energies */
    std::vector<double> Energies()  const {return _Energies;}
    /** @brief return vector of x data */
    std::vector<double> xData()     const {return _xData;}
    /** @brief return vector of x errors */
    std::vector<double> xError()    const {return _xError;}
    /** @brief return vector of F2Gamma values */
    std::vector<double> F2Gamma()   const {return _F2Gamma;}
    /** @brief return vector of y errors */
    std::vector<double> yError()    const {return _yError;}
    /** @brief return Name of the used LHAPDF Set */
    std::string NameLHAPDFSet()     const {return _NameLHAPDFSet;}
    ///@}


    /** @brief set value of ErrorDef  */
    void setErrorDef(double def) {_ErrorDef = def;}


    /**
     * @name Initial PDFs
     */
    ///@{
    /**
     * @brief Main function to calculate the IntialPDFs of the form 
     * 
     * @param x
     * @param Q
     * @param params: vector with 9 parameters: 
     *        K_s1 (0), A_g1 (1), B_g1 (2), AN_d1 (3), A_d1 (4), B_d1 (5), AN_u1 (6), A_u1 (7), B_u1 (8)
     * @param outputAN_g1: bool to decide if AN_g1 should be output to terminal. Default false. 
     */
    std::map<int, double> InitialPDFsMain0(double                const& x,
                                           double                const& Q,
                                           std::vector<double>   const& params,
                                           bool                  const& outputAN_g1 = false) const;

    /**
     * @brief InitialPDFs with nine params; 3 for up- and down-quark, 2 for gluon and 1 for strange-quark
     * Uses an * x**a * (1-x)**b for up and down, 0 for charm, bottom and top, an * x**a * (1-x)**b for gluon,
     * where an is calculated using the momentum sum rule and K/2 * (x*u + x*d) for the strange quark
     * 
     * @param x
     * @param Q
     * @param params: vector with 9 parameters: 
     *        K_s1 (0), A_g1 (1), B_g1 (2), AN_d1 (3), A_d1 (4), B_d1 (5), AN_u1 (6), A_u1 (7), B_u1 (8)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMain0)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_9gdus  (double                const& x,
                                              double                const& Q,
                                              std::vector<double>   const& params,
                                              bool                  const& returnParameters = false) const;

    /**
     * @brief InitialPDFs with eight params; 3 for up- and down-quark and 2 for gluon
     * Uses an * x**a * (1-x)**b for up and down, 0 for charm, bottom and top, an * x**a * (1-x)**b for gluon,
     * where an is calculated using the momentum sum rule
     * @param x
     * @param Q
     * @param params: vector with 9 parameters: 
     *        A_g1 (0), B_g1 (1), AN_d1 (2), A_d1 (3), B_d1 (4), AN_u1 (5), A_u1 (6), B_u1 (7)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMain0)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_8gdu  (double                const& x,
                                             double                const& Q,
                                             std::vector<double>   const& params,
                                             bool                  const& returnParameters = false) const;

    /**
     * @brief InitialPDFs with eight params; 3 for up- and down-quark and 2 for gluon
     * Uses an * x**a * (1-x)**b for up and down, 0 for charm, bottom and top, an * x**a * (1-x)**b for gluon,
     * where an is calculated using the momentum sum rule
     * @param x
     * @param Q
     * @param params: vector with 9 parameters: 
     *        A_g1 (0), B_g1 (1), AN_q1 (2), A_q1 (3), B_q1 (4)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMain0)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_5gq  (double                const& x,
                                            double                const& Q,
                                            std::vector<double>   const& params,
                                            bool                  const& returnParameters = false) const;

    /**
     * @brief InitialPDFs with nine params; 3 for gluon, up- and down-quark
     * Uses 1 / alpha_em * an * x**a * (1-x)**b for gluon, up and down; the GRV PDFs for the rest.
     * 
     * @param x
     * @param Q
     * @param params: vector with 9 parameters: 
     *        AN_g1 (0), A_g1 (1), B_g1 (2), AN_d1 (3), A_d1 (4), B_d1 (5), AN_u1 (6), A_u1 (7), B_u1 (8)
     * @param dist: LHAPDF data set
     */
    std::map<int, double> InitialPDFs_9gdu  (double                const& x,
                                             double                const& Q,
                                             std::vector<double>   const& params,
                                             LHAPDF::PDF*                 dist) const;

    /**
     * @brief InitialPDFs with two params; 3 for gluon
     * Uses 1 / alpha_em * an * x**a * (1-x)**b for gluon
     * 
     * @param x
     * @param Q
     * @param params: vector with 2 parameters: 
     *        AN_g1 (0), A_g1 (1), B_g1 (2)
     * @param dist: LHAPDF data set
     */
    std::map<int, double> InitialPDFs_3g    (double                const& x,
                                             double                const& Q,
                                             std::vector<double>   const& params,
                                             LHAPDF::PDF*                 dist) const;

    /**
     * @brief InitialPDFs with two params; 2 for gluon
     * Uses 1 / alpha_em * x**a * (1-x)**b for gluon
     * 
     * @param x
     * @param Q
     * @param params: vector with 2 parameters: 
     *        A_g1 (0), B_g1 (1)
     * @param dist: LHAPDF data set
     */
    std::map<int, double> InitialPDFs_2g    (double                const& x,
                                             double                const& Q,
                                             std::vector<double>   const& params,
                                             LHAPDF::PDF*                 dist) const;
    ///@}


private:
    std::vector<double> _Energies;
    std::vector<double> _xData;
    std::vector<double> _xError;
    std::vector<double> _F2Gamma;
    std::vector<double> _yError;
    double              _ErrorDef;
    std::string         _NameLHAPDFSet;
    LHAPDF::PDF*        _LHAPDFSet;
};