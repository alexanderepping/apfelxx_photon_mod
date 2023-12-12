/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

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
     * @param IncludedExpData: the experimental data, that is included from experimentalData.h
     * @param ErrorDef: default value = 1.
     */
    StructureFunctionsFcn(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                          std::vector<std::string>                                          const& IncludedExpData = IncludedExperimentalData,
                          double      const& ErrorDef = 1. );

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
     * @return returns the chi2 of the initial PDF with the given parameters, using Chi2PerExperiment()
     */
    virtual double operator()(std::vector<double> const& params) const;

    /**
     * @brief virtual function to return _ErrorDef
     */
    virtual double Up() const {return _ErrorDef;}

    /**
     * @brief function used in operator() and can be used on its own aswell
     * 
     * @param params: parameters for the initial PDFs
     * 
     * @return returns the Chi2PerExperiment of the initial PDF with the given parameters
     */
    std::map<std::string, double> Chi2PerExperiment(std::vector<double> const& params) const;

    /**
     * @brief function to call the correct InitialPDF 
     * 
     * @param x
     * @param Q
     * @param params: vector with parameters 
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMain0)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs(double              const& x,
                                      double              const& Q,
                                      std::vector<double> const& params,
                                      bool                const& returnParameters = false) const;


    /**
     * @brief function to combine the data from several experiments into only one vector
     * 
     * @param experimentalData: the experimental data in the format seen in experimentalData.h
     * @param IncludedExpData: the experimental data, that is included from experimentalData.h
     * @param dataName: "Name" of the data whose vectors should be combined (e.g. "xData" or "Q2Data")
     * 
     * @return vector with combined data
     */
    std::vector<double> combineData(std::map<std::string, std::map<std::string, std::vector<double>>> const& experimentalData,
                                    std::vector<std::string>                                          const& IncludedExpData,
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
    double MomentumSumRule0(std::vector<double> const& params,
                            double              const& totMom = totalMomentum) const;

    /**
     * @brief Calculate the normalization constant A_G_HAD for InitialPDFsMainSAL using the FG momentum sum rule
     * (see e.g. SAL, eq. 21 (https://arxiv.org/abs/hep-ph/0504003v2))
     * 
     * @param params: vector of parameters
     * @param Q
     * 
     * @return A_G_HAD
     * 
     */
    double MomentumSumRuleSAL(std::vector<double> const& params,
                            double              const& Q) const;


    /**
     * @name get functions
     * List of functions to return the different member data and data that is similar to that
     */
    ///@{
    /** @brief return experimentalData map */
    std::map<std::string, std::map<std::string, std::vector<double>>> ExperimentalData() const {return _ExperimentalData;}
    /** @brief return included experimental Data */
    std::vector<std::string> IncludedExpData() const {return _IncludedExpData;}
    /** @brief return vector of squared Energies */
    std::vector<double> Q2Data() const {return combineData(_ExperimentalData, _IncludedExpData, "Q2Data");}
    /** @brief return vector of x data */
    std::vector<double> xData() const {return combineData(_ExperimentalData, _IncludedExpData, "xData");}
    /** @brief return vector of F2Gamma values */
    std::vector<double> F2Gamma() const {return combineData(_ExperimentalData, _IncludedExpData, "F2Gamma");}
    /** @brief return vector of y errors */
    std::vector<double> F2GammaErr() const {return combineData(_ExperimentalData, _IncludedExpData, "F2GammaErr");}
    ///@}


    /** @brief set value of ErrorDef  */
    void setErrorDef(double def) {_ErrorDef = def;}


    /**
     * @name Initial PDFs
     */
    ///@{
    /**
     * @brief Main function to calculate the IntialPDFs of the form :
     * an * x**a * (1-x)**b for up and down, 0 for charm, bottom and top, an * x**a * (1-x)**b for gluon,
     * where an is calculated using the momentum sum rule and K/2 * (x*u + x*d) for the strange quark
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
     * @brief Main function to calculate the IntialPDFs of the form given in
     * SAL (https://arxiv.org/abs/hep-ph/0504003v2), 
     * where 0.3 (for the strange quark) is replaced by a free parameter
     * 
     * @param x
     * @param Q
     * @param params: vector with 8 parameters: 
     *        K_S (0), B_G_HAD(1), C_G_HAD(2), A_Q_HAD(3), B_Q_HAD(4), C_Q_HAD(5), A_Q_PL(6), B_Q_PL(7)
     * @param outputA_G_HAD: bool to decide if A_G_HAD should be output to terminal. Default false. 
     */
    std::map<int, double> InitialPDFsMainSAL(double                const& x,
                                             double                const& Q,
                                             std::vector<double>   const& params,
                                             bool                  const& outputA_G_HAD = false) const;

    /**
     * @brief InitialPDFs with 9 params; 3 for up- and down-quark, 2 for gluon and 1 for strange-quark
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
     * @brief InitialPDFs with 8 params; like InitialPDFs_9gdus, but the coefficient of strange quark is fixed.
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
     * @brief InitialPDFs with 6 params; like InitialPDFs_9gdus, but u=d.
     * @param x
     * @param Q
     * @param params: vector with 9 parameters: 
     *        K_s1 (0), A_g1 (1), B_g1 (2), AN_q1 (3), A_q1 (4), B_q1 (5)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMain0)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_6gqs  (double                const& x,
                                             double                const& Q,
                                             std::vector<double>   const& params,
                                             bool                  const& returnParameters = false) const;

    /**
     * @brief InitialPDFs with 5 params; like InitialPDFs_9gdus, but coefficient of strange quark is fixed and u=d.
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
     * @brief InitialPDFs with 8 params; has the form of InitialPDFsMainSAL
     * @param x
     * @param Q
     * @param params: vector with 8 parameters: 
     *        K_S (0), B_G_HAD(1), C_G_HAD(2), A_Q_HAD(3), B_Q_HAD(4), C_Q_HAD(5), A_Q_PL(6), B_Q_PL(7)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMainSAL)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_SAL8(double              const& x,
                                           double              const& Q,
                                           std::vector<double> const& params,
                                           bool                const& returnParameters = false) const;  

    /**
     * @brief InitialPDFs with 6 params; has the form of InitialPDFsMainSAL, but C_G_HAD is set to 3 and C_Q_HAD to 1
     * @param x
     * @param Q
     * @param params: vector with 6 parameters: 
     *        K_S (0), B_G_HAD(1), A_Q_HAD(2), B_Q_HAD(3), A_Q_PL(4), B_Q_PL(5)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMainSAL)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_SAL6(double              const& x,
                                           double              const& Q,
                                           std::vector<double> const& params,
                                           bool                const& returnParameters = false) const;  

    /**
     * @brief InitialPDFs with 5 params; has the form of InitialPDFsMainSAL, but C_G_HAD is set to 3, C_Q_HAD to 1 and K_S to 0.3.
     * The same form as originaly used by SAL with the exception that A_G_HAD is determined by the momentum sum rule
     * @param x
     * @param Q
     * @param params: vector with 5 parameters: 
     *        B_G_HAD(0), A_Q_HAD(1), B_Q_HAD(2), A_Q_PL(3), B_Q_PL(4)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMainSAL)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_SAL5(double              const& x,
                                           double              const& Q,
                                           std::vector<double> const& params,
                                           bool                const& returnParameters = false) const; 

    /**
     * @brief InitialPDFs with 4 params; has the form of InitialPDFsMainSAL, but C_G_HAD is set to 3, C_Q_HAD to 1, B_Q_PL to 0 and K_S to 0.3.
     * Same as SAL5, but B_Q_PL is set to zero. Same as used by Vadim
     * @param x
     * @param Q
     * @param params: vector with 4 parameters: 
     *        B_G_HAD(0), A_Q_HAD(1), B_Q_HAD(2), A_Q_PL(3)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMainSAL)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_SAL4(double              const& x,
                                           double              const& Q,
                                           std::vector<double> const& params,
                                           bool                const& returnParameters = false) const; 

    /**
     * @brief InitialPDFs with 3 params; has the form of InitialPDFsMainSAL, but C_G_HAD is set to 3, C_Q_HAD to 1, K_S to 0.5
     * and the pointlike parts of the quarks are set to zero. Like SAL5 w/out PL part.
     * @param x
     * @param Q
     * @param params: vector with 3 parameters: 
     *        B_G_HAD(0), A_Q_HAD(1), B_Q_HAD(2)
     * @param returnParameters: Boolean to either return the PDFs (false) 
     * or return the Parameters (in the form accepted by InitialPDFsMainSAL)(true). 
     * Default is false.
     */
    std::map<int, double> InitialPDFs_SAL3(double             const& x,
                                           double              const& Q,
                                           std::vector<double> const& params,
                                           bool                const& returnParameters = false) const;                                      
    ///@}


private:
    std::map<std::string, std::map<std::string, std::vector<double>>> _ExperimentalData;
    std::vector<std::string> _IncludedExpData;
    double                   _ErrorDef;
};