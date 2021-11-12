/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include "LHAPDF/LHAPDF.h"

#include "/usr/local/include/minuit-cpp/FCNBase.hh"

#include <string>
#include <map>
#include <vector>

/**
 * @brief Class derived from the minuit class FCNBase.
 * Calculates the StructureFunctions and returns chi2.
 * Multiple different initial functions for the PDFs 
 * are possible and can be configired inside 
 * InitialPDFs.h and configMinuit.h.
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