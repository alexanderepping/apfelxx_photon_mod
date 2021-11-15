/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include <string>
#include <map>
#include <vector>
#include <functional>


///////////////////////////////////////
// commonly changed
///////////////////////////////////////
/**
 * @brief enumerating the different Initial PDFs
 */
enum enumInitialPDFs { INITIALPDFS_9GDU_00, INITIALPDFS_9GDU_01, INITIALPDFS_3G, INITIALPDFS_2G };

/**
 * @brief Defining the name of the used InitialPDFs. see enumInitialPDFs
 */
const int usedInitialPDFs = INITIALPDFS_9GDU_00;

/**
 * @brief initial parameters for the PDFs
 */
const std::map<int, std::vector<double>> initialParams = {{INITIALPDFS_9GDU_00, {1., 1., 1., 1., 1., 1., 1., 1., 1.}},
                                                          {INITIALPDFS_9GDU_01, {1., 1., 1., 1., 1., 1., 1., 1., 1.}},
                                                          {INITIALPDFS_3G,      {1., 1., 1.}},
                                                          {INITIALPDFS_2G,      {1., 1.}}};

/**
 * @brief initial errors for the parameters for the PDFs
 */
const std::map<int, std::vector<double>> initialParamsErrors = {{INITIALPDFS_9GDU_00, {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_9GDU_01, {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_3G,      {0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_2G,      {0.1, 0.1}}};

/**
 * @brief lower bounds for the initial parameters for the PDFs
 */
const std::map<int, std::vector<double>> initialParamsLBounds = {{INITIALPDFS_9GDU_01, {0., -1., 0., 0., -1., 0., 0., -1., 0.}}};

/**
 * @brief names of the parameters
 */
const std::map<int, std::vector<std::string>> initialParamsNames = {{INITIALPDFS_9GDU_00, {"AN_g1", "A_g1", "B_g1", "AN_d1", "A_d1", "B_d1", "AN_u1", "A_u1", "B_u1"}},
                                                                    {INITIALPDFS_9GDU_00, {"K_s1", "A_g1", "B_g1", "AN_d1", "A_d1", "B_d1", "AN_u1", "A_u1", "B_u1"}},
                                                                    {INITIALPDFS_3G,      {"AN_g1", "A_g1", "B_g1"}},
                                                                    {INITIALPDFS_2G,      {"A_g1", "B_g1"}}};



///////////////////////////////////////
// choose which data is used
///////////////////////////////////////

/**
 * @brief Names of the included experimental data points (see experimentalData.h)
 */
const std::vector<std::string> IncludedExperimentalData = {"ALEPH"};

/**
 * @brief change which LHAPDF data set is used 
 */
#define GRVCustomSetLO



///////////////////////////////////////
// not changed that often
///////////////////////////////////////

#ifdef GRVCustomSetLO
/**
 * @name GRV LO values
 * Some values, usually defined in the LHAPDF .info file, on the GRV LO PDFs.
 */
///@{
/// @brief Name of the used LHAPDF data set
const std::string NameLHAPDFSet = "GRVCustomSetLO";

/// @brief perturbation order of Apfel, not including pointlike contributions 
const int    pto          = 0;
/// @brief reference for energy; mass of the z-boson 
const double Qref         = 91.1876;
/// @brief reference for the strong coupling constant at Qref 
const double asref        = 0.128;
/// @brief mass of the charm quark 
const double mc           = 1.5;
/// @brief mass of the bottom quark 
const double mb           = 4.5;
/// @brief mass of the top quark 
const double mt           = 100;
/// @brief initial energy; lowest energy 
const double Qin          = 1.295000e+00;
///@}
#endif //GRVCustomSetLO

#ifdef GRVCustomSetHO
/**
 * @name GRV HO values
 * Some values, usually defined in the LHAPDF .info file, on the GRV HO PDFs.
 */
///@{
/// @brief Name of the used LHAPDF data set
const std::string NameLHAPDFSet = "GRVCustomSetHO";

/// @brief perturbation order of Apfel, not including pointlike contributions 
const int    pto          = 1;
/// @brief reference for energy; mass of the z-boson 
const double Qref         = 91.1876;
/// @brief reference for the strong coupling constant at Qref 
const double asref        = 0.10902561127771493;
/// @brief mass of the charm quark 
const double mc           = 1.5;
/// @brief mass of the bottom quark 
const double mb           = 4.5;
/// @brief mass of the top quark 
const double mt           = 100;
/// @brief initial energy; lowest energy 
const double Qin          = 1.295000e+00;
///@}
#endif //GRVCustomSetHO