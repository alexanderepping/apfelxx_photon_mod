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
 * @brief Defining the name of the used InitialPDFs.
 * - testInitialPDFs            9: some InitialPDFs to test if everything works
 * - minimaltestInitialPDFs     2: testInitialPDFs with less params
 */
#define minimaltestInitialPDFs

/**
 * @brief initial parameters for the PDFs
 */
#ifdef testInitialPDFs
const std::vector<double> initialParams = {1., 1., 1., 1., 1., 1., 1., 1., 1.};
#endif
#ifdef minimaltestInitialPDFs
const std::vector<double> initialParams = {1., 1.};
#endif

/**
 * @brief initial errors for the parameters for the PDFs
 */
#ifdef testInitialPDFs
const std::vector<double> initialErrorParams = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
#endif
#ifdef minimaltestInitialPDFs
const std::vector<double> initialErrorParams = {0.1, 0.1};
#endif
/**
 * @brief names of the parameters
 */
#ifdef testInitialPDFs
const std::vector<std::string> ParamsNames = {"AN_glu1", "A_glu1", "B_glu1", "AN_dbar1", "A_dbar1", "B_dbar1", "AN_ubar1", "A_ubar1", "B_ubar1"};
#endif
#ifdef minimaltestInitialPDFs
const std::vector<std::string> ParamsNames = {"A_glu1", "B_glu1"};
#endif




///////////////////////////////////////
// not changed that often
///////////////////////////////////////
/**
 * @brief Names of the included experimental data points (see experimentalData.h)
 */
const std::vector<std::string> IncludedExperimentalData = {"ALEPH", "AMY"};

/**
 * @brief Name of the LHAPDF-Set used by exampleInitialPDF and as a comparison in StructureFunctionsFcn::operator() 
 */
const std::string NameLHAPDFSet = "GRVCustomSetLO";
#define GRVCustomSetLO

#ifdef GRVCustomSetLO
/**
 * @name GRV LO values
 * Some values, usually defined in the LHAPDF .info file, on the GRV LO PDFs.
 */
///@{
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