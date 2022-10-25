/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once 

#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <functional> 


///////////////////////////////////////
// commonly changed
///////////////////////////////////////
/**
 * @brief Defining the name of the used InitialPDFs. see enumInitialPDFs
 */
const std::string nameUsedInitialPDFs = "INITIALPDFS_SAL5";

/**
 * @brief change which perturbation order is used
 */
#define HO

/**
 * @brief change if ErrorPDFs should be calculated
 */
//#define ErrorPDFs

/**
 * @brief path of the file where data should be saved
 */
const std::string outputFile = "/home/alexander/Uni/apfelxx_photon_mod/plottingPython/dataInitialPDFs.txt";



///////////////////////////////////////
// experimental Data selection
///////////////////////////////////////

/**
 * @brief Names of the included experimental data points (see experimentalData.h).
 * All names: "ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2", "OPAL2_dop", "PLUTO", "TASSO", "TOPAZ", "TPC"
 */
//const std::vector<std::string> IncludedExperimentalData = {"ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "PLUTO", "TASSO", "TOPAZ"};
//const std::vector<std::string> IncludedExperimentalData = {"ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2", "PLUTO", "TASSO", "TOPAZ"};
const std::vector<std::string> IncludedExperimentalData = {"ALEPH1", "ALEPH2", "AMY", "DELPHI", "JADE", "L3", "OPAL1", "OPAL2_less", "PLUTO", "TASSO", "TOPAZ"};



///////////////////////////////////////
// not changed that often
///////////////////////////////////////
/**
 * @brief the total momentum / result of MomentumSumRule0
 */
const double totalMomentum = 1.;


/**
 * @brief change the derivative used in the Calculation of the Hessian.
 * The more function calls, the lower the numerical noise of derivative but also takes longer.
 *      - StandardCentralDifferences: two function calls per derivative
 *      - SevenPointLowNoise: six function calls per derivative
 */
#define StandardCentralDifferences

/**
 * @brief Delta Chi^2, used to calculate the error PDFs
 */
const double DeltaChi2 = 1.;

#ifdef LO
/**
 * @name LO parameters
 */
///@{
/// @brief perturbation order of Apfel
const int    pto          = 0;
/// @brief reference for energy; mass of the z-boson 
const double Qref         = 91.188;
/// @brief reference for the strong coupling constant at Qref 
const double asref        = 0.1179973;
/// @brief mass of the charm quark 
const double mc           = 1.3;
/// @brief mass of the bottom quark 
const double mb           = 4.5;
/// @brief mass of the top quark 
const double mt           = 174;
/// @brief initial energy; lowest energy 
const double Qin          = 1.3;
///@}
#endif //LO

#ifdef HO
/**
 * @name HO parameters
 */
///@{
/// @brief perturbation order of Apfel
const int    pto          = 1;
/// @brief reference for energy; mass of the z-boson 
const double Qref         = 91.188;
/// @brief reference for the strong coupling constant at Qref 
const double asref        = 0.1179973;
/// @brief mass of the charm quark 
const double mc           = 1.3;
/// @brief mass of the bottom quark 
const double mb           = 4.5;
/// @brief mass of the top quark 
const double mt           = 174;
/// @brief initial energy; lowest energy 
const double Qin          = 1.3;
///@}
#endif //HO


///////////////////////////////////////
// initialPDF parameter settings
///////////////////////////////////////
/**
 * @brief enumerating the different Initial PDFs
 */
enum enumInitialPDFs { // InitialPDFs related to InitialPDFsMain0
                       INITIALPDFS_9GDUS, 
                       INITIALPDFS_8GDU, 
                       INITIALPDFS_6GQS, 
                       INITIALPDFS_5GQ, 

                       // InitialPDFs related to InitialPDFsMainSAL
                       INITIALPDFS_SAL8, 
                       INITIALPDFS_SAL6, 
                       INITIALPDFS_SAL5,        // calculated parameters can be compared to the SAL input PDFs
                       INITIALPDFS_SAL4VADIM,   // same as SAL5 w/out B_Q_PL, used by Vadim
                       INITIALPDFS_SAL4,        // same as SAL6 w/out PL part
                       INITIALPDFS_SAL3};       // same as SAL5 w/out PL part

/**
 * @brief names of the InitialPDFs
 */
const std::map<std::string, int> initialPDFsNames = {{"INITIALPDFS_9GDUS",      INITIALPDFS_9GDUS},
                                                     {"INITIALPDFS_8GDU",       INITIALPDFS_8GDU},
                                                     {"INITIALPDFS_6GQS",       INITIALPDFS_6GQS},
                                                     
                                                     {"INITIALPDFS_SAL8",       INITIALPDFS_SAL8},
                                                     {"INITIALPDFS_SAL6",       INITIALPDFS_SAL6},
                                                     {"INITIALPDFS_SAL5",       INITIALPDFS_SAL5},
                                                     {"INITIALPDFS_SAL4VADIM",  INITIALPDFS_SAL4VADIM},
                                                     {"INITIALPDFS_SAL4",       INITIALPDFS_SAL4},
                                                     {"INITIALPDFS_SAL3",       INITIALPDFS_SAL3}};

/**
 * @brief Defining the name of the used InitialPDFs. see enumInitialPDFs
 */
const int usedInitialPDFs = initialPDFsNames.at(nameUsedInitialPDFs);

/**
 * @brief initial parameters for the PDFs
 */
const std::map<int, std::vector<double>> initialParams = {{INITIALPDFS_9GDUS,           {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_8GDU,            {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_6GQS,            {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_5GQ,             {0.5, 0.5, 0.5, 0.5, 0.5}},

                                                          {INITIALPDFS_SAL8,            {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_SAL6,            {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_SAL5,            {0.5, 0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_SAL4VADIM,       {0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_SAL4,            {0.5, 0.5, 0.5, 0.5}},
                                                          {INITIALPDFS_SAL3,            {0.5, 0.5, 0.5}}};

/**
 * @brief initial errors for the parameters for the PDFs
 */
const std::map<int, std::vector<double>> initialParamsErrors = {{INITIALPDFS_9GDUS,             {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_8GDU,              {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_6GQS,              {0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_5GQ,               {0.1, 0.1, 0.1, 0.1, 0.1}},

                                                                {INITIALPDFS_SAL8,              {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_SAL6,              {0.1, 0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_SAL5,              {0.1, 0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_SAL4VADIM,         {0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_SAL4,              {0.1, 0.1, 0.1, 0.1}},
                                                                {INITIALPDFS_SAL3,              {0.1, 0.1, 0.1}}};

/**
 * @brief upper bounds for the initial parameters for the PDFs.
 * bound for K_s1 comes from the fact that s should be smaller than 1/2(u+d).
 * 
 */
        // Main0:  K_s1 (0), A_g1 (1), B_g1 (2), AN_d1 (3), A_d1 (4), B_d1 (5), AN_u1 (6), A_u1 (7), B_u1 (8)
        // SAL:    K_S (0), B_G_HAD(1), C_G_HAD(2), A_Q_HAD(3), B_Q_HAD(4), C_Q_HAD(5), A_Q_PL(6), B_Q_PL(7)
const std::map<int, std::vector<double>> initialParamsUBounds = {{INITIALPDFS_9GDUS,            { 1.,  1.,  5., 40.,  1.,  5., 40.,  1.,  5.}},
                                                                 {INITIALPDFS_8GDU,             { 1.,  5., 40.,  1.,  5., 40.,  1.,  5.}},
                                                                 {INITIALPDFS_6GQS,             { 1.,  1.,  5., 40.,  1.,  5.}},
                                                                 {INITIALPDFS_5GQ,              { 1.,  5., 40.,  1.,  5.}},

                                                                 {INITIALPDFS_SAL8,             { 1.,  1.,  5., 40.,  1.,  5., 40., 40.}},
                                                                 {INITIALPDFS_SAL6,             { 1.,  1., 40.,  1., 40., 40.}},
                                                                 {INITIALPDFS_SAL5,             { 1., 40.,  1., 40., 40.}}, 
                                                                 {INITIALPDFS_SAL4VADIM,        { 1., 40.,  1., 40.}}, 
                                                                 {INITIALPDFS_SAL4,             { 1.,  1., 40.,  1.}},
                                                                 {INITIALPDFS_SAL3,             { 1., 40.,  1.}}};

/**
 * @brief lower bounds for the initial parameters for the PDFs.
 * bounds come from the fact that the gamma function only takes values bigger than zero,
 * the argument of the exponential integral shouldn't go to -inf and 
 * the PDFs should be positive.
 * Also, the B_Q_HAD are limited such that ExpInt and Exp in the MomentumSumRuleSAL don't get too big. 
 */
const std::map<int, std::vector<double>> initialParamsLBounds = {{INITIALPDFS_9GDUS,            { 0., -1.,  0.,  0., -1.,  0.,  0., -1.,  0.}},
                                                                 {INITIALPDFS_8GDU,             {-1.,  0.,  0., -1.,  0.,  0., -1.,  0.}},
                                                                 {INITIALPDFS_6GQS,             { 0., -1.,  0.,  0., -1.,  0.}},
                                                                 {INITIALPDFS_5GQ,              {-1.,  0.,  0., -1.,  0.}},

                                                                 {INITIALPDFS_SAL8,             { 0., -1.,  0.,  0., -1.,  0.,  0.,  0.1}},
                                                                 {INITIALPDFS_SAL6,             { 0., -1.,  0., -1.,  0.,  0.1}},
                                                                 {INITIALPDFS_SAL5,             {-1.,  0., -1.,  0.,  0.1}},
                                                                 {INITIALPDFS_SAL4VADIM,        {-1.,  0., -1.,  0.}},
                                                                 {INITIALPDFS_SAL4,             { 0., -1.,  0., -1.}},
                                                                 {INITIALPDFS_SAL3,             {-1.,  0., -1.}}};

/**
 * @brief names of the parameters
 */
const std::map<int, std::vector<std::string>> initialParamsNames = {{INITIALPDFS_9GDUS,                 {"K_s1", "A_g1", "B_g1", "AN_d1", "A_d1", "B_d1", "AN_u1", "A_u1", "B_u1"}},
                                                                    {INITIALPDFS_8GDU,                  {"A_g1", "B_g1", "AN_d1", "A_d1", "B_d1", "AN_u1", "A_u1", "B_u1"}},
                                                                    {INITIALPDFS_6GQS,                  {"K_s1", "A_g1", "B_g1", "AN_q1", "A_q1", "B_q1"}},
                                                                    {INITIALPDFS_5GQ,                   {"A_g1", "B_g1", "AN_q1", "A_q1", "B_q1"}},

                                                                    {INITIALPDFS_SAL8,                  {"K_S", "B_G_HAD", "C_G_HAD", "A_Q_HAD", "B_Q_HAD", "C_Q_HAD", "A_Q_PL", "B_Q_PL"}},
                                                                    {INITIALPDFS_SAL6,                  {"K_S", "B_G_HAD", "A_Q_HAD", "B_Q_HAD", "A_Q_PL", "B_Q_PL"}},
                                                                    {INITIALPDFS_SAL5,                  {"B_G_HAD", "A_Q_HAD", "B_Q_HAD", "A_Q_PL", "B_Q_PL"}},
                                                                    {INITIALPDFS_SAL4VADIM,             {"B_G_HAD", "A_Q_HAD", "B_Q_HAD", "A_Q_PL"}},
                                                                    {INITIALPDFS_SAL4,                  {"K_S", "B_G_HAD", "A_Q_HAD", "B_Q_HAD"}},
                                                                    {INITIALPDFS_SAL3,                  {"B_G_HAD", "A_Q_HAD", "B_Q_HAD"}}};