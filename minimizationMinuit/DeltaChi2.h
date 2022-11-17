
/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#pragma once

#include "/usr/local/include/minuit-cpp/FCNBase.hh"

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"

#include <Eigen/Eigenvalues>
#include <boost/math/special_functions/gamma.hpp>
#include <vector>

/**
 * @brief function to calculate Xi_p from nCTEQ15 (A2); inverse, incomplete gamma function of N/2 and p%
 * 
 * @param p: percentage p%, that fit with Chi2 > Xi_p genuinely describes given set of data
 * @param N: number of data points
 * 
 * @return Xi_p
 */
double XiP(int const& p, 
           int const& N);

/**
 * @brief function to calculate Xi90, rescaled, from nCTEQ15 (A3)
 * 
 * @param Nk: number of data points for experiment k
 * @param Chi2k0: minimum Chi2 for experiment k
 * 
 * @return rescaled Xi90
 */
double Xi90Rescaled(int    const& Nk,
                   double const& Chi2k0);

/**
 * @brief function to transform parameter set y to parameter set zTilde
 * 
 * @param y: parameter set y w/ y_i = a_i - a_i^0
 * @param EigenSolverHessian: object to calculate the eigenvalues and eigenvectors of the Hessian
 * 
 * @return vector zTilde: {z_i^Tilde}
 */
std::vector<double> YToZTilde(std::vector<double>                 const& y,
                              Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian);

/**
 * @brief function to transform parameter set zTilde back to parameter set y
 * 
 * @param zTilde: parameter set zTilde 
 * @param EigenSolverHessian: object to calculate the eigenvalues and eigenvectors of the Hessian
 * 
 * @return vector y: {y_i} w/ y_i = a_i - a_i^0
 */
std::vector<double> ZTildeToY(std::vector<double>                 const& zTilde,
                              Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian);

std::map<std::string, double> CalculateChi2k(StructureFunctionsFcn           const& StructureFunctions,
                                         Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                                         std::vector<double>                 const& finalParams,
                                         double                              const& deltaZ,
                                         int                                 const& i);

double FindZikPlusMinus(std::vector<double> const& chi2kData,
                        std::vector<double> const& zikData,
                        double              const& Xi90Rescaled);

double CalculateZiPlusMinus(StructureFunctionsFcn               const& StructureFunctions,
                            Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                            int                                 const& i,
                            std::vector<double>                 const& finalParams,
                            std::map<std::string, double>       const& Xi90RescaledMap,
                            std::string                         const& sign,
                            double                              const& deltaZstepsize);