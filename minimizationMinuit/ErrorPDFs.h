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
#include <vector>


/**
 * @brief function to calculate the second derivative of a StructureFunctionsFcn. Used to calculate Hessian Matrix.
 * 
 * @param function: function which should be dervived. In this case it gives the chi2 of the structure functions.
 * @param a0: parameters around which the function should be derived.
 * @param parameter1: index of the parameter, whith respect to which the function should be derived first.
 * @param parameter2: index of the parameter, whith respect to which the first derivative of the function should be derived.
 * @param h: step width of the derivative
 * 
 * @return value of the second derivative of the function at a0 in direction of the two parameters
 */
std::map<std::string, double> SecondDerivativeMap(StructureFunctionsFcn const& function,
                                               std::vector<double>   const& a0,    
                                               int                   const& parameter1, 
                                               int                   const& parameter2, 
                                               double                const& h=0.005);

/**
 * @brief function to calculate the second derivative of a StructureFunctionsFcn. Used to calculate Hessian Matrix.
 * 
 * @param function: function which should be dervived. In this case it gives the chi2 of the structure functions.
 * @param a0: parameters around which the function should be derived.
 * @param parameter1: index of the parameter, whith respect to which the function should be derived first.
 * @param parameter2: index of the parameter, whith respect to which the first derivative of the function should be derived.
 * @param h: step width of the derivative
 * 
 * @return value of the second derivative of the function at a0 in direction of the two parameters
 */
double SecondDerivative(StructureFunctionsFcn const& function,
                        std::vector<double>   const& a0,    
                        int                   const& parameter1, 
                        int                   const& parameter2, 
                        double                const& h=0.005);


/**
 * @brief function to calculate the Hessian Matrix of a function.
 * 
 * @param function: function for which to calculate the hessian matrix. In this case it gives the chi2 of the structure functions.
 * @param a0: parameters around which the function should be derived.
 * @param h: step width of the derivative
 * 
 * @return Hessian Matrix
 */
std::map<std::string, Eigen::MatrixXd> CalculateHessianMap(StructureFunctionsFcn const& function,
                                                        std::vector<double>   const& a0,
                                                        double                const& h=0.005);

/**
 * @brief function to calculate the Hessian Matrix of a function.
 * 
 * @param function: function for which to calculate the hessian matrix. In this case it gives the chi2 of the structure functions.
 * @param a0: parameters around which the function should be derived.
 * @param h: step width of the derivative
 * 
 * @return Hessian Matrix
 */
Eigen::MatrixXd CalculateHessian(StructureFunctionsFcn const& function,
                                 std::vector<double>   const& a0,
                                 double                const& h=0.005);