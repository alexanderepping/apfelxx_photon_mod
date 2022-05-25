/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"

#include <Eigen/Eigenvalues>
#include <vector>
#include <functional>



double SecondDerivative(StructureFunctionsFcn const& function,
                        std::vector<double>   const& a0,    
                        int                   const& parameter1, 
                        int                   const& parameter2, 
                        double                const& h) 
{
    std::function<std::vector<double>(std::vector<double> const&, int const&, int const&)> a1 = [&] (std::vector<double> const& a0, int const& parameter, int const& k) -> std::vector<double>
    {
        std::vector<double> a1 = a0;
        a1[parameter] += k*h;
        return a1;
    };

    std::function<double(std::vector<double> const&)> FirstDerivative = [&] (std::vector<double> const& a0) -> double 
        {
            return ((function(a1(a0, parameter1, 1)) - function(a1(a0, parameter1, -1))) 
              + 2 * (function(a1(a0, parameter1, 2)) - function(a1(a0, parameter1, -2))) 
              + 3 * (function(a1(a0, parameter1, 3)) - function(a1(a0, parameter1, -3)))) / (28 * h);
        };

    return ((FirstDerivative(a1(a0, parameter2, 1)) - FirstDerivative(a1(a0, parameter2, -1))) 
      + 2 * (FirstDerivative(a1(a0, parameter2, 2)) - FirstDerivative(a1(a0, parameter2, -2))) 
      + 3 * (FirstDerivative(a1(a0, parameter2, 3)) - FirstDerivative(a1(a0, parameter2, -3)))) / (28 * h);
};


Eigen::MatrixXd CalculateHessian(StructureFunctionsFcn const& function,
                                 std::vector<double>   const& a0,
                                 double                const& h) 
{
    const int numberParams = a0.size();
    Eigen::MatrixXd Hessian(numberParams, numberParams);
    // Eigen::MatrixXd Hessian();
    // Eigen::Matrix<double, numberParams, numberParams> Hessian = Eigen::Matrix<double, numberParams, numberParams>::Zero();

    for (int i=0; i<a0.size(); i++)
        for (int j=i; j<a0.size(); j++)
            Hessian(i, j) = 0.5 * SecondDerivative(function, a0, i, j, h);

    for (int i=0; i<a0.size(); i++)
        for (int j=0; j<i; j++)
            Hessian(i, j) = Hessian(j, i);

    return Hessian;
};


