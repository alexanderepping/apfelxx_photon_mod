
/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "DeltaChi2.h"
#include "configMinuit.h"
#include "HelperFunctions.h"

#include <boost/math/special_functions/gamma.hpp>
#include <vector>

double XiP(int const& p, 
           int const& N)
{
    return boost::math::gamma_p_inv(N/2., p/100.);
}

double Xi90Rescaled(int    const& Nk,
                    double const& Chi2k0)
{
    return XiP(90, Nk) / XiP(50, Nk) * Chi2k0;
}

std::vector<double> YToZTilde(std::vector<double>                 const& y,
                              Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian)
{
    std::vector<double> zTilde(NumberOfFreeParams, 0);

    for (int i=0; i<NumberOfFreeParams; i++)
        for (int j=0; j<NumberOfFreeParams; j++)
            // z_i^Tilde = sqrt(lambda_i) * Sum_j [ y_j * V_j^(i) ]
            zTilde[i] += std::sqrt(std::real(EigenSolverHessian.eigenvalues()[i])) * y[j] * std::real(EigenSolverHessian.eigenvectors()(j,i));

    return zTilde;
}

std::vector<double> ZTildeToY(std::vector<double>                 const& zTilde,
                              Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian)
{
    std::vector<double> y(NumberOfFreeParams, 0);

    for (int i=0; i<NumberOfFreeParams; i++)
        for (int j=0; j<NumberOfFreeParams; j++)
            // y_i = Sum_j [ 1/sqrt(lambda_j) * V_i^(j) * z_i^Tilde ]
            y[i] += std::real(EigenSolverHessian.eigenvectors()(i,j)) * zTilde[j] / std::sqrt(std::real(EigenSolverHessian.eigenvalues()[j]));

    return y;
}

std::map<std::string, double> CalculateChi2k(StructureFunctionsFcn           const& StructureFunctions,
                                         Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                                         std::vector<double>                 const& finalParams,
                                         double                              const& deltaZ,
                                         int                                 const& i) 
{
        // make the zTilde vectors, but there is only a variation in one ziTilde direction
        // the variation deltaZ is in positive and negative direction
        std::vector<double> zTilde(NumberOfFreeParams, 0.);
        zTilde[i]  = deltaZ;

        // calculate the {y_i}+-
        std::vector<double> yParams  = ZTildeToY(zTilde, EigenSolverHessian);

        // define the vectors, in which the {a_i}^(+-)
        std::vector<double> aParams(NumberOfFreeParams, 0);

        // a_i^(+-) = a_i^0 + y_i^(+-)
        for (int j=0; j<NumberOfFreeParams; j++)
            aParams[j]  = finalParams[j] + yParams[j];
        
        DebugVectorDoubles(aParams, true, true, 1); //debug

        // calculating the chi2k for the plus and minus parameters
        std::map<std::string, double> chi2k  = StructureFunctions.Chi2PerExperiment(aParams);

        return chi2k;
}

double FindZikPlusMinus(std::vector<double> const& chi2kData,
                        std::vector<double> const& zikData,
                        double              const& Xi90Rescaled) 
{
    int j = 0; // the index of the biggest chi2k, which is still smaller than Xi90Rescaled

    for(int i=0; i<chi2kData.size(); i++)
        if (chi2kData[i] >= chi2kData[j] && chi2kData[i] < Xi90Rescaled)
            j = i;
    
    DebugString("chi2k: " + std::to_string(chi2kData[j]) + ", Xi90Rescaled: " + std::to_string(Xi90Rescaled),false); //debug
    DebugString(", z: " + std::to_string(zikData[j]) + ", Index: " + std::to_string(j) + "/" + std::to_string(chi2kData.size() -1), true, false); //debug

    return zikData[j];
}

double CalculateZiPlusMinus(StructureFunctionsFcn               const& StructureFunctions,
                            Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                            int                                 const& i,
                            std::vector<double>                 const& finalParams,
                            std::map<std::string, double>       const& Xi90RescaledMap,
                            std::string                         const& sign,
                            double                              const& deltaZstepsize)
{
    double precision = 0.000001;
    double deltaZ = 0.;
    double stepsize = deltaZstepsize;

    std::vector<std::string> DataSetK = StructureFunctions.IncludedExpData();
    bool FoundZiPlusMinus = false;

    while(stepsize >= precision)
    {
        std::vector<std::string> DataSetKTemp;
        deltaZ += stepsize;

        while(not FoundZiPlusMinus)
        {
            std::map<std::string, double> chi2k;

            DebugString("debug: before CalculateChi2k, deltaZ = " + std::to_string(deltaZ)); //debug

            if (sign == "plus") //debug
            {
                DebugString("", true, false, 1); //debug
                DebugString("deltaZ_" + std::to_string(i) + "^+ = " + std::to_string(deltaZ) + ", stepsize = " + std::to_string(stepsize), true, true, 1); //debug
            }
            else if (sign == "minus") //debug
            {
                DebugString("", true, false, 1); //debug
                DebugString("deltaZ_" + std::to_string(i) + "^- = " + std::to_string(deltaZ) + ", stepsize = " + std::to_string(stepsize), true, true, 1); //debug
            }

            if (sign == "plus")
                chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, +deltaZ, i);
            else if (sign == "minus")
                chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, -deltaZ, i);

            for (std::string DataSet : DataSetK)
            {
                    DebugString("  " + DataSet + ": " + "chi2k: " + std::to_string(chi2k.at(DataSet)) + ", Xi: " + std::to_string(Xi90RescaledMap.at(DataSet)), true, true, 1); //debug

                if (chi2k.at(DataSet) >= Xi90RescaledMap.at(DataSet))
                {
                    DataSetKTemp.push_back(DataSet);
                    FoundZiPlusMinus = true;
                    
                    DebugString("      Found chi2k > Xi90Rescaled : " + DataSet, true, true, 1); //debug
                }
            }

            if (not FoundZiPlusMinus)
                deltaZ += stepsize;
        }

        deltaZ -= stepsize;
        stepsize = stepsize/10.;
        FoundZiPlusMinus = false;
        if (DataSetKTemp.size() > 0)
            DataSetK = DataSetKTemp;
    }

    return deltaZ;
}