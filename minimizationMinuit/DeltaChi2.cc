
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
    return 2. * boost::math::gamma_p_inv(N/2., p/100.);
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

std::vector<double> FindZTildeLimit(int                                 const& j,
                                    std::vector<double>                 const& finalParams,
                                    Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian)
{
    std::vector<double> zTildeMax(NumberOfFreeParams, 0);

    for (int i=0; i<NumberOfFreeParams; i++)
    {
        // calculated from nCTEQ15 (2.24), but only zTilde_j with the specific j (given to this function) is not zero
        zTildeMax[i] = (initialParamsLBounds.at(usedInitialPDFs)[i] - finalParams[i]) * std::sqrt(std::real(EigenSolverHessian.eigenvalues()[j])) / std::real(EigenSolverHessian.eigenvectors()(i,j));

        DebugString("a_" + std::to_string(j) + " = " + std::to_string(initialParamsLBounds.at(usedInitialPDFs)[i]) + " for zTilde_" + std::to_string(i) + " = " + std::to_string(zTildeMax[i])); //debug

        std::vector<double> zTildep(NumberOfFreeParams, 0.);
        std::vector<double> zTildem(NumberOfFreeParams, 0.);
        zTildep[i]  = zTildeMax[j]+1;
        zTildem[i]  = zTildeMax[j]-1;
        std::vector<double> yp = ZTildeToY(zTildep, EigenSolverHessian);
        std::vector<double> ym = ZTildeToY(zTildem, EigenSolverHessian);

        for (int j=0; j<NumberOfFreeParams; j++)
        {
            yp[j] += finalParams[j];
            ym[j] += finalParams[j];
        }

        DebugString("  a_i+ = ", false); //debug
        DebugVectorDoubles(yp, true, false); //debug
        DebugString("  a_i- = ", false); //debug
        DebugVectorDoubles(ym, true, false); //debug
    }

    return zTildeMax;
}

void FindDeltaChi2Limit(std::vector<double>                const& finalParams,
                       Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                       int                                 const& localRequiredVerbosity)
{
    // loop through all the different functions fk-
    for (int k=0; k<NumberOfFreeParams; k++)
    {
        DebugString("k = " + std::to_string(k), true, true, localRequiredVerbosity); //debug
        
        // loop throug all the parameters
        for (int i=0; i<NumberOfFreeParams; i++)
        {
            // calculated from nCTEQ15 (2.32)
            double DeltaChi2Max = std::pow(-(initialParamsLBounds.at(usedInitialPDFs)[i] - finalParams[i]) / std::real(EigenSolverHessian.eigenvectors()(i,k)), 2.) * std::real(EigenSolverHessian.eigenvalues()[k]);
            DebugString("a_" + std::to_string(i) + " = " + std::to_string(initialParamsLBounds.at(usedInitialPDFs)[i]) + " for DeltaChi2" + " < " + std::to_string(DeltaChi2Max), true, true, localRequiredVerbosity); //debug

            std::vector<std::string> tempS = {"-a_i-", "-a_i ", "-a_i+", "+a_i-", "+a_i ", "+a_i+"};
            std::vector<std::vector<double>> temp = {{}, {}, {}, {}, {}, {}};

            for (int j=0; j<NumberOfFreeParams; j++)
            {
                temp[0].push_back(finalParams[j] - std::sqrt((DeltaChi2Max-1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j)));
                temp[1].push_back(finalParams[j] - std::sqrt((DeltaChi2Max+0.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j)));
                temp[2].push_back(finalParams[j] - std::sqrt((DeltaChi2Max+1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j)));
                temp[3].push_back(finalParams[j] + std::sqrt((DeltaChi2Max-1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j)));
                temp[4].push_back(finalParams[j] + std::sqrt((DeltaChi2Max+0.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j)));
                temp[5].push_back(finalParams[j] + std::sqrt((DeltaChi2Max+1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j)));
            }

            for (int j=0; j<NumberOfFreeParams; j++)
            {
                DebugString("    "+tempS[j]+" = ", false, true, localRequiredVerbosity); //debug
                DebugVectorDoubles(temp[j], true, false, localRequiredVerbosity); //debug
            }
        }
    }
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
    double deltaZMax = 0.;
    double stepsize = deltaZstepsize;

    std::vector<std::string> DataSetK = StructureFunctions.IncludedExpData();
    bool FoundZiPlusMinus = false;


    // find the deltaZ, for which the lower bound of one parameter is violated
    DebugString("Finding the ZTildeLimit for i=" + std::to_string(i) + ", " + sign + ": ", true, true, 1); //debug
    std::vector<double> ZTildeLimit = FindZTildeLimit(i, finalParams, EigenSolverHessian);
    for (double zTilde : ZTildeLimit)
        if ((sign == "plus" && zTilde > 0.) || (sign == "minus" && zTilde < 0.))
            if (deltaZMax == 0. || deltaZMax > std::abs(zTilde))
                deltaZMax = std::abs(zTilde);
    DebugString("ZTildeLimit = " + std::to_string(deltaZMax) + "\n", true, true, 1); //debug


    while(stepsize >= precision)
    {
        std::vector<std::string> DataSetKTemp;
        deltaZ += stepsize;

        while(not FoundZiPlusMinus && deltaZ < deltaZMax)
        {
            std::map<std::string, double> chi2k;

            DebugString("debug: before CalculateChi2k, deltaZ = " + std::to_string(deltaZ)); //debug
            if (sign == "plus")
                chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, +deltaZ, i);
            else if (sign == "minus")
                chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, -deltaZ, i);

            if (sign == "plus") //debug
            {
                DebugString("", true, false, 1); //debug
                DebugString("deltaZ_" + std::to_string(i) + "^+ = " + std::to_string(deltaZ) + ", deltaZMax = " + std::to_string(deltaZMax) + ", stepsize = " + std::to_string(stepsize), true, true, 1); //debug
            }
            else if (sign == "minus") //debug
            {
                DebugString("", true, false, 1); //debug
                DebugString("deltaZ_" + std::to_string(i) + "^- = " + std::to_string(deltaZ) + ", deltaZMax = " + std::to_string(deltaZMax) + ", stepsize = " + std::to_string(stepsize), true, true, 1); //debug
            }

            for (std::string DataSet : DataSetK)
            {
                    DebugString("  " + DataSet + ": " + "chi2k: " + std::to_string(chi2k.at(DataSet)) + ", Xi: " + std::to_string(Xi90RescaledMap.at(DataSet)), true, true, 1); //debug

                if (chi2k.at(DataSet) >= Xi90RescaledMap.at(DataSet))
                {
                    DataSetKTemp.push_back(DataSet);
                    FoundZiPlusMinus = true;
                    
                    DebugString("Found chi2k > Xi90Rescaled : " + DataSet, true, true, 1); //debug
                }
            }
            std::cout << std::endl; //debug

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