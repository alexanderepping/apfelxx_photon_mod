
/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "DeltaChi2.h"
#include "configMinuit.h"

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
    const int n = y.size();
    std::vector<double> zTilde(n, 0);

    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            // z_i^Tilde = sqrt(lambda_i) * Sum_j [ y_j * V_j^(i) ]
            // here was an error with the indizes of the eigenvectors
            zTilde[i] += std::sqrt(std::real(EigenSolverHessian.eigenvalues()[i])) * y[j] * std::real(EigenSolverHessian.eigenvectors()(j,i));

    return zTilde;
}

std::vector<double> ZTildeToY(std::vector<double>                 const& zTilde,
                              Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian)
{
    const int n = zTilde.size();
    std::vector<double> y(n, 0);

    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            // y_i = Sum_j [ 1/sqrt(lambda_j) * V_i^(j) * z_i^Tilde ]
            // here was an error with the indizes of the eigenvectors
            y[i] += std::real(EigenSolverHessian.eigenvectors()(i,j)) * zTilde[j] / std::sqrt(std::real(EigenSolverHessian.eigenvalues()[j]));

    return y;
}

std::vector<double> FindZTildeLimit(int                                 const& j,
                                    std::vector<double>                 const& finalParams,
                                    Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian)
{
    int size = finalParams.size();
    std::vector<double> zTildeMax(size, 0);

    for (int i=0; i<size; i++)
    {
        // calculated from nCTEQ15 (2.24), but only zTilde_j with the specific j (given to this function) is not zero
        zTildeMax[i] = (initialParamsLBounds.at(usedInitialPDFs)[i] - finalParams[i]) * std::sqrt(std::real(EigenSolverHessian.eigenvalues()[j])) / std::real(EigenSolverHessian.eigenvectors()(i,j));

        //debug beginning
        std::cout << "§  a_" << std::to_string(j) << " = " << std::to_string(initialParamsLBounds.at(usedInitialPDFs)[i]) << " for zTilde_" << std::to_string(i) << " = " << std::to_string(zTildeMax[i]) << std::endl; //debug
        /*
        std::vector<double> zTildep(size, 0.);
        zTildep[i]  = zTildeMax[j]+1;
        std::vector<double> yp = ZTildeToY(zTildep, EigenSolverHessian);
        std::cout << "  a_i+ = ";
        for (int j=0; j<finalParams.size(); j++)
            std::cout << std::to_string(yp[j] + finalParams[j]) << ", ";
        std::cout << std::endl;
        std::vector<double> zTildem(size, 0.);
        zTildem[i]  = zTildeMax[j]-1;
        std::vector<double> ym = ZTildeToY(zTildem, EigenSolverHessian);
        std::cout << "  a_i- = ";
        for (int j=0; j<finalParams.size(); j++)
            std::cout << std::to_string(ym[j] + finalParams[j]) << ", ";
        std::cout << std::endl;
        */
        //debug end
    }

    return zTildeMax;
}

int FindDeltaChi2Limit(std::vector<double>                 const& finalParams,
                       Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian)
{
    int size = finalParams.size();

    // loop through all the different functions fk-
    for (int k=0; k<size; k++)
    {
        std::cout << "k = " << std::to_string(k) << std::endl;
        
        // loop throug all the parameters
        for (int i=0; i<size; i++)
        {
            // calculated from nCTEQ15 (2.32)
            double DeltaChi2Max = std::pow(-(initialParamsLBounds.at(usedInitialPDFs)[i] - finalParams[i]) / std::real(EigenSolverHessian.eigenvectors()(i,k)), 2.) * std::real(EigenSolverHessian.eigenvalues()[k]);
            std::cout << "§  a_" << std::to_string(i) << " = " << std::to_string(initialParamsLBounds.at(usedInitialPDFs)[i]) << " for DeltaChi2" << " < " << std::to_string(DeltaChi2Max) << std::endl; //debug

            std::cout << "§     -a_i- = ";
            for (int j=0; j<size; j++)
                std::cout << std::to_string(finalParams[j] - std::sqrt((DeltaChi2Max-1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j))) << ", ";
            std::cout << std::endl;
            std::cout << "§     -a_i  = ";
            for (int j=0; j<size; j++)
                std::cout << std::to_string(finalParams[j] - std::sqrt((DeltaChi2Max+0.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(j,k))) << ", ";
            std::cout << std::endl;
            std::cout << "§     -a_i+ = ";
            for (int j=0; j<size; j++)
                std::cout << std::to_string(finalParams[j] - std::sqrt((DeltaChi2Max+1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j))) << ", ";
            std::cout << std::endl;

            std::cout << "§     +a_i- = ";
            for (int j=0; j<size; j++)
                std::cout << std::to_string(finalParams[j] + std::sqrt((DeltaChi2Max-1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j))) << ", ";
            std::cout << std::endl;
            std::cout << "§     +a_i  = ";
            for (int j=0; j<size; j++)
                std::cout << std::to_string(finalParams[j] + std::sqrt((DeltaChi2Max+0.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(j,k))) << ", ";
            std::cout << std::endl;
            std::cout << "§     +a_i+ = ";
            for (int j=0; j<size; j++)
                std::cout << std::to_string(finalParams[j] + std::sqrt((DeltaChi2Max+1.) / std::real(EigenSolverHessian.eigenvalues()[k])) * std::real(EigenSolverHessian.eigenvectors()(k,j))) << ", ";
            std::cout << std::endl;
        }
    }
    return 0;
}


std::map<std::string, double> CalculateChi2k(StructureFunctionsFcn           const& StructureFunctions,
                                         Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                                         std::vector<double>                 const& finalParams,
                                         double                              const& deltaZ,
                                         int                                 const& i) 
{
        // make the zTilde vectors, but there is only a variation in one ziTilde direction
        // the variation deltaZ is in positive and negative direction
        std::vector<double> zTilde(finalParams.size(), 0.);
        zTilde[i]  = deltaZ;

        // calculate the {y_i}+-
        std::vector<double> yParams  = ZTildeToY(zTilde, EigenSolverHessian);

        // define the vectors, in which the {a_i}^(+-)
        std::vector<double> aParams(finalParams.size(), 0);

        // a_i^(+-) = a_i^0 + y_i^(+-)
        for (int j=0; j<finalParams.size(); j++)
            aParams[j]  = finalParams[j] + yParams[j];
        
        for (int j=0; j<finalParams.size(); j++) //debug
            std::cout << std::to_string(aParams[j]) << ", "; //debug
        std::cout << std::endl; //debug


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
    
    std::cout << "chi2k: " << std::to_string(chi2kData[j]) << ", Xi90Rescaled: " << std::to_string(Xi90Rescaled); //debug
    std::cout << ", z: " << std::to_string(zikData[j]) << ", Index: " << std::to_string(j) << "/" << std::to_string(chi2kData.size() -1) << std::endl; //debug

    return zikData[j];
}

std::vector<double> CalculateZikPlusMinus(StructureFunctionsFcn               const& StructureFunctions,
                                                       Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                                                       int                                 const& i,
                                                       std::vector<double>                 const& finalParams,
                                                       std::map<std::string, double>       const& Xi90RescaledMap,
                                                       std::string                         const& sign,
                                                       double                              const& deltaZmax,
                                                       double                              const& deltaZstepsize)
{
    std::map<std::string, std::vector<double>> chi2kMap;

    // setting up the chi2kMap
    for (std::string DataSet : StructureFunctions.IncludedExpData())
        chi2kMap[DataSet] = {};

    // here we save the z, which lead to the chi2k, saved in this chi2kMap
    chi2kMap["z"] = {};

    std::cout << "deltaZ: " << std::endl; //debug

    // for (a lot of different deltaZ)
    // here, maybe variable deltaZ?
    // for that, I need to record them and change the FindZikPlusMinusFunction
    // the maximum for "minus" is 2.119172
    for (double deltaZ=0.000001; deltaZ<=deltaZmax; deltaZ+=deltaZstepsize)
    {

        std::cout << "deltaZ: " << std::to_string(deltaZ) << std::endl; //debug
        std::map<std::string, double> chi2k;
        if (sign == "plus")
            chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, +deltaZ, i);
        else if (sign == "minus")
            chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, -deltaZ, i);


        for (std::string DataSet : StructureFunctions.IncludedExpData())
        {
            std::cout << " " << std::setw(12) << DataSet << ": "; //debug
            std::cout << "chi2k: " << std::to_string(chi2k.at(DataSet)); //debug
            std::cout << ", Xi: " << std::to_string(Xi90RescaledMap.at(DataSet)) << std::endl; //debug

            chi2kMap[DataSet].push_back(chi2k.at(DataSet));
        }

        if (sign == "plus")
            chi2kMap["z"].push_back(+deltaZ);
        else if (sign == "minus")
            chi2kMap["z"].push_back(-deltaZ);
    }
    std::cout << std::endl; //debug

    // vectors, in which the z_i^{(k)+} and z_i^{(k)-} for all experiments are saved
    std::vector<double> zik;

    for (std::string DataSet : StructureFunctions.IncludedExpData())
    {
        std::cout << "DataSet: " << DataSet << std::endl; //debug

        // calculate the z_i^{(k)+} and z_i^{(k)-} for specific experiment (see nCTEQ15, (A4))
        zik.push_back(FindZikPlusMinus(chi2kMap.at(DataSet), chi2kMap.at("z"), Xi90RescaledMap.at(DataSet)));
    }
        std::cout << std::endl; //debug

    return zik;
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
    std::cout << "§  Finding the ZTildeLimit for i=" << std::to_string(i) << ", " << sign << ": " << std::endl; //debug
    std::vector<double> ZTildeLimit = FindZTildeLimit(i, finalParams, EigenSolverHessian);
    for (double zTilde : ZTildeLimit)
        if ((sign == "plus" && zTilde > 0.) || (sign == "minus" && zTilde < 0.))
            if (deltaZMax == 0. || deltaZMax > std::abs(zTilde))
                deltaZMax = std::abs(zTilde);
    std::cout << "§  ZTildeLimit = " << std::to_string(deltaZMax) << std::endl << std::endl; //debug


    while(stepsize >= precision)
    {
        std::vector<std::string> DataSetKTemp;
        deltaZ += stepsize;

        while(not FoundZiPlusMinus && deltaZ < deltaZMax)
        {
            std::map<std::string, double> chi2k;

            std::cout << "debug: before CalculateChi2k, deltaZ = " << std::to_string(deltaZ) << std::endl; //debug
            if (sign == "plus")
                chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, +deltaZ, i);
            else if (sign == "minus")
                chi2k = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, -deltaZ, i);

            if (sign == "plus") //debug
                std::cout << "\n§  deltaZ_" << std::to_string(i) << "^+ = " << std::to_string(deltaZ) << ", deltaZMax = " << std::to_string(deltaZMax) << ", stepsize = " << std::to_string(stepsize) << std::endl; //debug
            else if (sign == "minus") //debug
                std::cout << "\n§  deltaZ_" << std::to_string(i) << "^- = " << std::to_string(deltaZ) << ", deltaZMax = " << std::to_string(deltaZMax) << ", stepsize = " << std::to_string(stepsize) << std::endl; //debug

            for (std::string DataSet : DataSetK)
            {
                    std::cout << "§   " << std::setw(12) << DataSet << ": "; //debug
                    std::cout << "chi2k: " << std::to_string(chi2k.at(DataSet)); //debug
                    std::cout << ", Xi: " << std::to_string(Xi90RescaledMap.at(DataSet)) << std::endl; //debug

                if (chi2k.at(DataSet) >= Xi90RescaledMap.at(DataSet))
                {
                    DataSetKTemp.push_back(DataSet);
                    FoundZiPlusMinus = true;
                    
                    std::cout << "§  Found chi2k > Xi90Rescaled : " << DataSet << std::endl; //debug
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