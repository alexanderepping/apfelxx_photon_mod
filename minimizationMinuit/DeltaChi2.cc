
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
            zTilde[i] += std::sqrt(std::real(EigenSolverHessian.eigenvalues()[i])) * y[j]* std::real(EigenSolverHessian.eigenvectors()(i,j));

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
            y[i] += std::real(EigenSolverHessian.eigenvectors()(j,i)) * zTilde[j] / std::sqrt(std::real(EigenSolverHessian.eigenvalues()[i]));

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
        std::vector<double> zTilde(finalParams.size(), 0.);
        zTilde[i]  = deltaZ;

        // calculate the {y_i}+-
        std::vector<double> yParams  = ZTildeToY(zTilde, EigenSolverHessian);

        // define the vectors, in which the {a_i}^(+-)
        std::vector<double> aParams(finalParams.size(), 0);

        // a_i^(+-) = a_i^0 + y_i^(+-)
        for (int j=0; j<finalParams.size(); j++)
            aParams[j]  = finalParams[j] + yParams[j];
        
        for (int j=0; j<finalParams.size(); j++)
            std::cout << std::to_string(aParams[j]) << ", ";
        std::cout << std::endl;


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

std::vector<std::vector<double>> CalculateZikPlusMinus(StructureFunctionsFcn               const& StructureFunctions,
                                                       Eigen::EigenSolver<Eigen::MatrixXd> const& EigenSolverHessian,
                                                       int                                 const& i,
                                                       std::vector<double>                 const& finalParams,
                                                       std::map<std::string, double>       const& Xi90RescaledMap)
{
    std::map<std::string, std::map<std::string, std::vector<double>>> chi2kMap;

    // setting up the chi2kMap
    for (std::string DataSet : StructureFunctions.IncludedExpData())
    {
        chi2kMap[DataSet]["+"] = {};
        chi2kMap[DataSet]["-"] = {};
    }
    // here we save the z, which lead to the chi2k, saved in this chi2kMap
    chi2kMap["z"]["+"] = {};
    chi2kMap["z"]["-"] = {};

    std::cout << "deltaZ: " << std::endl; //debug

    // for (a lot of different deltaZ)
    // here, maybe variable deltaZ?
    // for that, I need to record them and change the FindZikPlusMinusFunction
    //for (double deltaZ=0.0001; deltaZ<0.03; deltaZ+=0.0005)
    //for (double deltaZ=0.000001; deltaZ<2.119173; deltaZ+=0.05)
    for (double deltaZ=2.11917; deltaZ<3.; deltaZ+=0.000101)
    //for (double deltaZ=0.000001; deltaZ<0.050002; deltaZ+=0.01)
    {

        std::cout << "deltaZ: " << std::to_string(deltaZ) << std::endl; //debug
        std::map<std::string, double> chi2kPlus  = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, +deltaZ, i);
        std::map<std::string, double> chi2kMinus = CalculateChi2k(StructureFunctions, EigenSolverHessian, finalParams, -deltaZ, i);


        for (std::string DataSet : StructureFunctions.IncludedExpData())
        {
            std::cout << " " << std::setw(12) << DataSet << ": "; //debug
            std::cout << "+: " << std::to_string(chi2kPlus.at(DataSet)); //debug
            std::cout << ", -: " << std::to_string(chi2kMinus.at(DataSet)); //debug
            std::cout << ", Xi: " << std::to_string(Xi90RescaledMap.at(DataSet)) << std::endl; //debug
            // std::cout << "\tchi2kPlus:  " << std::to_string(chi2kPlus.at(DataSet)) << "(" << DataSet << ")" << std::endl; //debug
            // std::cout << "\tchi2kMinus: " << std::to_string(chi2kMinus.at(DataSet)) << "(" << DataSet << ")" << std::endl; //debug

            chi2kMap[DataSet]["+"].push_back(chi2kPlus.at(DataSet));
            chi2kMap[DataSet]["-"].push_back(chi2kMinus.at(DataSet));
        }

        chi2kMap["z"]["+"].push_back(+deltaZ);
        chi2kMap["z"]["-"].push_back(-deltaZ);
    }
    std::cout << std::endl; //debug

    // vectors, in which the z_i^{(k)+} and z_i^{(k)-} for all experiments are saved
    std::vector<double> zikPlus;
    std::vector<double> zikMinus;

    for (std::string DataSet : StructureFunctions.IncludedExpData())
    {
        std::cout << "DataSet: " << DataSet << std::endl; //debug

        // calculate the z_i^{(k)+} and z_i^{(k)-} for specific experiment (see nCTEQ15, (A4))
        zikPlus.push_back(FindZikPlusMinus(chi2kMap.at(DataSet).at("+"), chi2kMap.at("z").at("+"), Xi90RescaledMap.at(DataSet)));
        zikMinus.push_back(FindZikPlusMinus(chi2kMap.at(DataSet).at("-"), chi2kMap.at("z").at("-"), Xi90RescaledMap.at(DataSet)));
    }
        std::cout << std::endl; //debug
        //std::cout << "size of zikPlus: " << std::to_string(zikPlus.size()) << std::endl; //debug

    return {zikPlus, zikMinus};
}