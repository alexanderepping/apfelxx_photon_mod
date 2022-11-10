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




std::map<std::string, double> SecondDerivativeMap(StructureFunctionsFcn const& function,
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

#ifdef StandardCentralDifferences
    std::function<std::map<std::string, double>(std::vector<double> const&)> FirstDerivative = [&] (std::vector<double> const& a0) -> std::map<std::string, double> 
        {
            std::map<std::string, double> map1 = function.Chi2PerExperiment(a1(a0, parameter1, 1));
            std::map<std::string, double> map2 = function.Chi2PerExperiment(a1(a0, parameter1, -1));
            std::map<std::string, double> map3;

            for (std::string DataSet : function.IncludedExpData())
                map3[DataSet] = (map1.at(DataSet) - map2.at(DataSet)) / (2 * h);
            map3["All"] = (map1.at("All") - map2.at("All")) / (2 * h);

            return map3;
        };

    std::map<std::string, double> map1 = FirstDerivative(a1(a0, parameter2, 1));
    std::map<std::string, double> map2 = FirstDerivative(a1(a0, parameter2, -1));
    std::map<std::string, double> map3;

    for (std::string DataSet : function.IncludedExpData())
        map3[DataSet] = (map1.at(DataSet) - map2.at(DataSet)) / (2 * h);
    map3["All"] = (map1.at("All") - map2.at("All")) / (2 * h);

    return map3;
#endif //StandardCentralDifferences

#ifdef SevenPointLowNoise
    std::function<std::map<std::string, double>(std::vector<double> const&)> FirstDerivative = [&] (std::vector<double> const& a0) -> std::map<std::string, double> 
        {
            std::map<std::string, double> map1 = function.Chi2PerExperiment(a1(a0, parameter1, +1));
            std::map<std::string, double> map2 = function.Chi2PerExperiment(a1(a0, parameter1, +2));
            std::map<std::string, double> map3 = function.Chi2PerExperiment(a1(a0, parameter1, +3));
            std::map<std::string, double> map4 = function.Chi2PerExperiment(a1(a0, parameter1, -1));
            std::map<std::string, double> map5 = function.Chi2PerExperiment(a1(a0, parameter1, -2));
            std::map<std::string, double> map6 = function.Chi2PerExperiment(a1(a0, parameter1, -3));
            std::map<std::string, double> map7;

            for (std::string DataSet : function.IncludedExpData())
                map7[DataSet] = ((map1.at(DataSet) - map4.at(DataSet)) 
                           + 2 * (map2.at(DataSet) - map5.at(DataSet)) 
                           + 3 * (map3.at(DataSet) - map6.at(DataSet))) / (28 * h);
            map7["All"] = ((map1.at("All") - map4.at("All")) 
                       + 2 * (map2.at("All") - map5.at("All")) 
                       + 3 * (map3.at("All") - map6.at("All"))) / (28 * h);

            return map7;
        };

    std::map<std::string, double> map1 = FirstDerivative(a1(a0, parameter2, +1));
    std::map<std::string, double> map2 = FirstDerivative(a1(a0, parameter2, +2));
    std::map<std::string, double> map3 = FirstDerivative(a1(a0, parameter2, +3));
    std::map<std::string, double> map4 = FirstDerivative(a1(a0, parameter2, -1));
    std::map<std::string, double> map5 = FirstDerivative(a1(a0, parameter2, -2));
    std::map<std::string, double> map6 = FirstDerivative(a1(a0, parameter2, -3));
    std::map<std::string, double> map7;

    for (std::string DataSet : function.IncludedExpData())
        map7[DataSet] = ((map1.at(DataSet) - map4.at(DataSet)) 
                    + 2 * (map2.at(DataSet) - map5.at(DataSet)) 
                    + 3 * (map3.at(DataSet) - map6.at(DataSet))) / (28 * h);
    map7["All"] = ((map1.at("All") - map4.at("All")) 
                + 2 * (map2.at("All") - map5.at("All")) 
                + 3 * (map3.at("All") - map6.at("All"))) / (28 * h);

    return map7;
#endif //SevenPointLowNoise
};



double SecondDerivative(StructureFunctionsFcn const& function,
                        std::vector<double>   const& a0,    
                        int                   const& parameter1, 
                        int                   const& parameter2, 
                        double                const& h) 
{
    return SecondDerivativeMap(function, a0, parameter1, parameter2, h).at("All");
};



std::map< std::string, Eigen::MatrixXd> CalculateHessianMap(StructureFunctionsFcn const& function,
                                                         std::vector<double>   const& a0,
                                                         double                const& h) 
{
    const int NumberOfFreeParams = NumberOfFreeParams;
    std::map<std::string, Eigen::MatrixXd> Hessian;

    for (std::string DataSet : function.IncludedExpData())
        Hessian[DataSet] = Eigen::MatrixXd(NumberOfFreeParams, NumberOfFreeParams);
    Hessian["All"] = Eigen::MatrixXd(NumberOfFreeParams, NumberOfFreeParams);

        //Hessian[DataSet] = Eigen::MatrixXd::Zero(NumberOfFreeParams, NumberOfFreeParams);
    //Hessian["All"] = Eigen::MatrixXd::Zero(NumberOfFreeParams, NumberOfFreeParams);

    for (int i=0; i<NumberOfFreeParams; i++)
    {
        for (int j=i; j<NumberOfFreeParams; j++)
        {
            std::map<std::string, double> SD = SecondDerivativeMap(function, a0, i, j, h);

            for (std::string DataSet : function.IncludedExpData())
            {
                Hessian.at(DataSet)(i, j) = 0.5 * SD.at(DataSet);
                Hessian.at(DataSet)(j, i) = 0.5 * SD.at(DataSet);
            }

            Hessian.at("All")(i, j) = 0.5 * SD.at("All");
            Hessian.at("All")(j, i) = 0.5 * SD.at("All");
        }
    }

    return Hessian;
};



Eigen::MatrixXd CalculateHessian(StructureFunctionsFcn const& function,
                                 std::vector<double>   const& a0,
                                 double                const& h) 
{
    return CalculateHessianMap(function, a0, h).at("All");
};
