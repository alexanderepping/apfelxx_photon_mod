/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include "StructureFunctionsFcn.h"
#include "configMinuit.h"
#include "HelperFunctions.h"

#include <Eigen/Eigenvalues>
#include <vector>
#include <functional>




std::map<std::string, double> SecondDerivativeMap(StructureFunctionsFcn const& function,
                                               std::vector<double>   const& a0,    
                                               int                   const& parameter1, 
                                               int                   const& parameter2, 
                                               double                const& h) 
{
    DebugString("      in SecondDerivativeMap - debug 01"); //debug
    std::function<std::vector<double>(std::vector<double> const&, int const&, int const&)> a1 = [&] (std::vector<double> const& a0, int const& parameter, int const& k) -> std::vector<double>
    {
        std::vector<double> a1 = a0;
        a1[parameter] += k*h;
        return a1;
    };

#ifdef StandardCentralDifferences
/**
 * FirstDerivative = [X²({a_1,...,a_i+h,...a_n}) - X²({a_1,...,a_i-h,...a_n})]/[2h], 
 * where a is a1(a0, j, +1) or a1(a0, j, -1) and i=parameter1 and j=parameter2
*/
    std::function<std::map<std::string, double>(std::vector<double> const&)> FirstDerivative = [&] (std::vector<double> const& a) -> std::map<std::string, double> 
        {
            DebugString("        in FirstDerivative - debug 01"); //debug
            std::map<std::string, double> map1 = function.Chi2PerExperiment(a1(a, parameter1, 1));
            DebugString("        in FirstDerivative - debug 02"); //debug
            std::map<std::string, double> map2 = function.Chi2PerExperiment(a1(a, parameter1, -1));
            DebugString("        in FirstDerivative - debug 03"); //debug
            std::map<std::string, double> map3;

            DebugString("        in FirstDerivative - debug 04"); //debug
            for (std::string DataSet : function.IncludedExpData())
                map3[DataSet] = (map1.at(DataSet) - map2.at(DataSet)) / (2 * h);
            map3["All"] = (map1.at("All") - map2.at("All")) / (2 * h);
            DebugString("        in FirstDerivative - debug 05"); //debug

            return map3;
        };

    DebugString("      in SecondDerivativeMap - debug 02"); //debug
    // map1 = [X²({a0_1,...,a0_i+h,...,a0_j+h,...a0_n}) - X²({a0_1,...,a0_i-h,...,a0_j+h,...a0_n})]/[2h], where i=parameter1 and j=parameter2
    std::map<std::string, double> map1 = FirstDerivative(a1(a0, parameter2, 1));
    DebugString("      in SecondDerivativeMap - debug 03"); //debug
    // map2 = [X²({a0_1,...,a0_i+h,...,a0_j-h,...a0_n}) - X²({a0_1,...,a0_i-h,...,a0_j-h,...a0_n})]/[2h], where i=parameter1 and j=parameter2
    std::map<std::string, double> map2 = FirstDerivative(a1(a0, parameter2, -1));
    DebugString("      in SecondDerivativeMap - debug 04"); //debug
    std::map<std::string, double> map3;

    DebugString("      in SecondDerivativeMap - debug 05"); //debug
    for (std::string DataSet : function.IncludedExpData())
        map3[DataSet] = (map1.at(DataSet) - map2.at(DataSet)) / (2 * h);
    map3["All"] = (map1.at("All") - map2.at("All")) / (2 * h);
    DebugString("      in SecondDerivativeMap - debug 06"); //debug

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
    std::map<std::string, Eigen::MatrixXd> Hessian;

    for (std::string DataSet : function.IncludedExpData())
        Hessian[DataSet] = Eigen::MatrixXd(NumberOfFreeParams, NumberOfFreeParams);
    Hessian["All"] = Eigen::MatrixXd(NumberOfFreeParams, NumberOfFreeParams);

    for (int i=0; i<NumberOfFreeParams; i++)
    {
        DebugString("in first for-loop - i = "+std::to_string(i)); //debug
        for (int j=i; j<NumberOfFreeParams; j++)
        {
            DebugString("    in second for-loop - j = "+std::to_string(j)+", before SD"); //debug
            std::map<std::string, double> SD = SecondDerivativeMap(function, a0, i, j, h);
            DebugString("    in second for-loop - j = "+std::to_string(j)+", after SD"); //debug

            for (std::string DataSet : function.IncludedExpData())
            {
                Hessian.at(DataSet)(i, j) = 0.5 * SD.at(DataSet);
                Hessian.at(DataSet)(j, i) = 0.5 * SD.at(DataSet);
            }

            DebugString("    in second for-loop - j = "+std::to_string(j)+", after saving SD"); //debug

            Hessian.at("All")(i, j) = 0.5 * SD.at("All");
            Hessian.at("All")(j, i) = 0.5 * SD.at("All");
        }
    }

    return Hessian;
};



Eigen::MatrixXd CalculateHessianMatrix(StructureFunctionsFcn const& function,
                                       std::vector<double>   const& a0,
                                       double                const& h) 
{
    return CalculateHessianMap(function, a0, h).at("All");
};
