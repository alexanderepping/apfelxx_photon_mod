/**
 * apfelxx_photon_mod 2021
 * 
 * Author: Alexander Epping: a_eppi01@uni-muenster.de
 * GitHub: https://github.com/alexanderepping/apfelxx_photon_mod
 */

#include <string>
#include <map>
#include <vector>

/**
 * @brief Map of experimental values.
 * The data is taken from the papers referenced in 
 * "NLO photon parton parametrization using ee and ep data" by  W. Slominski, H. Abramowicz and A. Levy
 * (https://arxiv.org/abs/hep-ph/0504003v2).
 */
const std::map<std::string, std::map<std::string, std::vector<double>>> experimentalData = {
{"ALEPH1",     {{"xData",      {0.0425, 0.14, 0.30, 0.60, 0.0645, 0.195, 0.385, 0.695, 0.19, 0.50, 0.81}}, 
                {"Q2Data",     {9.9, 9.9, 9.9, 9.9, 20.7, 20.7, 20.7, 20.7, 284.0, 284.0, 284.0}}, 
                {"F2Gamma",    {0.30, 0.40, 0.41, 0.27, 0.36, 0.34, 0.56, 0.45, 0.65, 0.70, 1.28}}, 
                {"F2GammaErr", {0.03, 0.07, 0.10, 0.16, 0.05, 0.12, 0.11, 0.12, 0.14, 0.25, 0.37}}}}, 

{"ALEPH2",     {{"xData",      {0.0065, 0.0224, 0.0563, 0.1137, 0.1958, 0.3026, 0.4349, 0.6037, 0.0211, 0.0656, 0.1381, 0.2359, 0.3556, 0.4959, 0.6535, 0.8478}}, 
                {"Q2Data",     {17.3, 17.3, 17.3, 17.3, 17.3, 17.3, 17.3, 17.3, 67.2, 67.2, 67.2, 67.2, 67.2, 67.2, 67.2, 67.2}}, 
                {"F2Gamma",    {0.430, 0.270, 0.350, 0.350, 0.390, 0.460, 0.400, 0.180, 0.570, 0.430, 0.470, 0.500, 0.600, 0.660, 0.650, 0.660}}, 
                {"F2GammaErr", {0.116, 0.051, 0.044, 0.032, 0.038, 0.045, 0.150, 0.236, 0.147, 0.047, 0.050, 0.044, 0.052, 0.086, 0.137, 0.150}}}}, 

{"AMY",        {{"xData",      {0.0700, 0.2500, 0.5000, 0.2500, 0.5000, 0.7500, 0.3100, 0.6500}}, 
                {"Q2Data",     {6.8, 6.8, 6.8, 73.0, 73.0, 73.0, 390.0, 390.0}}, 
                {"F2Gamma",    {0.337, 0.302, 0.322, 0.650, 0.600, 0.650, 0.940, 0.820}}, 
                {"F2GammaErr", {0.053, 0.049, 0.097, 0.100, 0.160, 0.110, 0.250, 0.190}}}}, 

{"DELPHI",     {{"xData",      {0.0405, 0.1465, 0.3205, 0.6375}}, 
                {"Q2Data",     {12.0, 12.0, 12.0, 12.0}}, 
                {"F2Gamma",    {0.210, 0.410, 0.450, 0.450}}, 
                {"F2GammaErr", {0.067, 0.064, 0.071, 0.149}}}}, 

{"JADE",       {{"xData",      {0.05, 0.15, 0.30, 0.50, 0.75, 0.20, 0.45, 0.75}}, 
                {"Q2Data",     {24.0, 24.0, 24.0, 24.0, 24.0, 100.0, 100.0, 100.0}}, 
                {"F2Gamma",    {0.51, 0.29, 0.34, 0.59, 0.23, 0.52, 0.75, 0.90}}, 
                {"F2GammaErr", {0.15, 0.12, 0.10, 0.12, 0.12, 0.23, 0.22, 0.27}}}}, 

{"L3",         {{"xData",      {0.0035, 0.0075, 0.015, 0.025, 0.04, 0.075, 0.0075, 0.015, 0.03, 0.05, 0.08, 0.15, 0.055, 0.15, 0.25, 0.055, 0.15, 0.25, 0.40, 0.055, 0.15, 0.25, 0.40}}, 
                {"Q2Data",     {1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,  10.8, 10.8, 10.8, 15.3, 15.3, 15.3, 15.3, 23.1, 23.1, 23.1, 23.1}}, 
                {"F2Gamma",    {0.184, 0.179, 0.176, 0.191, 0.193, 0.185, 0.307, 0.282, 0.263, 0.278, 0.270, 0.252, 0.300, 0.350, 0.300, 0.370, 0.420, 0.42, 0.35, 0.40, 0.440, 0.470, 0.44}}, 
                {"F2GammaErr", {0.050, 0.023, 0.017, 0.009, 0.012, 0.027, 0.096, 0.047, 0.023, 0.015, 0.023, 0.047, 0.04, 0.04, 0.11, 0.04, 0.05, 0.09, 0.16, 0.05, 0.06, 0.06, 0.13}}}}, 

{"OPAL1",      {{"xData",      {0.0460, 0.1870, 0.4660, 0.0715, 0.2305, 0.4230, 0.6790, 0.2000, 0.4500, 0.7000, 0.0600, 0.1750, 0.4250, 0.0600, 0.1750, 0.4250, 0.0750, 0.2250, 0.4750, 0.7000, 0.0750, 0.2250, 0.4750, 0.7000, 0.0044, 0.0132, 0.0300, 0.0700, 0.0132, 0.0300, 0.0700, 0.1500}}, 
                {"Q2Data",     {7.5, 7.5, 7.5, 14.7, 14.7, 14.7, 14.7, 135.0, 135.0, 135.0, 9.0, 9.0, 9.0, 14.5, 14.5, 14.5, 30.0, 30.0, 30.0, 30.0, 59.0, 59.0, 59.0, 59.0, 1.86, 1.86, 1.86, 1.86, 3.76, 3.76, 3.76, 3.76}}, 
                {"F2Gamma",    {0.280, 0.320, 0.380, 0.380, 0.410, 0.410, 0.540, 0.650, 0.730, 0.720, 0.330, 0.290, 0.390, 0.370, 0.420, 0.390, 0.320, 0.520, 0.410, 0.460, 0.370, 0.440, 0.480, 0.510, 0.270, 0.220, 0.200, 0.230, 0.350, 0.290, 0.320, 0.320}}, 
                {"F2GammaErr", {0.036, 0.132, 0.214, 0.130, 0.063, 0.114, 0.314, 0.342, 0.113, 0.816, 0.067, 0.064, 0.310, 0.163, 0.149, 0.125, 0.117, 0.139, 0.219, 0.418, 0.286, 0.106, 0.184, 0.500, 0.076, 0.054, 0.092, 0.054, 0.085, 0.067, 0.073, 0.085}}}}, 

/* "old"/incomplete/wrong OPAL2 data
{"OPAL2",      {{"xData",      {0.0460, 0.1870, 0.4660, 0.0715, 0.2305, 0.4230, 0.6790, 0.0017, 0.0086, 0.0433, 0.2202, 0.0041, 0.0185, 0.0828, 0.3709, 0.0304, 0.1364, 0.5209, 0.0029, 0.0162, 0.0885, 0.4841, 0.0528, 0.1843, 0.5957, 0.0044, 0.0222, 0.1098, 0.5438, 0.2750, 0.5500, 0.8000}}, 
                {"Q2Data",     {5.90, 5.90, 5.9, 14.7, 14.70, 14.70, 14.70, 1.90, 1.90, 1.90, 1.90, 3.70, 3.70, 3.70, 3.70, 8.90, 8.90, 8.90, 10.70, 10.70, 10.70, 10.70, 17.50, 17.50, 17.50, 17.8, 17.80, 17.80, 17.80, 780.00, 780.00, 780.00}}, 
                {"F2Gamma",    {0.224, 0.352, 0.348, 0.325, 0.465, 0.446, 0.409, 0.269, 0.177, 0.179, 0.227, 0.269, 0.232, 0.259, 0.296, 0.221, 0.308, 0.379, 0.362,   0.263, 0.275, 0.351, 0.273, 0.375, 0.501, 0.428, 0.295, 0.336, 0.430, 0.930, 0.870, 0.970}}, 
                {"F2GammaErr", {0.018, 0.030, 0.090, 0.029, 0.038, 0.051, 0.102, 0.043, 0.019, 0.010, 0.016, 0.057, 0.026, 0.016, 0.032, 0.034, 0.018, 0.028, 0.073,  0.035, 0.032, 0.028, 0.048, 0.030, 0.038, 0.094, 0.038, 0.044, 0.035, 0.172, 0.180, 0.286}}}}, 
*/

{"OPAL2",      {{"xData",      {0.0460, 0.1870, 0.4660, 0.0715, 0.2305, 0.4230, 0.6790, 0.0017, 0.0086, 0.0433, 0.2202, 0.0041, 0.0185, 0.0828, 0.3709, 0.0304, 0.1364, 0.5209, 0.0029, 0.0162, 0.0885, 0.4841, 0.0528, 0.1843, 0.5957, 0.0044, 0.0222, 0.1098, 0.5438, 0.2750, 0.5500, 0.8000}}, 
                {"Q2Data",     {5.90,   5.90,   5.9,    14.7,   14.70,  14.70,  14.70,  1.90,   1.90,   1.90,   1.90,   3.70,   3.70,   3.70,   3.70,   8.90,   8.90,   8.90,   10.70,  10.70,  10.70,  10.70,  17.50,  17.50,  17.50,  17.8,   17.80,  17.80,  17.80,  780.00, 780.00, 780.00}},   
                {"F2Gamma",    {0.224,  0.352,  0.348,  0.325,  0.465,  0.446,  0.409,  0.269,  0.177,  0.179,  0.227,  0.269,  0.232,  0.259,  0.296,  0.221,  0.308,  0.379,  0.362,  0.263,  0.275,  0.351,  0.273,  0.375,  0.501,  0.428,  0.295,  0.336,  0.430,  0.930,  0.870,  0.970}}, 
                {"F2GammaErr", {0.026,  0.037,  0.132,  0.056,  0.045,  0.058,  0.124,  0.043,  0.019,  0.010,  0.016,  0.057,  0.026,  0.016,  0.032,  0.034,  0.018,  0.028,  0.073,  0.035,  0.032,  0.028,  0.048,  0.030,  0.038,  0.094,  0.038,  0.044,  0.035,  0.172,  0.180,  0.286}}}}, 

{"OPAL2_less", {{"xData",      {0.0460, 0.1870, 0.4660, 0.0715, 0.2305, 0.4230, 0.6790, 0.2750, 0.5500, 0.8000}}, 
                {"Q2Data",     {5.90,   5.90,   5.9,    14.7,   14.70,  14.70,  14.70,  780.00, 780.00, 780.00}},   
                {"F2Gamma",    {0.224,  0.352,  0.348,  0.325,  0.465,  0.446,  0.409,  0.930,  0.870,  0.970}}, 
                {"F2GammaErr", {0.026,  0.037,  0.132,  0.056,  0.045,  0.058,  0.124,  0.172,  0.180,  0.286}}}}, 

{"OPAL2_dop",  {{"xData",      {0.175,  0.425,  0.175,  0.425,  0.175,  0.425,  0.725,  0.175,  0.425,  0.725}}, 
                {"Q2Data",     {12.1,   12.1,   19.9,   19.9,   39.7,   39.7,   39.7,   76.4,   76.4,   76.4}},   
                {"F2Gamma",    {0.38,   0.43,   0.39,   0.49,   0.47,   0.63,   0.65,   0.55,   0.68,   0.73}}, 
                {"F2GammaErr", {0.032,  0.032,  0.032,  0.022,  0.022,  0.032,  0.067,  0.032,  0.022,  0.045}}}}, 

{"PLUTO",      {{"xData",      {0.063, 0.240, 0.535, 0.100, 0.305, 0.620, 0.145, 0.385, 0.720, 0.0535, 0.123, 0.2465, 0.4045, 0.570, 0.745, 0.175, 0.375, 0.625, 0.825}}, 
                {"Q2Data",     {2.4, 2.4, 2.4, 4.3, 4.3, 4.3, 9.2, 9.2, 9.2, 5.3, 5.3, 5.3, 5.3, 5.3, 5.3, 45.0, 45.0, 45.0, 45.0}}, 
                {"F2Gamma",    {0.204, 0.272, 0.222, 0.256, 0.295, 0.336, 0.354, 0.402, 0.492, 0.245, 0.307, 0.277, 0.329, 0.439, 0.361, 0.480, 0.550, 0.890, 0.870}}, 
                {"F2GammaErr", {0.053, 0.049, 0.064, 0.066, 0.048, 0.067, 0.093, 0.067, 0.101, 0.063, 0.078, 0.049, 0.061, 0.084, 0.093, 0.177, 0.132, 0.183, 0.274}}}}, 

{"TASSO",      {{"xData",      {0.11, 0.30, 0.50, 0.70, 0.90}}, 
                {"Q2Data",     {23.0, 23.0, 23.0, 23.0, 23.00}}, 
                {"F2Gamma",    {0.366, 0.670, 0.722, 0.693, 0.407}}, 
                {"F2GammaErr", {0.112, 0.153, 0.172, 0.176, 0.235}}}}, 

{"TOPAZ",      {{"xData",      {0.0430, 0.1380, 0.0850, 0.2400, 0.5550, 0.1900, 0.4550, 0.7850}}, 
                {"Q2Data",     {5.1, 5.1, 16.0, 16.0, 16.0, 80.0, 80.0, 80.0}}, 
                {"F2Gamma",    {0.33, 0.29, 0.60, 0.56, 0.46, 0.68, 0.83, 0.53}}, 
                {"F2GammaErr", {0.054, 0.042, 0.100, 0.098, 0.162, 0.265, 0.226, 0.216}}}}, 

{"TPC",        {{"xData",      {0.04, 0.118, 0.2295, 0.4515, 0.11, 0.279, 0.5495}}, 
                {"Q2Data",     {2.8, 2.8, 2.8, 2.8, 5.1, 5.1, 5.1}}, 
                {"F2Gamma",    {0.134, 0.234, 0.198, 0.16, 0.224, 0.373, 0.3}}, 
                {"F2GammaErr", {0.023, 0.040, 0.050, 0.040, 0.042, 0.077, 0.061}}}}};



/**
 * @brief Map of experimental values.
 * For easier writing the values can be printed by plottingPython/calculateExperimentalData.py.
 * The Q2Data, F2Gamma and F2GammaErr values are taken from 
 * "The Photon Structure from Deep Inelastic Electron-Photon Scattering" by R. Nisius, Appendix D
 * (https://arxiv.org/abs/hep-ex/9912049v1).
 */
const std::map<std::string, std::map<std::string, std::vector<double>>> experimentalDataOld = {
{"ALEPH",      {{"Q2Data",     {9.9, 9.9, 9.9, 9.9, 20.7, 20.7, 20.7, 20.7, 284, 284, 284}},
                {"xData",      {0.0425, 0.14, 0.3, 0.6, 0.0645, 0.195, 0.385, 0.695, 0.19, 0.5, 0.81}},
                {"xError",     {0.0375, 0.06, 0.1, 0.2, 0.0555, 0.075, 0.115, 0.195, 0.16, 0.15, 0.16}},
                {"F2Gamma",    {0.3, 0.4, 0.41, 0.27, 0.36, 0.34, 0.56, 0.45, 0.65, 0.7, 1.28}},
                {"F2GammaErr", {0.03, 0.07, 0.1, 0.16, 0.05, 0.12, 0.11, 0.12, 0.14, 0.25, 0.37}}}},
{"AMY",        {{"Q2Data",     {6.8, 6.8, 6.8, 73, 73, 73, 390, 390}},
                {"xData",      {0.07, 0.25, 0.5, 0.25, 0.5, 0.75, 0.31, 0.65}},
                {"xError",     {0.055, 0.125, 0.125, 0.125, 0.125, 0.125, 0.19, 0.15}},
                {"F2Gamma",    {0.337, 0.302, 0.322, 0.65, 0.6, 0.65, 0.94, 0.82}},
                {"F2GammaErr", {0.053, 0.049, 0.097, 0.1, 0.16, 0.14, 0.25, 0.19}}}}};
