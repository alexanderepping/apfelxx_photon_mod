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
 * For easier writing the values can be printed by plottingPython/calculateExperimentalData.py.
 * The Energies, F2Gamma and yError values are taken from 
 * "The Photon Structure from Deep Inelastic Electron-Photon Scattering" by R. Nisius, Appendix D
 * (https://arxiv.org/abs/hep-ex/9912049v1).
 */
const std::map<std::string, std::map<std::string, std::vector<double>>> experimentalData = {
{"ALEPH",      {{"Energies", {9.9, 9.9, 9.9, 9.9, 20.7, 20.7, 20.7, 20.7, 284, 284, 284}},
                {"xData",    {0.0425, 0.14, 0.3, 0.6, 0.0645, 0.195, 0.385, 0.695, 0.19, 0.5, 0.81}},
                {"xError",   {0.0375, 0.06, 0.1, 0.2, 0.0555, 0.075, 0.115, 0.195, 0.16, 0.15, 0.16}},
                {"F2Gamma",  {0.3, 0.4, 0.41, 0.27, 0.36, 0.34, 0.56, 0.45, 0.65, 0.7, 1.28}},
                {"yError",   {0.03, 0.07, 0.1, 0.16, 0.05, 0.12, 0.11, 0.12, 0.14, 0.25, 0.37}}}},
{"AMY",        {{"Energies", {6.8, 6.8, 6.8, 73, 73, 73, 390, 390}},
                {"xData",    {0.07, 0.25, 0.5, 0.25, 0.5, 0.75, 0.31, 0.65}},
                {"xError",   {0.055, 0.125, 0.125, 0.125, 0.125, 0.125, 0.19, 0.15}},
                {"F2Gamma",  {0.337, 0.302, 0.322, 0.65, 0.6, 0.65, 0.94, 0.82}},
                {"yError",   {0.053, 0.049, 0.097, 0.1, 0.16, 0.14, 0.25, 0.19}}}}};
