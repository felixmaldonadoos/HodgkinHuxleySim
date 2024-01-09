#include "Models.h"
#include "gatingfunctions.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

void HH(const std::vector<double>& y, std::vector<double>& dxdt, double t)
{
    //
//
//    /*
//     * \brief Defines the O.DE system proposed by Hodgking and Huxley as a class
//     *
//     * The default constructor receives a std::vector<double> that contains the parameters of the model in the order
//     *   * Membrane capacitance
//     *   * Induced current on axon, 0 means no external current
//     *   * Na conductances
//     *   * Na Nernst Potential
//     *   * K Conductance
//     *   * K Nernst Potential
//     *   * Leakage conductance (Due to a Cl current)
//     *   * Leakage Nernst potential (Due to a Cl current)
//     * Default values: 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42
//     * Initial conditions: -68, 0.1, 0, 0
//     * So parameter[0] contains membrane capacitance and so on
//     * https://github.com/Daniel-M/Hodgking-Huxley/blob/master/include/HH_Model_class.hpp
//     */
//    dxdt[0] = (1 / x[0]) * (x[1] - x[2] * pow(y[1], 3) * y[2] * (y[0] - x[3]) - x[4] * pow(y[3], 4) * (y[0] - x[5]) - x[6] * (y[0] - x[7]));
//    dxdt[1] = alpha_m(y[0]) * (1 - y[1]) - beta_m(y[0]) * y[1];
//    dxdt[2] = alpha_h(y[0]) * (1 - y[2]) - beta_h(y[0]) * y[2];
//    dxdt[3] = alpha_n(y[0]) * (1 - y[3]) - beta_n(y[0]) * y[3];

    std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables

    dxdt[0] = (1 / p[0]) * (p[1] - p[2] * pow(y[1], 3) * y[2] * (y[0] - p[3]) - p[4] * pow(y[3], 4) * (y[0] - p[5]) - p[6] * (y[0] - p[7]));
    dxdt[1] = alpha_m(y[0]) * (1 - y[1]) - beta_m(y[0]) * y[1];
    dxdt[2] = alpha_h(y[0]) * (1 - y[2]) - beta_h(y[0]) * y[2];
    dxdt[3] = alpha_n(y[0]) * (1 - y[3]) - beta_n(y[0]) * y[3];
}



