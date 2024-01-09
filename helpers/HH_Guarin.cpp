#pragma once
#include "gatingfunctions.h"
#include <cmath>
#include <algorithm>
#include <numeric>


std::vector<double> _t; 
double _C_m = 1.0;
double _g_Na = 120.0;
double _g_K = 36.0;
double _g_L = 0.3;
double _E_Na = 50.0;
double _E_K = 77;
double _E_L = -54.387;

double I_Na_v2(double V, double m, double h) {

    //Sodium membrane current(in uA / cm ^ 2)
    return _g_Na * pow(m,3) * h * (V - _E_Na);
}

double I_K_v2(double V, double n) {
    //Potassium membrane current(in uA / cm ^ 2)
    return _g_K * pow(n,4) * (V - _E_K);
}

double I_L_v2(double V) {

    //Leak membrane current(in uA / cm ^ 2)

    return _g_L * (V - _E_L);
}

double I_inj(double t_query) {
    // Find the index of the closest time point
    int index = 0;
    double min_diff = std::numeric_limits<double>::max();
    for (size_t i = 0; i < _t.size(); ++i) {
        double diff = std::abs(_t[i] - t_query);
        if (diff < min_diff) {
            min_diff = diff;
            index = i;
        }
    }
}

double alpha_n(double V)
{
    /*! alpha_n
    * The alpha gating function for the n gate
    *
    */
    return (0.01 * (V + 50.0)) / (1.0 - exp(-(V + 50.0) / 10.0));
}

/*! \def beta_n
    * \brief The beta gating function for the n gate
    *
    */
double beta_n(double V)
{
    return 0.125 * exp(-(V + 60.0) / 80.0);
}

/*! \def alpha_m
    * \brief The alpha gating function for the m gate
    *
    */
double alpha_m(double V)
{
    return (0.1 * (V + 35.0)) / (1 - exp(-(V + 35.0) / 10.0));
}

/*! \def beta_m
    * \brief The beta gating function for the m gate
    *
    */
double beta_m(double V)
{
    return 4.0 * exp(-(V + 60.0) / 18.0);
}

/*! \def alpha_h
    * \brief The alpha gating function for the h gate
    *
    */
double alpha_h(double V)
{
    return 0.07 * exp(-(V + 60.0) / 20.0);
}

/*! \def beta_h
    * \brief The beta gating function for the h gate
    *
    */
double beta_h(double V)
{
    return 1.0 / (1.0 + exp(-(V + 30.0) / 10.0));
}


void HH_Guarin_Model(const std::vector<double>&y, std::vector<double>&dxdt, double t)
{
    /* data stored in y: 
    * y[0] = Voltage 
    * y[1] = m 
    * y[2] = h
    * y[3] = n
    */

    dxdt[0] = (I_inj(t)) - I_Na_v2(y[0], y[1], y[2]) - I_K_v2(y[0], y[3]) - I_L_v2(y[0]) / _C_m;
    dxdt[1] = alpha_m(y[0]) * (1.0 - y[1]) - beta_m(y[0]) * y[1];
    dxdt[2] = alpha_h(y[0]) * (1.0 - y[2]) - beta_h(y[0]) * y[2];
    dxdt[3] = alpha_n(y[0]) * (1.0 - y[3]) - beta_n(y[0]) * y[3];
}

