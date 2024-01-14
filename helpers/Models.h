#pragma once
#include "InjectionCurrent.h"
#include <iostream>
#include <vector>

class Models { public: 

    /* injection current */
    std::vector<double> stepCurrent;
    /* model constants */
    std::vector<double> t;
    std::vector<double> time_space;

    /* i can do 2 things:
    * 1. case-switch statement and pass in a string f_value_to_change and f_double new_value. gui has a vector of values that cab
    * 2. when gui value is set, depending on the button
    */

    static double C_m;
    static double g_Na;
    static double g_K;
    static double g_L;
    static double E_Na;
    static double E_K;
    static double E_L;
    /* dont forget ab INjectioncurrent::stepper */


    void SetC_m(double value);
    void Setg_Na(double value);
    void Setg_K(double value);
    void Setg_L(double value);
    void SetE_Na(double value);
    void SetE_K(double value);
    void SetE_L(double value);

    std::vector<double> y_store; 

    //Sodium membrane current(in uA / cm ^ 2)
    static double I_Na(double V, double m, double h);

    static double I_K(double V, double n);

    static double I_L(double V);
    //external current(in uA / cm ^ 2)
    static double I_inj(double t_query);
    
   /* void HH_Complete(const std::vector<double>& y, std::vector<double>& dxdt, double t);*/

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

    static void HH_Basic(const std::vector<double>& y, std::vector<double>& dxdt, double t);
};
