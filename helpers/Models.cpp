#include "Models.h"
#include "gatingfunctions.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>


double Models::C_m = 1.0;
double Models::g_Na = 120.0;
double Models::g_K = 36.0;
double Models::g_L = 0.3;
double Models::E_Na = 50.0;
double Models::E_K = -74;
double Models::E_L = -54.387;

void Models::Models::SetC_m(double value) {

}

void Models::Setg_Na(double value) {

}

void Models::Setg_K(double value) {

}

void Models::Setg_L(double value) {

}

void Models::SetE_Na(double value) {

}

void Models::SetE_K(double value) {

}

void Models::SetE_L(double value) {

}


double Models::I_Na(double V, double m, double h) {

    //Sodium membrane current(in uA / cm ^ 2)
    return g_Na * pow(m, 3) * h * (V - E_Na);
}

double Models::I_K(double V, double n) {
    //Potassium membrane current(in uA / cm ^ 2)
    return g_K * pow(n, 4) * (V - E_K);
}

double Models::I_L(double V) {

    //Leak membrane current(in uA / cm ^ 2)

    return g_L * (V - E_L);
}

double Models::I_inj(double t_query) {
    // t_query: current time 
    // Find the index of the closest time point
    int index = 0;
    std::vector<double> time_space_local = time_space;
    for (int i = 0; i < time_space.size(); ++i) {
        time_space_local[i] = time_space_local[i] - t_query;
    }

    double smallest = 10000000;
    int i_smallest = -1;
    for (int i = 0; i < time_space.size(); ++i) {
        if (time_space_local[i] < smallest) {
            i_smallest = i;
        }

    }
    // Calculate the index
    /*if (i_smallest < 0) {
        std::cout << "[ERROR] double I_inj(double t_query)\n";
        exit(1);
    }*/
    return stepCurrent[i_smallest];
}

//void Models::HH_Complete(const std::vector<double>& y, std::vector<double>& dxdt, double t)
//{
//    /* data stored in y:
//    * y[0] = Voltage
//    * y[1] = m
//    * y[2] = h
//    * y[3] = n
//    */
//    //I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) }; 
//
//    dxdt[0] = ((I_inj(t)) - I_store[0] - I_store[1] - I_store[3]) / C_m;
//    dxdt[1] = GatingFunctions::alpha_m(y[0]) * (1.0 - y[1]) - GatingFunctions::beta_m(y[0]) * y[1];
//    dxdt[2] = GatingFunctions::alpha_h(y[0]) * (1.0 - y[2]) - GatingFunctions::beta_h(y[0]) * y[2];
//    dxdt[3] = GatingFunctions::alpha_n(y[0]) * (1.0 - y[3]) - GatingFunctions::beta_n(y[0]) * y[3];
//}


void Models::HH_Basic(const std::vector<double>& y, std::vector<double>& dxdt, double t) {
    dxdt[0] = (Models::I_inj(t)) - Models::I_Na(y[0], y[1], y[2]) - Models::I_K(y[0], y[3]) - Models::I_L(y[0]) / C_m;
    dxdt[1] = GatingFunctions::alpha_m(y[0]) * (1 - y[1]) - GatingFunctions::beta_m(y[0]) * y[1];
    dxdt[2] = GatingFunctions::alpha_h(y[0]) * (1 - y[2]) - GatingFunctions::beta_h(y[0]) * y[2];
    dxdt[3] = GatingFunctions::alpha_n(y[0]) * (1 - y[3]) - GatingFunctions::beta_n(y[0]) * y[3];
}