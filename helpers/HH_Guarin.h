#pragma once
#include <vector>

std::vector<double> _t;
double _C_m = 1.0;
double _g_Na = 120.0;
double _g_K = 36.0;
double _g_L = 0.3;
double _E_Na = 50.0;
double _E_K = 77;
double _E_L = -54.387;

double I_Na_v2(double V, double m, double h);
double I_K_v2(double V, double n);
double I_L_v2(double V);
double I_inj(double t_query);

/* gating functions */
double alpha_n(double V);
double beta_n(double V);
double alpha_m(double V);
double beta_m(double V);
double alpha_h(double V);
double beta_h(double V);

void HH_Guarin_Model(const std::vector<double>& y, std::vector<double>& dxdt, double t);

