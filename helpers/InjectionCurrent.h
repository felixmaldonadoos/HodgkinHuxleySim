#pragma once
#include <vector>

class InjectionCurrent {
public:
    static std::vector<double> Step(const std::vector<double>& t, const double t0, const double duration, const double max = 10.0);
    static std::vector<double> TwoStep(const std::vector<double>& t, double t0, double t1, double duration, double v1 = 10.0, double v2 = 50.0);
    static std::vector<double> RampUp(const std::vector<double>& t, double t0, double tend, double min = 10.0, double max = 30.0);
    
};