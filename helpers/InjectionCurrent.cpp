#include "InjectionCurrent.h"
std::vector<double> InjectionCurrent::Step(const std::vector<double>& t, const double t0, const double duration, const double max) {
    /*double max = 10.0*/
    std::vector<double> out(t.size(), 0);
    for (int i = 0; i < t.size(); ++i) {
        if (t[i] >= t0 && t[i] <= t0 + duration) {
            out[i] = max;
        }
    }
    return out;
}

std::vector<double> InjectionCurrent::TwoStep(const std::vector<double>& t, double t0, double t1, double duration, double v1, double v2) {
    // double v1 = 10.0, double v2 = 50.0
    std::vector<double> out(t.size(), 0);
    for (int i = 0; i < t.size(); ++i) {
        if (t[i] >= t0 && t[i] <= t0 + duration) {
            out[i] = v1;
        }
        if (t[i] >= t1 && t[i] <= t1 + duration) {
            out[i] = v2;
        }
    }
    return out;
}

std::vector<double> InjectionCurrent::RampUp(const std::vector<double>& t, double t0, double tend, double min, double max) {
    // double min = 10.0, double max = 30.0
    std::vector<double> out(t.size(), 0);
    std::vector<double> ramp;
    double step = (max - min) / (tend - t0);
    for (double value = min; value <= max; value += step) {
        ramp.push_back(value);
    }

    int rampIndex = 0;
    for (int i = 0; i < t.size(); ++i) {
        if (t[i] >= t0 && t[i] < tend) {
            out[i] = ramp[rampIndex++];
        }
    }
    return out;
}