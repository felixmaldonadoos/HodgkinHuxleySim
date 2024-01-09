#pragma once

#include <iostream>
#include <fstream>
#include <vector>
//#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <sciplot/sciplot.hpp>
//#include "Models.h"
#include <cmath>
//#include "matplotlibcpp.h"
#include "InjectionCurrent.h"
#include "gatingfunctions.h"
//#include "HH_Guarin.h"

//typedef boost::array< double, 4 > state_type;
//typedef boost::array< double, 8 > state_type_params;

//std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables
std::vector<double> y = { -65, 0.05, 0.6, 0.32 }; // initial conditions
std::vector<double> y_model = { -68, 0.1, 0, 0 }; // will be used by the model



/* safety */
int MAX_SAMPLES = 5000;

/* used to bind function to odeint*/
std::vector<double> dxdt;
double t; 
int n_samples = 0;
/* run model only */
double time_start = 0.0;
double time_end = 1000.0;
double t_now = time_start;
double dt = 0.1;

std::vector<double> time_space;
std::vector<std::vector<double>> y_store;
std::vector<double> t_store    = { t_now };
std::vector<double> V_store    = { y[0] };
std::vector<double> I_store    = { y[1],  y[2] , y[3] }; // Na, K, L
std::vector<double> I_Na_store = { y[1]}; 
std::vector<double> I_K_store  = { y[2] }; 
std::vector<double> I_L_store  = { y[3] }; 

/* injection current */
std::vector<double> stepCurrent;

/* to do: make this a seperate file */

std::vector<double> _t;
double _C_m = 1.0;
double _g_Na = 120.0;
double _g_K = 36.0;
double _g_L = 0.3;
double _E_Na = 50.0;
double _E_K = -74;
double _E_L = -54.387;

double I_Na(double V, double m, double h) {

    //Sodium membrane current(in uA / cm ^ 2)
    return _g_Na * pow(m, 3) * h * (V - _E_Na);
}

double I_K(double V, double n) {
    //Potassium membrane current(in uA / cm ^ 2)
    return _g_K * pow(n, 4) * (V - _E_K);
}

double I_L(double V) {

    //Leak membrane current(in uA / cm ^ 2)

    return _g_L * (V - _E_L);
}

double I_inj(double t_query) {
    // t_query: current time 
    // Find the index of the closest time point
    int index = 0;

    for (int i = 0; i < time_space.size(); ++i) {
        time_space[i] = time_space[i] - t_query;
    }

    double smallest = 10000000;
    int i_smallest = -1;
    for (int i = 0; i < time_space.size(); ++i) {
        if (time_space[i] < smallest) {
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

void Model(const std::vector<double>& y, std::vector<double>& dxdt, double t)
{
    /* data stored in y:
    * y[0] = Voltage
    * y[1] = m
    * y[2] = h
    * y[3] = n
    */
    //I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) }; 

    dxdt[0] = ((I_inj(t)) - I_store[0] - I_store[1] - I_store[3]) / _C_m;
    dxdt[1] = GatingFunctions::alpha_m(y[0]) * (1.0 - y[1]) - GatingFunctions::beta_m(y[0]) * y[1];
    dxdt[2] = GatingFunctions::alpha_h(y[0]) * (1.0 - y[2]) - GatingFunctions::beta_h(y[0]) * y[2];
    dxdt[3] = GatingFunctions::alpha_n(y[0]) * (1.0 - y[3]) - GatingFunctions::beta_n(y[0]) * y[3];
}

void HH(const std::vector<double>& y, std::vector<double>& dxdt, double t){
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

    //std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables

    //dxdt[0] = (1 / p[0]) * (p[1] - p[2] * pow(y[1], 3) * y[2] * (y[0] - p[3]) - p[4] * pow(y[3], 4) * (y[0] - p[5]) - p[6] * (y[0] - p[7]));
    dxdt[0] = (I_inj(t)) - I_Na(y[0], y[1], y[2]) - I_K(y[0], y[3]) - I_L(y[0]) / _C_m;
    dxdt[1] = GatingFunctions::alpha_m(y[0]) * (1 - y[1]) - GatingFunctions::beta_m(y[0]) * y[1];
    dxdt[2] = GatingFunctions::alpha_h(y[0]) * (1 - y[2]) - GatingFunctions::beta_h(y[0]) * y[2];
    dxdt[3] = GatingFunctions::alpha_n(y[0]) * (1 - y[3]) - GatingFunctions::beta_n(y[0]) * y[3];
}

void write_model(const std::vector<double>& y, const double t)
{
    I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) };

    //++n_samples;
    //if (n_samples == MAX_SAMPLES) { 
    //    exit(1); 
    //}
    t_now = t_now + dt; 
    t_store.push_back(t_now);
    V_store.push_back(y[0]);
    y_store.push_back(y);

    I_Na_store.push_back(I_store[0]);
    I_K_store.push_back(I_store[1]);
    I_L_store.push_back(I_store[2]);

    std::cout << t_now  << ',' << I_store[0] << ',' << I_store[1] << ',' << I_store[2] << std::endl;
    /* to do: update y_model to keep y */
}

void write_model_csv() {

    // Open a CSV file for writing
    std::ofstream file("C:/Users/famal/Documents/Projects/HodkinHuxleySim/extras/output.csv");

    // Check if file is open
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    // Iterate over the data and write to the file
    for (const auto& row : y_store) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) file << ","; // Add comma after each element except the last
        }
        file << "\n"; // End of row
    }

    // Close the file
    file.close();
}

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> linspaced;

    if (num == 0) {
        return linspaced;
    }
    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i) {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // Ensure that end is exactly at 'end'

    return linspaced;
}

void MyPlot() {
    https://sciplot.github.io/tutorials/
    // Create a vector with values from 0 to pi divived into 200 uniform intervals for the x-axis

    // Create a Plot object
    sciplot::Plot2D plot;
    std::cout << "size timespace: " << time_space.size() << "\nsize V_store: " << V_store.size() << "\nsize t_store: " << t_store.size() << std::endl;
    std::cout << "size I_Na: " << I_Na_store.size() << "\nsize I_K: " << I_K_store.size() << "\nsize I_L_store: " << I_L_store.size() << std::endl;

    // Set the x and y labels
    plot.xlabel("x");
    plot.ylabel("y");

    // Set the x and y ranges
    plot.xrange(time_start, time_end);
    std::min_element(V_store.begin(), V_store.end());
    plot.yrange(*std::min_element(V_store.begin(), V_store.end()), 
        *std::max_element(V_store.begin(), V_store.end()));

    // Set the legend to be on the bottom along the horizontal
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    // Plot sin(i*x) from i = 1 to i = 6
    plot.drawCurve(t_store, V_store).label("V (mV)");
    /*plot.drawCurve(x, std::sin(2.0 * x)).label("sin(2x)");
    plot.drawCurve(x, std::sin(3.0 * x)).label("sin(3x)");
    plot.drawCurve(x, std::sin(4.0 * x)).label("sin(4x)");
    plot.drawCurve(x, std::sin(5.0 * x)).label("sin(5x)");
    plot.drawCurve(x, std::sin(6.0 * x)).label("sin(6x)");*/

    // Create figure to hold plot
    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };

    // Show the plot in a pop-up window
    canvas.show();

    /* currents */
    sciplot::Plot2D plot_currents;
    std::cout << "size I_Na: " << I_Na_store.size() << "\nsize I_K: " << I_K_store.size() << "\nsize I_L_store: " << I_L_store.size() << std::endl;

    // Set the x and y labels
    plot_currents.xlabel("x");
    plot_currents.ylabel("y");

    // Set the x and y ranges
    plot_currents.xrange(time_start, time_end);
    std::min_element(I_Na_store.begin(), I_Na_store.end());
    plot_currents.yrange(-5, +5);

    plot_currents.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    plot_currents.drawCurve(t_store, I_Na_store).label("I_Na");
    plot_currents.drawCurve(t_store, I_K_store).label("I_K");
    plot_currents.drawCurve(t_store, I_L_store).label("I_L");


    // Create figure to hold plot
    sciplot::Figure fig_currents = { {plot_currents} };
    // Create canvas to hold figure
    sciplot::Canvas canvas_currents = { {fig_currents} };

    // Show the plot in a pop-up window
    canvas_currents.show();


}

int main(int argc, char** argv)
{
    std::cout << "Starting program.\n";
    time_start = 0.0;
    time_end = 500.0;
    time_space = linspace(time_start, time_end, (int)time_end/dt); // 500 (ms) / 0.1 (dt, step) = 5000 samples
    //dt = time_end / MAX_SAMPLES;
    const double time_start_injection = time_start + 50.0;
    const double time_duration_injection = 50.0;
    const double max_injection_amplitude = 10.0;
    
    /* 
    initial conditions 
    * y[0] = Voltage (mV)
    * y[1] = m
    * y[2] = h
    * y[3] = n
    */
    y = { -65, 0.05, 0.6, 0.32 };
    stepCurrent = InjectionCurrent::Step(time_space, time_start_injection, time_duration_injection, max_injection_amplitude);
    //boost::numeric::odeint::integrate(Model, y, time_start, time_end, dt, write_model);
    
    /* prepare initial ion concentrations with initial conditions */
    I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) };

    /* solve model over time domain */
    boost::numeric::odeint::runge_kutta4<std::vector<double>> stepper;
    for (double t = time_start; t < time_end; t += dt) {
        stepper.do_step(HH, y, t, dt);
        write_model(y, t);
    }

    /* basic plot */
    MyPlot();    
    
    return 0;
}
