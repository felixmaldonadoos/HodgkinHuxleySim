#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "gatingfunctions.h"

typedef boost::array< double, 4 > state_type;
typedef boost::array< double, 8 > state_type_params;

std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables
std::vector<double> y = { -68, 0.1, 0, 0 }; // initial conditions

std::vector<std::vector<double>> y_store;

void HH(const std::vector<double>& y, std::vector<double>& dxdt, double t)
{   
    dxdt[0] = (1 / p[0]) * (p[1] - p[2] * pow(y[1], 3) * y[2] * (y[0] - p[3]) - p[4] * pow(y[3], 4) * (y[0] - p[5]) - p[6] * (y[0] - p[7]));
    dxdt[1] = alpha_m(y[0]) * (1 - y[1]) - beta_m(y[0]) * y[1];
    dxdt[2] = alpha_h(y[0]) * (1 - y[2]) - beta_h(y[0]) * y[2];
    dxdt[3] = alpha_n(y[0]) * (1 - y[3]) - beta_n(y[0]) * y[3];
}

void write_model(const std::vector<double>& y, const double t)
{
    std::cout << t << ',' << y[0] << ',' << y[1] << ',' << y[2] << ',' << y[3] << std::endl;
    y_store.push_back(y);
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

int main(int argc, char** argv)
{
    /* run model only */
    double time_start = 0.0;
    double time_end = 1;
    int n_points = 50;
    double h_step = time_end / n_points;
    boost::numeric::odeint::integrate(HH, y, 0.0, 10.0, 0.1, write_model);
    std::cout << "Saving to file 'output.csv'\n";
    std::cout << "Using gating functions from `gatingfunctions.cpp`!\n";
    //write_model_csv();
}

//
//
//void HHModel(const HH_params& x, HH_InitConds& dxdt, double t)
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
//
//    //dxdt[0] = sigma * (x[1] - x[0]); // voltage 
//    //dxdt[1] = R * x[0] - x[1] - x[0] * x[2]; //  current 
//    //dxdt[2] = -b * x[2] + x[0] * x[1]; // gating value 
//    //dxdt[3] = -b * x[2] + x[0] * x[1]; // Injection current
//    }
//
//    void write_HHModel(const HH_params& x, const double t)
//    {
//        std::cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << std::endl;
//    }
//
