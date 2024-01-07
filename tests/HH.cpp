#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>


typedef boost::array< double, 4 > state_type;
typedef boost::array< double, 8 > state_type_params;

std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables
std::vector<double> y = { -68, 0.1, 0, 0 }; // initial conditions

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
}

int main(int argc, char** argv)
{

    double time_start = 0.0; 
    double time_end = 1; 
    int n_points = 50; 
    double h_step = time_end / n_points;
    boost::numeric::odeint::integrate(HH, y, 0.0, 10.0, 0.1, write_model);
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
