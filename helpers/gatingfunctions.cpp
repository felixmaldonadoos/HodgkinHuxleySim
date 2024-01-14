#include "gatingfunctions.h"

struct Expression {
    Expression() {
        table.add_variable("x", x);
        table.add_constants();
        expr.register_symbol_table(table);
    }

    bool set(const std::string& _str) {
        str = _str;
        return valid = parser.compile(_str, expr);
    }

    double eval(double _x) {
        x = _x;
        return expr.value();
    }

    ImVec4 color;
    bool valid;
    std::string str;
    exprtk::symbol_table<double> table;
    exprtk::expression<double> expr;
    exprtk::parser<double> parser;
    double x;
};

double GatingFunctions::alpha_n(double V)
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
double GatingFunctions::beta_n(double V)
{
    return 0.125 * exp(-(V + 60.0) / 80.0);
}

/*! \def alpha_m
    * \brief The alpha gating function for the m gate
    *
    */
double GatingFunctions::alpha_m(double V)
{
    return (0.1 * (V + 35.0)) / (1 - exp(-(V + 35.0) / 10.0));
}

/*! \def beta_m
    * \brief The beta gating function for the m gate
    *
    */
double GatingFunctions::beta_m(double V)
{
    return 4.0 * exp(-(V + 60.0) / 18.0);
}

/*! \def alpha_h
    * \brief The alpha gating function for the h gate
    *
    */
double GatingFunctions::alpha_h(double V)
{
    return 0.07 * exp(-(V + 60.0) / 20.0);
}

/*! \def beta_h
    * \brief The beta gating function for the h gate
    *
    */
double GatingFunctions::beta_h(double V)
{
    return 1.0 / (1.0 + exp(-(V + 30.0) / 10.0));
}


