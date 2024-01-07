// Demo:   graph.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   6/7/2021

#include "App.h"
#include <exprtk.hpp>
#include <iostream>
#include <chrono>
#include <thread>
#include <imgui_stdlib.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

template <typename T>
static inline T remap(T x, T x0, T x1, T y0, T y1)
{
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

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

struct ScrollingBuffer {
    int MaxSize;
    int Offset;
    ImVector<ImVec2> Data;
    ScrollingBuffer(int max_size = 2000) {
        MaxSize = max_size;
        Offset = 0;
        Data.reserve(MaxSize);
    }
    void AddPoint(float x, float y) {
        if (Data.size() < MaxSize)
            Data.push_back(ImVec2(x, y));
        else {
            Data[Offset] = ImVec2(x, y);
            Offset = (Offset + 1) % MaxSize;
        }
    }
    void Erase() {
        if (Data.size() > 0) {
            Data.shrink(0);
            Offset = 0;
        }
    }
};

// utility structure for realtime plot
struct RollingBuffer {
    float Span;
    ImVector<ImVec2> Data;
    RollingBuffer() {
        Span = 10.0f;
        Data.reserve(2000);
    }
    void AddPoint(float x, float y) {
        float xmod = fmodf(x, Span);
        if (!Data.empty() && xmod < Data.back().x)
            Data.shrink(0);
        Data.push_back(ImVec2(xmod, y));
    }
};

void Sparkline(const char* id, const float* values, int count, float min_v, float max_v, int offset, const ImVec4& col, const ImVec2& size) {
    ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(0, 0));
    if (ImPlot::BeginPlot(id, size, ImPlotFlags_CanvasOnly)) {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoDecorations, ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxesLimits(0, count - 1, min_v, max_v, ImGuiCond_Always);
        ImPlot::SetNextLineStyle(col);
        ImPlot::SetNextFillStyle(col, 0.25);
        ImPlot::PlotLine(id, values, count, 1, 0, ImPlotLineFlags_Shaded, offset);
        ImPlot::EndPlot();
    }
    ImPlot::PopStyleVar();
}

struct ImGraph : App {

    Expression expr;
    ImPlotRect limits;

    using App::App;

    void Start() override {
        expr.set("0.25*sin(2*pi*5*x)+0.5");
        expr.color = ImVec4(1, 0.75f, 0, 1);
    }

    void Update() override {

        ImGui::SetNextWindowSize(GetWindowSize());
        ImGui::SetNextWindowPos({ 0,0 });

        ImGui::Begin("ImGraph", nullptr, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize);
        bool valid = expr.valid;
        if (!valid)
            ImGui::PushStyleColor(ImGuiCol_FrameBg, { 1,0,0,1 });
        if (ImGui::InputText("f(x)", &expr.str, ImGuiInputTextFlags_EnterReturnsTrue))
            expr.set(expr.str);
        if (!valid)
            ImGui::PopStyleColor();
        ImGui::SameLine();
        ImGui::ColorEdit4("##Color", &expr.color.x, ImGuiColorEditFlags_NoInputs);

        ImGui::BulletText("Move your mouse to change the data!");
        ImGui::BulletText("This example assumes 60 FPS. Higher FPS requires larger buffer size.");
        static ScrollingBuffer sdata1, sdata2;
        static RollingBuffer   rdata1, rdata2;
        ImVec2 mouse = ImGui::GetMousePos();
        static float t = 0;
        t += ImGui::GetIO().DeltaTime;
        sdata1.AddPoint(t, mouse.x * 0.0005f);
        rdata1.AddPoint(t, mouse.x * 0.0005f);
        sdata2.AddPoint(t, mouse.y * 0.0005f);
        rdata2.AddPoint(t, mouse.y * 0.0005f);

        static float history = 10.0f;
        ImGui::SliderFloat("History", &history, 1, 30, "%.1f s");
        rdata1.Span = history;
        rdata2.Span = history;

        static ImPlotAxisFlags flags = ImPlotAxisFlags_NoTickLabels;
        /* Voltage */
        if (ImPlot::BeginPlot("Voltage", ImVec2(-1, 150))) {
            ImPlot::SetupAxes("time (ms) ", "V (uV)", flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, t - history, t, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1);
            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::PlotShaded("Mouse X", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), -INFINITY, 0, sdata1.Offset, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &sdata2.Data[0].x, &sdata2.Data[0].y, sdata2.Data.size(), 0, sdata2.Offset, 2 * sizeof(float));
            ImPlot::EndPlot();
        }

        /* Current  */
        if (ImPlot::BeginPlot("Current", ImVec2(-1, 150))) {
            ImPlot::SetupAxes("time (ms) ", "I (uA)", flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, history, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1);
            ImPlot::PlotLine("Mouse X", &rdata1.Data[0].x, &rdata1.Data[0].y, rdata1.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::EndPlot();
        }

        /* Gating Value*/
        if (ImPlot::BeginPlot("Gating Value", ImVec2(-1, 150))) {
            ImPlot::SetupAxes("t (ms)", "Gating Value", flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, history, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1);
            ImPlot::PlotLine("Mouse X", &rdata1.Data[0].x, &rdata1.Data[0].y, rdata1.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::EndPlot();
        }

        /* Injection current (I_inj) */
        if (ImPlot::BeginPlot("##Injection current (I_inj)", ImVec2(-1, 150))) {
            ImPlot::SetupAxes("t (ms)", "I_inj (uA/cm**2)", flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, history, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1);
            ImPlot::PlotLine("Mouse X", &rdata1.Data[0].x, &rdata1.Data[0].y, rdata1.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::EndPlot();
        }
    ImGui::End();
    }
};

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

typedef boost::array< double, 3 > state_type;

void lorenz(const state_type& x, state_type& dxdt, double t)
{
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

void write_lorenz(const state_type& x, const double t)
{
    std::cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << std::endl;
}



struct push_back_state_and_time
{
    /*! Stores each set of points during the integration step.
     *
     * This structre stores the time value and the variable set in one file, this is
     * the libboost odeint observer so the vectors x and t contain each time step value for time,
     * voltaje, and gating probabilities h, m, and n.
     *
     * */

    std::vector< std::vector<double> >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time(std::vector< std::vector<double> >& states, std::vector< double >& times)
        : m_states(states), m_times(times) { }

    void operator()(const std::vector<double>& x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};


/*!
 *  \brief Gating functions for the model
 *
 *  The Hodgking-Huxley model it's endowed with a set of auxiliary gating functions
 *  Those functions proposed on the original paper of the 1952 called
 *  _A quantitative description of membrane current and it's applications to conduction and excitation in nerve_
 *   published at J.physiol. (1952) 117, 500-544 it's free online
 *
 * In other books like _Nonlinear dynamics in physiology and medicine_ from Anne Beuter et al.
 * or _Mathematical aspects of the Hodgking-Huxley Neural theory_ from Janne Cronin
 * are other formulations for the gating functions
 */


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


class hh_model
{
private:

    std::vector<double> parameters;

public:

    hh_model(std::vector<double> params) : parameters(params) { }

    void operator()(const std::vector<double>& y, std::vector<double>& f, const double /* t */)
    {
        f[0] = (1 / parameters[0]) * (parameters[1] - parameters[2] * pow(y[1], 3) * y[2] * (y[0] - parameters[3]) - parameters[4] * pow(y[3], 4) * (y[0] - parameters[5]) - parameters[6] * (y[0] - parameters[7]));
        f[1] = alpha_m(y[0]) * (1 - y[1]) - beta_m(y[0]) * y[1];
        f[2] = alpha_h(y[0]) * (1 - y[2]) - beta_h(y[0]) * y[2];
        f[3] = alpha_n(y[0]) * (1 - y[3]) - beta_n(y[0]) * y[3];
    }
};

int main(int argc, char const* argv[])
{
    state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions
    boost::numeric::odeint::integrate(lorenz, x, 0.0, 25.0, 0.1, write_lorenz);

    ///* solve HH model */
    ////HH_InitConds y = { -68, 0.1, 0, 0}; // initial conditions
    //double time_start = 0.0;
    //double time_final = 25.0; 
    //double h_step = time_final / 100; // 100 points
    //boost::numeric::odeint::integrate(HHModel, y, time_start, time_final, h_step, write_HHModel);

    /* run GUI */
    //ImGraph app("Hodkin-Huxley Simulation by Felix A. Maldonado", 640, 480, argc, argv);
    //app.Run();
    return 0;
}
