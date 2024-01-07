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

bool PLOT_MODEL = true; 

/* stores all data points for plotting */
std::vector<std::vector<double>> y_store;

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

        /* check if plot only still model */
        if (PLOT_MODEL) 
        {
          /*  static float xs1[1001], ys1[1001];
            for (int i = 0; i < 1001; ++i) {
                xs1[i] = i * 0.001f;
                ys1[i] = 0.5f + 0.5f * sinf(50 * (xs1[i] + (float)ImGui::GetTime() / 10));
            }
            static double xs2[20], ys2[20];
            for (int i = 0; i < 20; ++i) {
                xs2[i] = i * 1 / 19.0f;
                ys2[i] = xs2[i] * xs2[i];
            }*/
            int n = y_store.size();
            static double y1[150], y2[150], y3[150], y4[150];
            static double x[150];
            if (n == 150) { // * bug: try to prevent seg fault* 

                static double ys1[150];

                for (int i = 0; i < n; ++i) {
                    x[i] = i;
                    ys1[i] = 0.5f + 0.5f * sinf(50 * (x[i] + ImGui::GetTime() / 10));
                    y1[i] = y_store[i][1];
                    y2[i] = y_store[i][2];
                    y3[i] = y_store[i][3];
                    y4[i] = y_store[i][4];
                    //x[i] = 0;
                    //y1[i] = 0.5f + 0.5f * sinf(50 * (xs1[i] + (float)ImGui::GetTime() / 10));
                }
                if (ImPlot::BeginPlot("Line Plots")) {
                    ImPlot::SetupAxes("x", "y");
                    ImPlot::PlotLine("f(x)", x, y1, n);
                    ImPlot::PlotLine("g(x)", x, ys1, n);
                 /*
                    ImPlot::PlotLine("f(x)", x, y2, n);
                    ImPlot::PlotLine("f(x)", x, y3, n);
                    ImPlot::PlotLine("f(x)", x, y4, n);*/
                    ImPlot::EndPlot();
                }
            }
            else {
                std::cout << "size does not match.";
            }
        }
        else 
        {
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
        }
    ImGui::End();
    }
};

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
    y_store.push_back(y);
}



void run_model() {
    /* run model only */
    double time_start = 0.0;
    double time_end = 30;
    double h_step = 0.1;
    boost::numeric::odeint::integrate(HH, y, time_start, time_end, h_step, write_model);
    std::cout << "finished solving model. Printing all datapoints\n";
    for (const auto y : y_store) {
        std::cout << y[0] << "," << y[1] << "," << y[2] << std::endl;
    }
}
int main(int argc, char const* argv[])
{
    /* 1 to plot model only (still) and 0 if you want to plot all channels (not ready) */
    PLOT_MODEL = true; 

    /* run model only */
    run_model();

    /* run GUI only */
    ImGraph app("Hodkin-Huxley Simulation by Felix A. Maldonado", 640, 480, argc, argv);
    app.Run();
    return 0;
}
