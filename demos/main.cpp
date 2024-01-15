// Demo:   graph.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   6/7/2021

#include "App.h"
#include <exprtk.hpp>
#include <iostream>
#include <imgui_stdlib.h>
#include <sciplot/sciplot.hpp>
#include "gatingfunctions.h"
#include "InjectionCurrent.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/array.hpp>
#include <thread>
#include <boost/numeric/odeint.hpp>
#include "Models.h"

/*================= INPUT ====================*/
bool bPlotModel = false;
/* ===== simulation setup =====*/
double time_start = 0.0; // will not put as constants becuase eventually
double time_end = 100.0; // will be user-defined in GUI
double dt = 0.1; //  will be user-defined in GUI
double h_step = 0.1;
//std::vector<double> y_init = { -61, 1.0, 0.6, 0.32 }; // initial conditions V,m, h, n (original)
std::vector<double> y_init = { -61, 0.0, 0.6, 0.32 }; // initial conditions V,m, h, n
std::vector<double> y = y_init; // initial conditions V,m, h, n
std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables
bool bApplyStim = false;
bool bApplyStim_10 = false;


/* used to bind function to odeint*/
std::vector<double> dxdt;
double t;
/* run model only */

/* plotting */
static float history = 10.0f;

/*=================== MAIN ====================*/

const double time_start_injection = time_start + 10.0; // 15 ms delay
const double time_duration_injection = 5.0;
const double max_injection_amplitude = 10.0;

/* model constants */
std::vector<double> _t;
double C_m = 1.0;
double g_Na = 120.0;
double g_K = 36.0;
double g_L = 0.3;
double E_Na = 50.0;
double E_K = -74;
double E_L = -54.387;

/*=============================================*/

/* logging purposes*/
double t_now  = time_start;
int n_samples;
int n_samples_now = 0; 
std::vector<double> time_space;
std::vector<std::vector<double>> y_store;
std::vector<double> t_store = { t_now };
std::vector<double> V_store = { y[0] };
std::vector<double> I_store = { y[1],  y[2] , y[3] }; // Na, K, L
std::vector<double> I_Na_store = { y[1] };
std::vector<double> I_K_store = { y[2] };
std::vector<double> I_L_store = { y[3] };

/* control flow of model */
bool bModelDone = true; 
bool bIsThreadRunning = false;
/* injection current */
std::vector<double> stepCurrent;

/* solve model over time domain */
boost::numeric::odeint::runge_kutta4<std::vector<double>> stepper;

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

// utility structure for realtime plot
struct ScrollingBuffer2 {
    int MaxSize;
    int Offset;
    ImVector<ImVec2> Data;

    ScrollingBuffer2(int max_size = n_samples) {
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
struct RollingBuffer_v1 {
    float Span;
    ImVector<ImVec2> Data;
    RollingBuffer_v1() {
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

struct ScrollingBuffer {
    int MaxSize;
    int Offset;
    std::vector<double> Data;
    ScrollingBuffer(int max_size = n_samples) {
        MaxSize = max_size;
        Offset = 0;
        Data.reserve(MaxSize);
    }
    void AddPoint(double x) {
        if (Data.size() < MaxSize)
            Data.push_back(x);
        else {
            Data[Offset] = x;
            Offset = (Offset + 1) % MaxSize;
        }
    }
    void Erase() {
        if (Data.size() > 0) {
            Data.erase(Data.begin());
            Offset = 0;
        }
    }
};

static ScrollingBuffer2 V_buff(n_samples), m_buff(n_samples), h_buff(n_samples), n_buff(n_samples);

double I_Na(double V, double m, double h) {

    //Sodium membrane current(in uA / cm ^ 2)
    return g_Na * pow(m, 3) * h * (V - E_Na);
}

double I_K(double V, double n) {
    //Potassium membrane current(in uA / cm ^ 2)
    return g_K * pow(n, 4) * (V - E_K);
}

double I_L(double V) {

    //Leak membrane current(in uA / cm ^ 2)

    return g_L * (V - E_L);
}

double I_inj(double t_query) {
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
    
    return stepCurrent[i_smallest];
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

void ModelWrite(/*const std::vector<double>& y, const double t*/)
{
    I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) };
    t_store.push_back(t_now);
    V_store.push_back(y[0]);
    y_store.push_back(y);

    I_Na_store.push_back(I_store[0]);
    I_K_store.push_back(I_store[1]);
    I_L_store.push_back(I_store[2]);

    /*V_buff.AddPoint(0.0, 0.0);
    m_buff.AddPoint(t_now, (float)y[1]);
    h_buff.AddPoint(t_now, (float)y[2]);
    n_buff.AddPoint(t_now, (float)y[3]);*/
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
    //if (I_inj(t) > 0.0) {
    //    std::cout << "Stimulus applies: " << t << "time / " << I_inj(t) << " uA\n";
    //}



    dxdt[0] = (1 / p[0]) * (p[1] - p[2] * std::pow(y[1], 3) * y[2] * (y[0] - p[3]) - p[4] * std::pow(y[3], 4) * (y[0] - p[5]) - p[6] * (y[0] - p[7]));
    //dxdt[0] = ((I_inj(t)) - I_store[0] - I_store[1] - I_store[3]) / C_m;
    dxdt[1] = GatingFunctions::alpha_m(y[0]) * (1.0 - y[1]) - GatingFunctions::beta_m(y[0]) * y[1];
    dxdt[2] = GatingFunctions::alpha_h(y[0]) * (1.0 - y[2]) - GatingFunctions::beta_h(y[0]) * y[2];
    dxdt[3] = GatingFunctions::alpha_n(y[0]) * (1.0 - y[3]) - GatingFunctions::beta_n(y[0]) * y[3];
}

void ModelReset() {
    /* put initial parameters 
    * 
    * eventually user input will control this
    * 
    */
    std::cout << "{ModelReset()] Startting.\n";
    y = y_init; // initial conditions
    std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables
    double t_now = time_start;
    n_samples_now = 0;
    time_space = linspace(time_start, time_end, n_samples);
    t_now = time_start; 
    std::vector<std::vector<double>> y_store;
    t_store = { t_now };
    V_store = { y[0] };
    I_store = { y[1],  y[2] , y[3] }; // Na, K, L
    I_Na_store = { y[1] };
    I_K_store = { y[2] };
    I_L_store = { y[3] };
}

void HH(const std::vector<double>& y, std::vector<double>& dxdt, double t) {
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

    std::vector<double> p = { 0.01, 0.0, 1.20, 55.16, 0.36, -72.14, 0.003, -49.42 }; // constant variables

    dxdt[0] = (1 / p[0]) * (p[1] - p[2] * std::pow(y[1], 3) * y[2] * (y[0] - p[3]) - p[4] * std::pow(y[3], 4) * (y[0] - p[5]) - p[6] * (y[0] - p[7]));
    dxdt[1] = GatingFunctions::alpha_m(y[0]) * (1 - y[1]) - GatingFunctions::beta_m(y[0]) * y[1];
    dxdt[2] = GatingFunctions::alpha_h(y[0]) * (1 - y[2]) - GatingFunctions::beta_h(y[0]) * y[2];
    dxdt[3] = GatingFunctions::alpha_n(y[0]) * (1 - y[3]) - GatingFunctions::beta_n(y[0]) * y[3];
}

void PlotModel2() {

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
    plot.grid();


    plot.drawCurve(t_store, V_store).label("V (mV)"); // voltage plot

    // vertical line 
    std::vector<int> v_line(I_L_store.size(), 0);
    const double vline_amp = -75; // to do: make this dynamic
    const int idx_vline = (int)time_start_injection / dt;
    v_line[idx_vline] = *std::max_element(V_store.begin(), V_store.end());
    v_line[idx_vline + 1] = vline_amp; // extend vline to create filled effect


    plot.drawCurveFilled(t_store, v_line);

    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };

    // Show the plot in a pop-up window
    canvas.show();

}

void write_model(const std::vector<double>& y, const double t)
{
    I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) };

    //++n_samples;
    //if (n_samples == MAX_SAMPLES) { 
    //    exit(1); 
    //}
    t_now = t_now + dt;
    t_store.push_back(t);
    //t_store.push_back(t_now); // og
    V_store.push_back(y[0]);
    //V_buff.AddPoint(t_now, (float)y[0]);
}

void SolveModel2() {
    double dt = 0.1;
    n_samples = (int)(time_end / dt);
    size_t steps = boost::numeric::odeint::integrate(HH, y, time_start, time_end, dt, write_model);

    /*for (int i = 0; i < V_store.size(); ++i) {
        V_buff.AddPoint((float)t_store[i], (float)V_store[i]);
    }*/
}

void HandleModelThread() {
    if (!bIsThreadRunning) {
        std::cout << "Starting SolveModel thread.\n";
        bIsThreadRunning = true;

        /* create */
        std::thread t(SolveModel2);

        /* wait for thread to finish */
        t.join();
        std::cout << "Ending SolveModel thread.\n";
        bIsThreadRunning = false;
    }
    else {
        std::cout << "SolveModel thread already running.\n";
    }
}

void ModelInit() {
    n_samples = (int)time_end / dt; // 500 (ms) / 0.1 (dt, step) = 5000 samples
    time_space = linspace(time_start, time_end, n_samples); 

    y = y_init; // initial conditions
    stepCurrent = InjectionCurrent::Step(time_space, time_start_injection, time_duration_injection, max_injection_amplitude);

    /* prepare initial ion concentrations with initial conditions */
    I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) };

    std::cout << "[ModelInit()] Initialized model.";
    
}

int n_stims = 0; 
void ModelDoStep() {

    if (bApplyStim && n_stims < 10) { 
        p[2] = 10.0; 
        n_stims += 1; 
    }
    else { 
        p[2] = 0.0; 
        n_stims = 0;
        bApplyStim = false; 
    }

    stepper.do_step(Model, y, t, dt); // basic model
    n_samples_now += 1; 
    t_now += dt; 
}

void run_model() {

}

struct RollingBuffer {
    float Span;
    ImVector<ImVec2> Data;
    RollingBuffer() {
        Span = 10.0f;
        Data.reserve(n_samples);
    }
    void AddPoint(float x, float y) {
        float xmod = fmodf(x, Span);
        if (!Data.empty() && xmod < Data.back().x)
            Data.shrink(0);
        Data.push_back(ImVec2(xmod, y));
    }
};

struct ImGraph : App {

    Expression expr;
    ImPlotRect limits;

    /* thread params */
    bool bIsThreadRunning = false;

    using App::App;

    void Start() override {
        std::cout << "Starting App.\n bIsThreadRunning = false.\n";
        expr.set("0.25*sin(2*pi*5*x)+0.5");
        expr.color = ImVec4(1, 0.75f, 0, 1);
        bIsThreadRunning = false;
        ModelInit();
    }

    void InitializePlot() {
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
    }

    void PlotModel() {
        ImGui::BulletText("Move your mouse to change the data!");
        ImGui::BulletText("This example assumes 60 FPS. Higher FPS requires larger buffer size.");

        if (ImGui::Button("Apply stim: 5.0"))
        {
            bApplyStim = true;
            std::cout << "Stim applied!\n";
        }

        if (ImGui::Button("Apply stim: 10.0"))
        {
            bApplyStim = true;
            bApplyStim_10 = true;
            std::cout << "Stim applied!\n";
        }

        if (n_samples_now < n_samples) {
            ModelDoStep(); // does n_samples_now += 1; 
            ModelWrite();
        }
        else { ModelReset(); }
        
     

        static ScrollingBuffer2 V_buff, m_buff, h_buff, n_buff;
        static ScrollingBuffer2 step_current_buff; 

        static float t = 0;
        //t += ImGui::GetIO().DeltaTime;
        t += dt;

        /* add points */
        V_buff.AddPoint(t, (float)y[0]);
        m_buff.AddPoint(t, (float)y[1]);
        h_buff.AddPoint(t, (float)y[2]);
        n_buff.AddPoint(t, (float)y[3]);
        step_current_buff.AddPoint(t, (float)p[2]/10);

        /* control history with a slider*/
        static float history = (int)time_end/100 + 5;
        ImGui::SliderFloat("History", &history, 1, 100, "%.1f s");
        const float span = history;

        /* set flags if any*/
        static ImPlotAxisFlags flags;

        if (V_buff.Data.Size < 1) {
            std::cout << "V_buff empty. Waiting for model to be finished.\n";
            return;
        }

        if (ImPlot::BeginPlot("##Scrolling", ImVec2(-1, 350))) {

            ImPlot::SetupAxis(ImAxis_X1, "X1", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisLimits(ImAxis_X1, t - history, t, ImGuiCond_Always);

            ImPlot::SetupAxis(ImAxis_Y1, "Y1", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -70.0, 50.0);

            ImPlot::SetupAxis(ImAxis_Y2, "Y2", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisLimits(ImAxis_Y2, 0, 10.0);

            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
            ImPlot::PlotLine("Voltage", &V_buff.Data[0].x, &V_buff.Data[0].y, V_buff.Data.size(), 0, V_buff.Offset, 2 * sizeof(float));
            

            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlot::PlotLine("Injection Current (*10)", &step_current_buff.Data[0].x, &step_current_buff.Data[0].y, step_current_buff.Data.size(), 0, step_current_buff.Offset, 2 * sizeof(float));
           
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("##Scrolling 2", ImVec2(-1, 350))) {

            ImPlot::SetupAxis(ImAxis_X2, "X2", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisLimits(ImAxis_X2, t - history, t, ImGuiCond_Always);

            ImPlot::SetupAxis(ImAxis_Y3, "Y3", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisLimits(ImAxis_Y3, 0, 1);

            /* plot m */
            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::SetAxes(ImAxis_X2, ImAxis_Y3);
            ImPlot::PlotLine("m", &m_buff.Data[0].x, &m_buff.Data[0].y, m_buff.Data.size(), 0, m_buff.Offset, 2 * sizeof(float));
            
            /* plot h */
            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::SetAxes(ImAxis_X2, ImAxis_Y3);
            ImPlot::PlotLine("h", &h_buff.Data[0].x, &h_buff.Data[0].y, h_buff.Data.size(), 0, h_buff.Offset, 2 * sizeof(float));
            
            /* plot n */
            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::SetAxes(ImAxis_X2, ImAxis_Y3);
            ImPlot::PlotLine("n", &n_buff.Data[0].x, &n_buff.Data[0].y, n_buff.Data.size(), 0, n_buff.Offset, 2 * sizeof(float));

            ImPlot::EndPlot();
        }
    }

    void Update() override {
        this->InitializePlot();
        PlotModel();
        ImGui::End();
    }
};

int main(int argc, char const* argv[])
{
    /* run model only */
    //run_model();
    ModelInit();
    /* run GUI only */

    bool bRunModel2 = false;
    if (bRunModel2) {
        HandleModelThread();
        PlotModel2();
    }
    else {
        ImGraph app("Hodkin-Huxley Simulation by Felix A. Maldonado", 900, 720, argc, argv);
        app.Run();
    }
    return 0;
}
