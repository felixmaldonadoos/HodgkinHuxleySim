// Demo:   graph.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   6/7/2021

#include "App.h"
#include <exprtk.hpp>
#include <iostream>
//#include <functional>
//#include <chrono>
//#include <thread>
#include <imgui_stdlib.h>
#include <sciplot/sciplot.hpp>
#include "gatingfunctions.h"
#include "InjectionCurrent.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "Models.h"

/*================= INPUT ====================*/
bool bPlotModel = false;
/* ===== simulation setup =====*/
double time_start = 0.0; // will not put as constants becuase eventually
double time_end = 50.0; // will be user-defined in GUI
double dt = 0.1; //  will be user-defined in GUI
double h_step = 0.1;
std::vector<double> y = { -65, 0.05, 0.6, 0.32 }; // initial conditions V,m, h, n


/* used to bind function to odeint*/
std::vector<double> dxdt;
double t;
/* run model only */

/* plotting */
static float history = 10.0f;

/*=================== MAIN ====================*/

const double time_start_injection = time_start + 10.0; // 15 ms delay
const double time_duration_injection = 30.0;
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

/* injection current */
std::vector<double> stepCurrent;
/* solve model over time domain */
boost::numeric::odeint::runge_kutta4<std::vector<double>> stepper;

//template <typename T>
//static inline T remap(T x, T x0, T x1, T y0, T y1)
//{
//    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
//}


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
    // Calculate the index
    /*if (i_smallest < 0) {
        std::cout << "[ERROR] double I_inj(double t_query)\n";
        exit(1);
    }*/
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
    t_now = t_now + dt;
    t_store.push_back(t_now);
    V_store.push_back(y[0]);
    y_store.push_back(y);

    I_Na_store.push_back(I_store[0]);
    I_K_store.push_back(I_store[1]);
    I_L_store.push_back(I_store[2]);
    
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

    dxdt[0] = ((I_inj(t)) - I_store[0] - I_store[1] - I_store[3]) / C_m;
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
    y = { -65, 0.05, 0.6, 0.32 }; // initial conditions
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

void ModelInit() {
    n_samples = (int)time_end / dt; // 500 (ms) / 0.1 (dt, step) = 5000 samples
    time_space = linspace(time_start, time_end, n_samples); 

    y = { -65, 0.05, 0.6, 0.32 }; // initial conditions
    stepCurrent = InjectionCurrent::Step(time_space, time_start_injection, time_duration_injection, max_injection_amplitude);

    /* prepare initial ion concentrations with initial conditions */
    I_store = { I_Na(y[0], y[1], y[2]), I_K(y[0], y[3]), I_L(y[0]) };

    std::cout << "[ModelInit()] Initialized model.";
    
}

void ModelDoStep() {
    stepper.do_step(Model, y, t, dt); // basic model
    n_samples_now += 1; 
    t_now += dt; 
}

void run_model() {

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

//struct ScrollingBuffer {
//    int MaxSize;
//    int Offset;
//    ImVector<ImVec1> Data;
//    ScrollingBuffer(int max_size = n_samples) {
//        MaxSize = max_size;
//        Offset = 0;
//        Data.reserve(MaxSize);
//    }
//    void AddPoint(double x) {
//        if (Data.size() < MaxSize)
//            Data.push_back(x);
//        else {
//            Data[Offset] = x;
//            Offset = (Offset + 1) % MaxSize;
//        }
//    }
//    void Erase() {
//        if (Data.size() > 0) {
//            Data.shrink(0);
//            Offset = 0;
//        }
//    }
//};

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

    using App::App;

    void Start() override {
        expr.set("0.25*sin(2*pi*5*x)+0.5");
        expr.color = ImVec4(1, 0.75f, 0, 1);

        /* initialize model */
        

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

        if (n_samples_now <= n_samples) {
            ModelDoStep(); // does n_samples_now += 1; 
            ModelWrite();
        }
        else {
            ModelReset();
        }
        //if (n_samples_now <= n_samples) { ModelDoStep(); ModelWrite(); else { ModelReset(); } }


        static ScrollingBuffer t_now_data, V_data, m_data, h_data, n_data;
        static ScrollingBuffer2 V;

        V.AddPoint((float)t_now, (float)y[0]);
        /* PRIORITY: will change, just want to go to sleep so bad but want to see it run. sorry.*/
        t_now_data.AddPoint(t_now);
        V_data.AddPoint(y[0]);
        m_data.AddPoint(y[1]);
        h_data.AddPoint(y[2]);
        n_data.AddPoint(y[3]);

        /* demo stuff */
        static ScrollingBuffer2 sdata1;
        ImVec2 mouse = ImGui::GetMousePos();
        static float t = 0;
        t += ImGui::GetIO().DeltaTime;
        sdata1.AddPoint(t, (float)y[0]);

        static float history = 10.0f;
        ImGui::SliderFloat("History", &history, 1, 30, "%.1f s");
        const float span = history;
        static ImPlotAxisFlags flags;

        /* end demo */
 
        /* STEP HERE */
                                                      
        //if (ImPlot::BeginPlot("##Scrolling", ImVec2(-1, 150))) {
        //    ImPlot::SetupAxes(nullptr, nullptr, flags, flags); // to do: add time/Volt in label
        //    ImPlot::SetupAxisLimits(ImAxis_X1, t - history, t, ImGuiCond_Always);
        //    ImPlot::SetupAxisLimits(ImAxis_Y1, -100.0f, -30.0f );
        //    ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
        //    //ImPlot::PlotLine("Mouse Y", &sdata2.Data[0].x, &sdata2.Data[0].y, sdata2.Data.size(), 0, sdata2.Offset, 2 * sizeof(float));
        //    ImPlot::PlotLine("Mouse Y", &V.Data[0].x, &V.Data[0].y, V.Data.size(), 0, V.Offset, 2 * sizeof(float));
        //    //ImPlot::PlotLine("Mouse Y", &t, (double)0.5, sdata2.Data.size(), 0, sdata2.Offset, 2 * sizeof(float));
        //   
        //    ImPlot::EndPlot();
        //}

        if (ImPlot::BeginPlot("##Scrolling", ImVec2(-1, 150))) {
            ImPlot::SetupAxes(nullptr, nullptr, flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, t - history, t, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -85.0, -50.0);
            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::PlotLine("Mouse X", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), 0, sdata1.Offset, 2 * sizeof(float));
            ImPlot::EndPlot();
        }
      
    }

    void Demo_RealtimePlots() {
        ImGui::BulletText("Move your mouse to change the data!");
        ImGui::BulletText("This example assumes 60 FPS. Higher FPS requires larger buffer size.");
        static ScrollingBuffer2 sdata1, sdata2;
        static RollingBuffer   rdata1, rdata2;
//+       static RollingBufferV  t_buff, V_buff;

        ImVec2 mouse = ImGui::GetMousePos();
        static float t = 0;
        t += ImGui::GetIO().DeltaTime;
        sdata1.AddPoint(t, mouse.x * 0.0005f);
        rdata1.AddPoint(t, mouse.x * 0.0005f);
        sdata2.AddPoint(t, mouse.y * 0.0005f);
        rdata2.AddPoint(t, mouse.y * 0.0005f);

        ImGui::SliderFloat("History", &history, 1, 30, "%.1f s");
        rdata1.Span = history;
        rdata2.Span = history;

        static ImPlotAxisFlags flags = ImPlotAxisFlags_NoTickLabels;

        if (ImPlot::BeginPlot("##Scrolling", ImVec2(-1, 150))) {
            ImPlot::SetupAxes(nullptr, nullptr, flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, t - history, t, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1);
            ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
            ImPlot::PlotShaded("Mouse X", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), -INFINITY, 0, sdata1.Offset, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &sdata2.Data[0].x, &sdata2.Data[0].y, sdata2.Data.size(), 0, sdata2.Offset, 2 * sizeof(float));
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("##Rolling", ImVec2(-1, 150))) {
            ImPlot::SetupAxes(nullptr, nullptr, flags, flags);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, history, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1);
            ImPlot::PlotLine("Mouse X", &t, &rdata1.Data[0].y, rdata1.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::EndPlot();
        }
    }

    void Update() override {

        this->InitializePlot();

        /* check if plot only still model */
        if (true) 
        {

            PlotModel();
        }
        else 
        {
            this->Demo_RealtimePlots();
        }
    ImGui::End();
    }
};

int main(int argc, char const* argv[])
{
    /* run model only */
    //run_model();
    ModelInit();
    /* run GUI only */
    ImGraph app("Hodkin-Huxley Simulation by Felix A. Maldonado", 640, 480, argc, argv);
    app.Run();
    return 0;
}
