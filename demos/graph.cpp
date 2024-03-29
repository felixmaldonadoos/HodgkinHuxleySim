// Demo:   graph.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   6/7/2021

#include "App.h"
#include <exprtk.hpp>
#include <iostream>
#include <imgui_stdlib.h>

template <typename T>
static inline T remap(T x, T x0, T x1, T y0, T y1) 
{ return y0 + (x - x0) * (y1 - y0) / (x1 - x0); }


struct Expression {
    Expression() {
        table.add_variable("x",x);
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

struct ImGraph : App {

    Expression expr;
    ImPlotRect limits;

    using App::App;

    void Start() override {
        expr.set("0.25*sin(2*pi*5*x)+0.5");
        expr.color = ImVec4(1,0.75f,0,1);
    }

    void Update() override {
   
        ImGui::SetNextWindowSize(GetWindowSize());
        ImGui::SetNextWindowPos({0,0});

        ImGui::Begin("ImGraph",nullptr,ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize);
        bool valid = expr.valid;
        if (!valid)
            ImGui::PushStyleColor(ImGuiCol_FrameBg, {1,0,0,1});
        if (ImGui::InputText("f(x)",&expr.str,ImGuiInputTextFlags_EnterReturnsTrue)) 
            expr.set(expr.str);
        if (!valid)
            ImGui::PopStyleColor();
        ImGui::SameLine();
        ImGui::ColorEdit4("##Color",&expr.color.x,ImGuiColorEditFlags_NoInputs);
        /*
        if (ImPlot::BeginPlot("##Plot",0,0,ImVec2(-1,-1),0,ImPlotAxisFlags_NoInitialFit,ImPlotAxisFlags_NoInitialFit)) 
        {
            limits = ImPlot::GetPlotLimits();
            if (valid) {
                ImPlot::SetNextLineStyle(expr.color);
                ImPlot::PlotLineG("##item",
                    [](int idx, void* data) {
                    auto& self = *(ImGraph*)data;
                    double x = remap((double)idx, 0.0, 9999.0, self.limits.X.Min, self.limits.X.Max);
                    double y = self.expr.eval(x);
                    return ImPlotPoint(x,y);
                    },
                    this,
                    10000);
            }
            ImPlot::EndPlot();
        }
        */
        
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
            ImPlot::PlotLine("Mouse X", &rdata1.Data[0].x, &rdata1.Data[0].y, rdata1.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 0, 2 * sizeof(float));
            ImPlot::EndPlot();
        }
        ImGui::End();
    }
};

int main(int argc, char const *argv[])
{
    ImGraph app("ImGraph",640,480,argc,argv);
    app.Run();
    return 0;
}
