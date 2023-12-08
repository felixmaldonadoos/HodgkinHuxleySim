// MIT License

// Copyright (c) 2023 Evan Pezent

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// Demo:   graph.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   6/7/2021

#include "App.h"
#include <exprtk.hpp>
#include <iostream>
#include <imgui_stdlib.h>

template <typename T>
static inline T remap(T x, T x0, T x1, T y0, T y1)
{
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}



#ifndef IMPLOT_DISABLE_OBSOLETE_FUNCTIONS
#define IMPLOT_DISABLE_OBSOLETE_FUNCTIONS
#endif

#include "implot.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _MSC_VER
#define sprintf sprintf_s
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define CHECKBOX_FLAG(flags, flag) ImGui::CheckboxFlags(#flag, (unsigned int*)&flags, flag)

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


void Demo_DigitalPlots() {
    ImGui::BulletText("Digital plots do not respond to Y drag and zoom, so that");
    ImGui::Indent();
    ImGui::Text("you can drag analog plots over the rising/falling digital edge.");
    ImGui::Unindent();

    static bool paused = false;
    static ScrollingBuffer dataDigital[2];
    static ScrollingBuffer dataAnalog[2];
    static bool showDigital[2] = { true, false };
    static bool showAnalog[2] = { true, false };

    char label[32];
    ImGui::Checkbox("digital_0", &showDigital[0]); ImGui::SameLine();
    ImGui::Checkbox("digital_1", &showDigital[1]); ImGui::SameLine();
    ImGui::Checkbox("analog_0", &showAnalog[0]);  ImGui::SameLine();
    ImGui::Checkbox("analog_1", &showAnalog[1]);

    static float t = 0;
    if (!paused) {
        t += ImGui::GetIO().DeltaTime;
        //digital signal values
        if (showDigital[0])
            dataDigital[0].AddPoint(t, sinf(2 * t) > 0.45);
        if (showDigital[1])
            dataDigital[1].AddPoint(t, sinf(2 * t) < 0.45);
        //Analog signal values
        if (showAnalog[0])
            dataAnalog[0].AddPoint(t, sinf(2 * t));
        if (showAnalog[1])
            dataAnalog[1].AddPoint(t, cosf(2 * t));
    }
    if (ImPlot::BeginPlot("##Digital")) {
        ImPlot::SetupAxisLimits(ImAxis_X1, t - 10.0, t, paused ? ImGuiCond_Once : ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1, 1);
        for (int i = 0; i < 2; ++i) {
            if (showDigital[i] && dataDigital[i].Data.size() > 0) {
                snprintf(label, sizeof(label), "digital_%d", i);
                ImPlot::PlotDigital(label, &dataDigital[i].Data[0].x, &dataDigital[i].Data[0].y, dataDigital[i].Data.size(), 0, dataDigital[i].Offset, 2 * sizeof(float));
            }
        }
        for (int i = 0; i < 2; ++i) {
            if (showAnalog[i]) {
                snprintf(label, sizeof(label), "analog_%d", i);
                if (dataAnalog[i].Data.size() > 0)
                    ImPlot::PlotLine(label, &dataAnalog[i].Data[0].x, &dataAnalog[i].Data[0].y, dataAnalog[i].Data.size(), 0, dataAnalog[i].Offset, 2 * sizeof(float));
            }
        }
        ImPlot::EndPlot();
    }
}

struct ImGraph : App {

    Expression expr;
    ImPlotRect limits;

    using App::App;

    void Start() override {
        expr.set("0.25*sin(2*pi*5*x)+0.5");
        expr.color = ImVec4(1, 0.75f, 0, 1);
    }

    void Update() override 
    {
        Demo_DigitalPlots();
    }
};
int main(int argc, char const* argv[])
{
    ImGraph app("MyWindow!", 640, 480, argc, argv);
    app.Run();
    return 0;
}