#pragma once
#include <iostream>
#include <exprtk.hpp>
#include "imgui.h"

class GatingFunctions { public: 

	static double alpha_n(double V);
	static double beta_n(double V);
	static double alpha_m(double V);
	static double beta_m(double V);
	static double alpha_h(double V);
	static double beta_h(double V);
};