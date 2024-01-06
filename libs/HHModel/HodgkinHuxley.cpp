#include <vector>
#include <cmath>
#include <algorithm>
//#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

//class HodgkinHuxley {
//
//public:
//    /* main model */
//    HodgkinHuxley(
//        const std::vector<double>& t, // constant while we add continuous feature
//        const std::vector<double>& inj_current,
//        double C_m = 1.0, 
//        double g_Na = 120.0, 
//        double g_K = 36.0, 
//        double g_L = 0.3,
//        double E_Na = 50.0, 
//        double E_K = -77, 
//        double E_L = -54.387): t(t), inj_current(inj_current), C_m(C_m), g_Na(g_Na), g_K(g_K), g_L(g_L),
//        E_Na(E_Na), E_K(E_K), E_L(E_L) {}
//    
//    // private methods and members 
//private:
//
//    /* main variables used in model */
//    std::vector<double> t;
//    std::vector<double> inj_current;
//    double C_m, g_Na, g_K, g_L, E_Na, E_K, E_L;
//
//    /* main methods */
//
//    double alpha_m(double V) {
//        return 0.1 * (V + 40.0) / (1.0 - exp(-(V + 40.0) / 10.0));
//    }
//
//    double beta_m(double V) {
//        return 4.0 * exp(-(V + 65.0) / 18.0);
//    }
//
//    double alpha_h(double V) {
//        // Channel gating kinetics. Functions of membrane voltage
//        return 0.07 * exp(-(V + 65.0) / 20.0);
//    }
//
//    double beta_h(double V) {
//        // Channel gating kinetics. Functions of membrane voltage
//        return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
//    }
//
//    double alpha_n(double V) {
//        //"""Channel gating kinetics. Functions of membrane voltage"""
//        return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
//    }
//
//    double beta_n(double V) {
//        // Channel gating kinetics. Functions of membrane voltage"""
//        return 0.125 * exp(-(V + 65) / 80.0);
//    }
//
//    double I_Na(double V, double m, double h) {
//        /* Sodium membrane current (in uA/cm^2) */
//        return g_Na * pow(m, 3) * h * (V - E_Na);
//    }
//
//    double I_K(double V, double n) {
//        /* Potassium membrane current (in uA/cm^2) */
//        return g_K * pow(n, 4) * (V - E_K);
//    }
//
//    double I_L(double V) {
//        /* Leak membrane current (in uA/cm^2) */
//        return g_L * (V - E_L);
//    }
//
//    double I_inj(double t) {
//        // Implementation depends on how 'inj_current' is defined
//        auto it = std::lower_bound(this->t.begin(), this->t.end(), t);
//        if (it != this->t.end()) {
//            size_t index = std::distance(this->t.begin(), it);
//            return inj_current[index];
//        }
//        return 0.0;
//    }
//
//    std::vector<double> dALLdt(std::vector<double> X, t) {
//       /* Integrate */
//
//        /* docs https://www.boost.org/doc/libs/1_82_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/integrate_functions.html*/
//        /*V, m, h, n = X
//
//        dVdt = (self.I_inj(t) - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V)) / self.C_m
//        dmdt = self.alpha_m(V) * (1.0 - m) - self.beta_m(V) * m
//        dhdt = self.alpha_h(V) * (1.0 - h) - self.beta_h(V) * h
//        dndt = self.alpha_n(V) * (1.0 - n) - self.beta_n(V) * n
//        return dVdt, dmdt, dhdt, dndt*/
//
//        
//    }
//
//    std::vector<double> results() {
//        /* Main demo for the Hodgkin Huxley neuron model */
//        
//        /* solve ode */
//        // X = odeint(self.dALLdt, [-65, 0.05, 0.6, 0.32], self.t, args=(self,))
//        
//        /*
//            V = X[:,0]
//            m = X[:,1]
//            h = X[:,2]
//            n = X[:,3]
//            i_na = self.I_Na(V, m, h)
//            i_k = self.I_K(V, n)
//            i_l = self.I_L(V)
//        */
//            
//    }
//};
//
//// Helper functions for injection_current class
//class InjectionCurrent {
//public:
//    static std::vector<double> step(const std::vector<double>& t, double t0, double duration, double max = 10) {
//        std::vector<double> out(t.size(), 0);
//        for (std::size_t i = 0; i < t.size(); ++i) {
//            if (t[i] >= t0 && t[i] <= t0 + duration) {
//                out[i] = max;
//            }
//        }
//        return out;
//    }
//};
//
//class Helpers {
//public: 
//    std::vector<double> CreateTimeVector(double start, double end, int size) {
//        /* create a time vector of size "size" of doubles. */
//        std::vector<double> t_vec(size);
//        for (int i = 0; i < size; ++i) {
//            t_vec[i] = start + 1.0;
//        }
//        return t_vec;
//    }
//};

int main() {

    
    /* preamble setup to run model */
    /* create time vector */
    //Helpers Helper;
    //std::vector<double> t = Helper.CreateTimeVector(0.0, 500.0, 500);

    typedef boost::array< double, 3 > state_type;


    /* run main model */
    //HodgkinHuxley Model(t);
}
