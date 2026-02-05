#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

// Includes de Neun (sin prefijos de carpeta, ya que están en /usr/local/Neun/0.4.0)
#include "HindmarshRoseModel.h"
#include "Euler.h"
#include "RungeKutta4.h"
#include "RungeKutta6.h"

// ==========================================
// 1. CONFIGURACIÓN DE PARÁMETROS (MACROS)
// ==========================================
#define HR_I    3.281
#define HR_R    0.0021
#define HR_S    1.0
#define HR_XR   -1.6
#define HR_A    1.0
#define HR_B    3.0
#define HR_C    1.0
#define HR_D    5.0
#define HR_VH   0.1

#define TRANSIENT_TIME 2000.0
#define MEASURE_TIME   5000.0

// ==========================================
// 2. LISTA DE PASOS DE TIEMPO
// ==========================================
const std::vector<double> DTS_TO_TEST = {
    0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 
    0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500, 
    0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500, 
    0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300, 
    0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100, 
    0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600, 
    0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600, 
    0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900, 
    0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017400, 0.018600, 
    0.020000, 0.021600, 0.023500, 0.025700, 0.028300, 0.028600, 0.028900, 0.029300, 
    0.029700, 0.030100, 0.030500, 0.030900, 0.031300, 0.031700, 0.032100, 0.032500, 
    0.032900, 0.033300, 0.033700, 0.034200, 0.034700, 0.035200, 0.035700, 0.036200, 
    0.036700, 0.037200, 0.037700, 0.038300, 0.038900, 0.039500, 0.040100, 0.040700, 
    0.041300, 0.042000, 0.042700, 0.043400, 0.044100, 0.044900, 0.045700, 0.046500, 
    0.047300, 0.048200, 0.049100, 0.050000, 0.051000, 0.052000, 0.053000, 0.054100, 
    0.055200, 0.056400, 0.057600, 0.058900, 0.060200, 0.061600, 0.063000, 0.064500, 
    0.066100, 0.067700, 0.069400, 0.071200, 0.073100, 0.075100, 0.077200, 0.079500, 
    0.081900, 0.084400, 0.087100, 0.090000, 0.093100, 0.096400, 0.100000
};

// ==========================================
// 3. MODELO COMPATIBLE CON CONCEPTOS (FIX)
// ==========================================

// Esta clase hereda de HindmarshRoseModel y añade los métodos get/set
// requeridos por el SystemConcept de Neun, aunque no se usen en la integración por buffers.
template <typename P>
struct CompliantHindmarshRose : public HindmarshRoseModel<P> {
    using typename HindmarshRoseModel<P>::variable;
    using typename HindmarshRoseModel<P>::parameter;
    using typename HindmarshRoseModel<P>::precission_t;

    // Métodos dummy para satisfacer static_assert(SystemConcept)
    precission_t get(variable v) const { return 0.0; }
    precission_t get(parameter p) const { return 0.0; }
    void set(variable v, precission_t val) {}
    void set(parameter p, precission_t val) {}
};

// ==========================================
// 4. LÓGICA DE SIMULACIÓN
// ==========================================

void setup_hr_params(double* params) {
    typedef HindmarshRoseModel<double> HR;
    params[HR::e]  = HR_I;
    params[HR::mu] = HR_R;
    params[HR::S]  = HR_S;
    params[HR::a]  = HR_A;
    params[HR::b]  = HR_B;
    params[HR::c]  = HR_C;
    params[HR::d]  = HR_D;
    params[HR::xr] = HR_XR;
    params[HR::vh] = HR_VH;
}

template <typename TIntegrator>
double measure_points_per_burst(double dt) {
    // USAMOS LA CLASE COMPATIBLE AQUÍ
    CompliantHindmarshRose<double> model;
    
    double vars[HindmarshRoseModel<double>::n_variables];
    double params[HindmarshRoseModel<double>::n_parameters];

    setup_hr_params(params);
    vars[HindmarshRoseModel<double>::x] = -1.6;
    vars[HindmarshRoseModel<double>::y] = -10.0;
    vars[HindmarshRoseModel<double>::z] = 0.0;

    // Transitorio
    int steps_transient = (int)(TRANSIENT_TIME / dt);
    for(int i=0; i<steps_transient; ++i) {
        TIntegrator::step(model, dt, vars, params);
    }

    // Medición
    int steps_measure = (int)(MEASURE_TIME / dt);
    std::vector<double> peak_times;
    double z_prev = vars[HindmarshRoseModel<double>::z];
    double z_curr = z_prev; 
    
    TIntegrator::step(model, dt, vars, params);
    z_curr = vars[HindmarshRoseModel<double>::z];

    for(int i=0; i<steps_measure; ++i) {
        double t_now = i * dt;
        double z_old_prev = z_prev;
        double z_old_curr = z_curr;

        TIntegrator::step(model, dt, vars, params);
        double z_next = vars[HindmarshRoseModel<double>::z];

        if (z_old_curr > z_old_prev && z_old_curr > z_next) {
            peak_times.push_back(t_now);
        }

        z_prev = z_old_curr;
        z_curr = z_next;
        
        if (peak_times.size() >= 6) break; 
    }

    if (peak_times.size() < 2) {
        std::cerr << "Warning: No periodic bursting detected for dt=" << dt << std::endl;
        return -1.0;
    }

    double sum_periods = 0;
    for (size_t k = 1; k < peak_times.size(); ++k) {
        sum_periods += (peak_times[k] - peak_times[k-1]);
    }
    double avg_period = sum_periods / (peak_times.size() - 1);

    return avg_period / dt;
}

void print_cpp_vector(const std::string& name, const std::vector<double>& v) {
    std::cout << "const std::vector<double> " << name << " = {";
    std::cout << std::fixed << std::setprecision(6);
    for(size_t i=0; i<v.size(); ++i) {
        std::cout << v[i];
        if (i < v.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl << std::endl;
}

int main() {
    std::cout << "// =========================================================" << std::endl;
    std::cout << "// TABLAS GENERADAS PARA MODELO HINDMARSH-ROSE (Neun 0.4.0)" << std::endl;
    std::cout << "// =========================================================" << std::endl << std::endl;

    // EULER
    std::vector<double> pts_euler;
    std::cout << "Generando tabla Euler..." << std::endl;
    for(double dt : DTS_TO_TEST) {
        double pts = measure_points_per_burst<Euler>(dt);
        pts_euler.push_back(pts);
    }
    print_cpp_vector("dts_euler", DTS_TO_TEST);
    print_cpp_vector("pts_euler", pts_euler);

    // RUNGE-KUTTA 4
    std::vector<double> pts_rk4;
    std::cout << "Generando tabla RK4..." << std::endl;
    for(double dt : DTS_TO_TEST) {
        double pts = measure_points_per_burst<RungeKutta4>(dt);
        pts_rk4.push_back(pts);
    }
    print_cpp_vector("dts_rk4", DTS_TO_TEST);
    print_cpp_vector("pts_rk4", pts_rk4);

    // RUNGE-KUTTA 6
    std::vector<double> pts_rk6;
    std::cout << "Generando tabla RK6..." << std::endl;
    for(double dt : DTS_TO_TEST) {
        double pts = measure_points_per_burst<RungeKutta6>(dt);
        pts_rk6.push_back(pts);
    }
    print_cpp_vector("dts_rk65", DTS_TO_TEST);
    print_cpp_vector("pts_rk65", pts_rk6);

    return 0;
}