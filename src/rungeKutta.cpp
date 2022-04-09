#include "rungeKutta.h"

// =====================================
// Based on section "Runge Kutta for two coupled 1st order differential equations"
//      https://www.phys.ksu.edu/personal/washburn/Teaching/Class%20Files/NQO/Tutorials/Tutorial4_RungeKutta4thOrder.pdf
//
//
// in order to open and plot the output data in Mathematica, use:
/*
table = import["path/out.csv", "Table", "FieldSeparators" -> ";"]
Listplot[table[[;; , ;; 2]], ImageSize -> Large]
Listplot[table[[;; , ;; 3 ;; 2]], ImageSize -> Large]
*/
//      
// to get the path, run "pwd" in the build directory.
//
// =====================================

std::vector<std::vector<float>> solveRungeKutta4th2eq
    (
    float (*f1)(float, float, float),   // Equation for dy1/dt
    float (*f2)(float, float, float),   // Equation for dy2/dt
    int n,                              // Number of steps
    float t0, float tf,                 // Range of values of the independent variable
    float y10, float y20                // Initial conditions
    )
{
    // Compute the step size
    float h = ( tf - t0 ) / (float)n;

    // initialize vectors for the solution of the different functions
    std::vector<float> t ( n + 1, t0 );
    std::vector<float> y1 ( n + 1, y10 );
    std::vector<float> y2 ( n + 1, y20 );
    // values of the independent variable t
    for (int i = 0; i <= n; ++i)
    {
        t[i] = t0 + i * h;
    }

    // solve the system of differential equations
    for (int i = 0; i < n; ++i)
    {
        // compute the different parameters in the method
        float k1 = f1( y1[i], y2[i], t[i] );
        float m1 = f2( y1[i], y2[i], t[i] );
        float k2 = f1( y1[i] + 0.5 * k1 * h, y2[i] + 0.5 * m1 * h, t[i] + 0.5 * h );
        float m2 = f2( y1[i] + 0.5 * k1 * h, y2[i] + 0.5 * m1 * h, t[i] + 0.5 * h );
        float k3 = f1( y1[i] + 0.5 * k2 * h, y2[i] + 0.5 * m2 * h, t[i] + 0.5 * h );
        float m3 = f2( y1[i] + 0.5 * k2 * h, y2[i] + 0.5 * m2 * h, t[i] + 0.5 * h );
        float k4 = f1( y1[i] + k3 * h, y2[i] + m3 * h, t[i] + h );
        float m4 = f2( y1[i] + k3 * h, y2[i] + m3 * h, t[i] + h );

        // compute the solution for the next point
        t[i+1] = t[i] + h;
        y1[i+1] = y1[i] + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6;
        y2[i+1] = y2[i] + (m1 + 2. * m2 + 2. * m3 + m4) * h / 6;
    }

    std::vector<std::vector<float>> solution { t, y1, y2 };
    return solution;
}

//=====================================================
//  Custom implementation of the 4th order Runge-Kutta method for my problem with
//  5 equations
//=====================================================

std::vector<std::vector<float>> solveRungeKutta4th5eqCustom
    (
    float (*f1)(float, float, float),   // Equation for dy1/dt
    float (*f2)(float, float, float),   // Equation for dy2/dt
    float (*f3)(float, float, float),   // Equation for dy3/dt
    float (*f4)(float, float, float),   // Equation for dy4/dt
    float (*f5)(float, float, float),   // Equation for dy5/dt
    int n,                              // Number of steps
    float lambda0, float lambdaf,       // Range of values of the independent variable
    float y10, float y20, float y30,    // Initial conditions
    float y40, float y50,
    float p_phi                         // Value of p_phi (constant)
    );

