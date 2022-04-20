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
    float (*f1)(float, float, float, float, float, float, float),   // Equation for dy1/dt
    float (*f2)(float, float, float, float, float, float, float),   // Equation for dy2/dt
    float (*f3)(float, float, float, float, float, float, float),   // Equation for dy3/dt
    float (*f4)(float, float, float, float, float, float, float),   // Equation for dy4/dt
    float (*f5)(float, float, float, float, float, float, float),   // Equation for dy5/dt
    int n,                              // Number of steps
    float t0, float tf,                 // Range of values of the independent variable
    float y10, float y20, float y30,    // Initial conditions
    float y40, float y50,
    float p_phi                         // Value of p_phi (constant)
    )
{
    // Compute the step size
    float h = ( tf - t0 ) / (float)n;

    // Initialize vectors for the solution of the different functions
    std::vector<float> t ( n + 1, t0 );
    std::vector<float> y1 ( n + 1, y10 );
    std::vector<float> y2 ( n + 1, y20 );
    std::vector<float> y3 ( n + 1, y30 );
    std::vector<float> y4 ( n + 1, y40 );
    std::vector<float> y5 ( n + 1, y50 );
    // Values of the independent variable t
    for (int i = 0; i <= n; ++i)
    {
        t[i] = t0 + i * h;
    }

    float minr = 1000.;

    // Solve the system of differential equations
    for (int i = 0; i < n; ++i)
    {
        // Compute the different parameters in the method
        float k1 = f1( t[i], y1[i], y2[i], y3[i], y4[i], y5[i], p_phi );
        float m1 = f2( t[i], y1[i], y2[i], y3[i], y4[i], y5[i], p_phi );
        float n1 = f3( t[i], y1[i], y2[i], y3[i], y4[i], y5[i], p_phi );
        float o1 = f4( t[i], y1[i], y2[i], y3[i], y4[i], y5[i], p_phi );
        float p1 = f5( t[i], y1[i], y2[i], y3[i], y4[i], y5[i], p_phi );

        float k2 = f1( t[i]+0.5*h, y1[i]+0.5*k1*h, y2[i]+0.5*m1*h, y3[i]+0.5*n1*h, y4[i]+0.5*o1*h, y5[i]+0.5*p1*h, p_phi );
        float m2 = f2( t[i]+0.5*h, y1[i]+0.5*k1*h, y2[i]+0.5*m1*h, y3[i]+0.5*n1*h, y4[i]+0.5*o1*h, y5[i]+0.5*p1*h, p_phi );
        float n2 = f3( t[i]+0.5*h, y1[i]+0.5*k1*h, y2[i]+0.5*m1*h, y3[i]+0.5*n1*h, y4[i]+0.5*o1*h, y5[i]+0.5*p1*h, p_phi );
        float o2 = f4( t[i]+0.5*h, y1[i]+0.5*k1*h, y2[i]+0.5*m1*h, y3[i]+0.5*n1*h, y4[i]+0.5*o1*h, y5[i]+0.5*p1*h, p_phi );
        float p2 = f5( t[i]+0.5*h, y1[i]+0.5*k1*h, y2[i]+0.5*m1*h, y3[i]+0.5*n1*h, y4[i]+0.5*o1*h, y5[i]+0.5*p1*h, p_phi );

        float k3 = f1( t[i]+0.5*h, y1[i]+0.5*k2*h, y2[i]+0.5*m2*h, y3[i]+0.5*n2*h, y4[i]+0.5*o2*h, y5[i]+0.5*p2*h, p_phi );
        float m3 = f2( t[i]+0.5*h, y1[i]+0.5*k2*h, y2[i]+0.5*m2*h, y3[i]+0.5*n2*h, y4[i]+0.5*o2*h, y5[i]+0.5*p2*h, p_phi );
        float n3 = f3( t[i]+0.5*h, y1[i]+0.5*k2*h, y2[i]+0.5*m2*h, y3[i]+0.5*n2*h, y4[i]+0.5*o2*h, y5[i]+0.5*p2*h, p_phi );
        float o3 = f4( t[i]+0.5*h, y1[i]+0.5*k2*h, y2[i]+0.5*m2*h, y3[i]+0.5*n2*h, y4[i]+0.5*o2*h, y5[i]+0.5*p2*h, p_phi );
        float p3 = f5( t[i]+0.5*h, y1[i]+0.5*k2*h, y2[i]+0.5*m2*h, y3[i]+0.5*n2*h, y4[i]+0.5*o2*h, y5[i]+0.5*p2*h, p_phi );

        float k4 = f1( t[i]+h, y1[i]+k3*h, y2[i]+m3*h, y3[i]+n3*h, y4[i]+o3*h, y5[i]+p3*h, p_phi );
        float m4 = f2( t[i]+h, y1[i]+k3*h, y2[i]+m3*h, y3[i]+n3*h, y4[i]+o3*h, y5[i]+p3*h, p_phi );
        float n4 = f3( t[i]+h, y1[i]+k3*h, y2[i]+m3*h, y3[i]+n3*h, y4[i]+o3*h, y5[i]+p3*h, p_phi );
        float o4 = f4( t[i]+h, y1[i]+k3*h, y2[i]+m3*h, y3[i]+n3*h, y4[i]+o3*h, y5[i]+p3*h, p_phi );
        float p4 = f5( t[i]+h, y1[i]+k3*h, y2[i]+m3*h, y3[i]+n3*h, y4[i]+o3*h, y5[i]+p3*h, p_phi );

        // Compute the solution for the next point
        t[i+1] = t[i] + h;
        y1[i+1] = y1[i] + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6;
        y2[i+1] = y2[i] + (m1 + 2. * m2 + 2. * m3 + m4) * h / 6;
        y3[i+1] = y3[i] + (n1 + 2. * n2 + 2. * n3 + n4) * h / 6;
        y4[i+1] = y4[i] + (o1 + 2. * o2 + 2. * o3 + o4) * h / 6;
        y5[i+1] = y5[i] + (p1 + 2. * p2 + 2. * p3 + p4) * h / 6;

        // std::cout << y3[i+1] << ' ' << y5[i+1] << '\n';

        // std::cout << y1[i+1] << ' ';
        if (y1[i+1] < minr)
            minr = y1[i+1];

        // Check if the ray has gone inside the horizon
        // The coordinate r is y1
        if (y1[i+1] <= 2.)
        {
            // std::cout << "dentro1\n";
            // std::cout << y1[i+1] << '\n';
            // Set the last value of r to be less than 2.
            y1[n] = 1.;
            // Break the loop
            break;
        }

        // Check if the r coordinate is large enough
        if (y1[i+1] > 100.)
        {
            // std::cout << "Large r\n";
            // Set the last values of the angles to be the solutions
            y3[n] = y3[i+1];
            y5[n] = y5[i+1];
            // Break the loop
            break;
        }

        // std::cout << y1[i] << ' ';
        // std::cout << y3[i] << ' ';
        // std::cout << y4[i] << ' ';
        // std::cout << "--- ";
    }

    // std::cout << "minimum r: " << minr << "\n";

    std::vector<std::vector<float>> solution { t, y1, y2, y3, y4, y5 };
    return solution;
}
