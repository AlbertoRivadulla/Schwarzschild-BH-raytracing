#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

// =====================================
// Based on section "Runge Kutta for two coupled 1st order differential equations"
//      https://www.phys.ksu.edu/personal/washburn/Teaching/Class%20Files/NQO/Tutorials/Tutorial4_RungeKutta4thOrder.pdf
//
//
// In order to open and plot the output data in Mathematica, export it to a file (example: out.csv) and use:
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
    );

// std::vector<std::vector<float>> solveRungeKutta4th5eq
//     (
//     float (*f1)(float, float, float),   // Equation for dy1/dt
//     float (*f2)(float, float, float),   // Equation for dy2/dt
//     float (*f3)(float, float, float),   // Equation for dy3/dt
//     float (*f4)(float, float, float),   // Equation for dy4/dt
//     float (*f5)(float, float, float),   // Equation for dy5/dt
//     int n,                              // Number of steps
//     float lambda0, float lambdaf,       // Range of values of the independent variable
//     float y10, float y20, float y30,    // Initial conditions
//     float y40, float y50   
//     );

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

#endif
