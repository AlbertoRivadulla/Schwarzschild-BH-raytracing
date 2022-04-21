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

std::vector<std::vector<double>> solveRungeKutta4th2eq
    (
    double (*f1)(double, double, double),   // Equation for dy1/dt
    double (*f2)(double, double, double),   // Equation for dy2/dt
    int n,                              // Number of steps
    double t0, double tf,                 // Range of values of the independent variable
    double y10, double y20                // Initial conditions
    );

// std::vector<std::vector<double>> solveRungeKutta4th5eq
//     (
//     double (*f1)(double, double, double),   // Equation for dy1/dt
//     double (*f2)(double, double, double),   // Equation for dy2/dt
//     double (*f3)(double, double, double),   // Equation for dy3/dt
//     double (*f4)(double, double, double),   // Equation for dy4/dt
//     double (*f5)(double, double, double),   // Equation for dy5/dt
//     int n,                              // Number of steps
//     double lambda0, double lambdaf,       // Range of values of the independent variable
//     double y10, double y20, double y30,    // Initial conditions
//     double y40, double y50   
//     );

//=====================================================
//  Custom implementation of the 4th order Runge-Kutta method for my problem with
//  5 equations
//=====================================================

std::vector<std::vector<double>> solveRungeKutta4th5eqCustom
    (
    double (*f1)(double, double, double, double, double, double, double),   // Equation for dy1/dt
    double (*f2)(double, double, double, double, double, double, double),   // Equation for dy2/dt
    double (*f3)(double, double, double, double, double, double, double),   // Equation for dy3/dt
    double (*f4)(double, double, double, double, double, double, double),   // Equation for dy4/dt
    double (*f5)(double, double, double, double, double, double, double),   // Equation for dy5/dt
    int n,                                                                  // Number of steps
    double t0, double tf,                                                   // Range of values of the independent variable
    double y10, double y20, double y30,                                     // Initial conditions
    double y40, double y50,
    double p_phi                                                            // Value of p_phi (constant)
    );

#endif
