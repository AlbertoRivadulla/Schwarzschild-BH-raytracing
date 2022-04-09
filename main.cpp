#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

#include "saveToFiles.h"
#include "rungeKutta.h"
#include "traceRays.h"

// define the functions that appear in the differential equations
float f1( float y1, float y2, float t )
{
    return y2;
}

float f2( float y1, float y2, float t )
{
    return - 0.01 * y2 - 9.81;
}

int main()
{
    // Get the time of the beginning of the execution
    // The return type of this is std::chrono::time_point<std::chrono::high_resolution_clock>
    auto t0 {std::chrono::high_resolution_clock::now()};

    // // Solve the system of differential equations
    // std::vector<std::vector<float>> solution = solveRungeKutta4th2eq(&f1, &f2, 200, 0., 1000., 10., 600. * std::sin(60. * M_PI / 180.));
    // // Save the data to a file
    // outputDataToFile(solution, "output.csv");

    // Compute a single frame
    computeFrame(PositionSph( 2., M_PI / 2., 0. ), DirectionSph(0., 0.), 200, 100, 45.);

    // Get the time of the end of the execution
    auto t1 {std::chrono::high_resolution_clock::now()};
    // Compute the difference between the two times
    // Get the number of miliseconds as a double
    std::chrono::duration<double, std::milli> duration {t1 - t0};
    // In order to get the duration in seconds, write simply
    // std::chrono::duration<double> ... ;
    // Print the duration of execution
    std::cout << "Duration of execution: " << duration.count() << " ms\n";

    return 0;
}
