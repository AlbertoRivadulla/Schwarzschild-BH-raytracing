#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

#include "ioFiles.h"
#include "rungeKutta.h"
#include "traceRays.h"

//============================================================
/* Notes
        - The input image for the background should be a .ppm file, with a depth of 8 bits.
          To obtain it from an input in.pdf, I use:
                convert in.pdf -depth 8 out.ppm
*/
//============================================================
/*
   Features to implement:
        - Improve mapDirectionToImage() to use linear interpolation between pixels

        - Solve the system of differential equations
*/
//============================================================

int main()
{
    // Get the time of the beginning of the execution
    // The return type of this is std::chrono::time_point<std::chrono::high_resolution_clock>
    auto t0 {std::chrono::high_resolution_clock::now()};

    // Read the data from the file
    // ImageRGB background = loadFromPpm("../resources/fg.ppm");
    ImageRGB background = loadFromPpm("../resources/fg2.ppm");

    // saveToPpm( "test.ppm", image );
    // std::cout << mapDirectionToImage( image, DirectionSph( M_PI/2., M_PI ) );

    // // Solve the system of differential equations
    // std::vector<std::vector<float>> solution = solveRungeKutta4th2eq(&f1, &f2, 200, 0., 1000., 10., 600. * std::sin(60. * M_PI / 180.));
    // // Save the data to a file
    // outputDataToFile(solution, "output.csv");

    // Compute a single frame
    // computeFrame(PositionSph( 40., M_PI / 2., 0. ), DirectionSph(0., 0.), 200, 100, 45., background);
    computeFrame(PositionSph( 40., M_PI / 2., 0. ), DirectionSph(0., 0.), 100, 50, 45., background);
    // computeFrame(PositionSph( 6., M_PI / 2., 0. ), DirectionSph(0., 0.), 20, 10, 45, background);
    // computeFrame(PositionSph( 20., M_PI / 2., M_PI ), DirectionSph(0., 0.), 400, 200, 80., background);

    // int width = 100;
    // int height = 50;
    // float thetaFov = 45.;
    // float z_0 = width / ( 2. * std::tan(thetaFov / 2. * M_PI / 180.) );
    // float x = -width / 2.; 
    // // Compute the x and y coordinates of the pixel
    // float y = -height / 2.; 
    // x = 0;
    // y = 0;
    // // Get the angles from these x, y and z coordinates
    // DirectionSph rayDir = computeAnglesFromxyz( x, y, -z_0 );
    // std::cout << rayDir << '\n';
    // PositionSph resultPos = traceRayBack( PositionSph( 4., M_PI / 2., 0. ), DirectionSph(0., 0.), rayDir );

    // Get the time of the end of the execution
    auto t1 {std::chrono::high_resolution_clock::now()};
    // Compute the difference between the two times
    // Get the number of miliseconds as a double
    // std::chrono::duration<double, std::milli> duration {t1 - t0};
    std::chrono::duration<double> duration {t1 - t0};
    // In order to get the duration in seconds, write
    // std::chrono::duration<double> ... ;
    // Print the duration of execution
    std::cout << "\n\nDuration of execution: " << duration.count() << " s\n";

    return 0;
}
