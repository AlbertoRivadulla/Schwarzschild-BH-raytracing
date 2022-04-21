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
        - Modify the function traceRayBack to be able to have a viewDir different from (0, 0)
          as not it is not correct.
        - The custom Runge-Kutta solver must check if the ray has gone far enough
          also when it goes around the BH several times. I think in this case the 
          integration ends when the maximum value of the independent variable is
          reached, but this should continue until the ray escapes.
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

    // Compute a single frame
    computeFrame(PositionSph( 40., M_PI / 2., 0. ), DirectionSph(0., 0.), 50, 25, 45., background);

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
