#include "traceRays.h"
#include "saveToFiles.h"

// Function to compute the angles for a point in spherical coordinates
DirectionSph computeAnglesFromxyz( float x, float y, float z )
{
    // Compute the radius
    float r = std::sqrt(x*x + y*y + z*z);

    // Compute the angle theta
    float theta = std::acos(z / r);

    // Compute the angle phi
    float phi = 0;
    if ( x > 0 )
        phi = atan( y / x );
    else if ( x < 0 && y >= 0 )
        phi = atan( y / x ) + M_PI;
    else if ( x < 0 && y < 0 )
        phi = atan( y / x ) - M_PI;
    else if ( x == 0 && y > 0 )
        phi = M_PI / 2;
    else if ( x == 0 && y < 0 )
        phi = - M_PI / 2;
    // else    // Case x == y == 0
    //     phi = 0;

    return DirectionSph( theta, phi );
}

//=====================================================
//  Propagate a ray backwards from the camera, in a given direction
//=====================================================
DirectionSph traceRayBack(PositionSph cameraPos, DirectionSph viewDir, DirectionSph rayDir)
/* Input:
   - Position of the camera (only r).
   - View direction of the camera, relative to the angular direction.
   - Direction of the incoming ray relative to the view direction of the camera.
*/
{
    // I still have to take into account the viewDir of the camera here!!
    // These should be some offsets for the direction angles of the ray with respect
    // to the camera.


    // Get the initial values of p_theta and p_phi from the direction of the ray
    // The conversions are writte in page 193 of my notes
    // float n_theta = std::cos( rayDir.theta - viewDir.theta );
    float n_theta = - std::sin( rayDir.theta - viewDir.theta ) * std::sin( rayDir.phi - viewDir.phi );
    // float n_phi   = - std::sin( rayDir.theta - viewDir.theta ) * std::sin( rayDir.phi - viewDir.phi );
    float n_phi   = - std::sin( rayDir.theta - viewDir.theta ) * std::cos( rayDir.phi - viewDir.phi );
    // Since I propagate the rays backwards, I change the sign of these
    float p_theta = - cameraPos.r * n_theta;
    float p_phi   = - cameraPos.r * std::sin(cameraPos.theta) * n_phi;

    // Get the initial value of p_r from the two above
    float p_r = std::sqrt( 1. - ( (n_theta * n_theta) + (n_phi * n_phi) ) / (1. - 2./cameraPos.r) );

    // Solve the differential equations backwards in time

    // Return the final value of the ray angle

    // return DirectionSph(0., 0.);
    return DirectionSph( -n_theta, -std::sin(cameraPos.theta) * n_phi );
}

//=====================================================
//  Propagate all the rays backwards in a given frame
//=====================================================
void computeFrame(PositionSph cameraPos, DirectionSph viewDir, int width, int height, float thetaFov)
{
    // Compute z_0 (distance from the camera point to the screen) from the field of view.
    // The field of view angle is given in degrees.
    float z_0 = width / ( 2. * std::tan(thetaFov / 2. * M_PI / 180.) );

    // Vector to store the angles obtained as a result of the propagation of the rays
    std::vector<float> angles ( width * height * 2, 0 );

    // Iterate through the entire frame
    for (int ix = 0; ix < width; ++ix)
    {
        for (int iy = 0; iy < height; ++iy)
        {
            // Compute the x and y coordinates of the pixel
            float x = -width / 2. + (float)ix; 
            float y = -height / 2. + (float)iy; 

            // Get the angles from these x, y and z coordinates
            DirectionSph rayDir = computeAnglesFromxyz( x, y, -z_0 );

            // Propagate the rays
            DirectionSph resultDir = traceRayBack( cameraPos, viewDir, rayDir );

            // Store these in the auxiliary vector
            angles[ 2 * (width * iy + ix) ] = resultDir.theta;
            angles[ 2 * (width * iy + ix) + 1 ] = resultDir.phi;
        }
    }

    /////////////////////////////////////////////////////////////
    // Test output

    // Convert the angles to the appropriate range
    std::vector<int> colors ( width * height * 3, 0 );
    for (int i = 0; i < width * height; ++i)
    {
        colors[3 * i] = (int)( angles[2 * i] / M_PI * 255. );
        colors[3 * i + 1] = (int)(( angles[2 * i + 1] + M_PI ) / ( 2. * M_PI ) * 255. );
    }

    // Store these values in an image
    saveToPpm("out.ppm", colors, width, height, 3);

    // for (int i = 0; i < width * height; ++i)
    // {
    //     std::cout << angles[2*i] << ' ' << angles[2*i+1] << '\n';
    //     std::cout << colors[3*i] << ' ' << colors[3*i+1] << ' ' << colors[ 3*i+2 ] << '\n';
    // }
}
