#include "traceRays.h"

// Output stream operator for the structs for the position and direction in
// spherical coordinates
std::ostream& operator<<(std::ostream& os, const DirectionSph dir)
{
    os << "Theta = " << dir.theta << ", Phi = " << dir.phi;
    return os;
}
std::ostream& operator<<(std::ostream& os, const PositionSph pos)
{
    os << "r = " << pos.r << ", Theta = " << pos.theta << ", Phi = " << pos.phi;
    return os;
}

// Function to compute the angles for a point in spherical coordinates
DirectionSph computeAnglesFromxyz( double x, double y, double z )
{
    // Compute the radius
    double r = std::sqrt(x*x + y*y + z*z);
    // Compute the angle theta
    double theta = std::acos(z / r);
    // Compute the angle phi
    double phi = 0;
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
    else    // Case x == y == 0
        phi = 0;

    return DirectionSph( theta, phi );
}

// Function to get the color value of an image in a given angle
PixelRGB mapDirectionToImage( const ImageRGB& image, const DirectionSph& dir )
{
    // Map the angle phi to the range (0, 2 Pi)
    double phi = std::fmod( dir.phi, 2. * M_PI );
    if (phi < 0.)
        phi = 2. * M_PI + phi;
    if ( std::isnan(phi) )
        phi = 0.;

    // Map the angle theta to the range (0, Pi)
    double theta = std::fmod( dir.theta, 2. * M_PI );
    if (theta < 0.)
        theta = 2. * M_PI + theta;
    if (theta > M_PI)
    {
        // If needed, the angle phi should be converted too
        // I think this never happens!
        theta -= M_PI;
        phi += M_PI;
        phi = std::fmod( phi, 2.*M_PI );
    }

    // Map the horizontal direction (0, w) to (0, 2 Pi) in phi
    int ix = (int)( phi / ( 2. * M_PI ) * (image.width - 1) );
    // Map the vertical direction (0, h) to (0, Pi) in theta
    int iy = (int)( theta / M_PI * (image.height - 1) );

    return image.pixels[ image.width * iy + ix ];
}

//=====================================================
//  Propagate a ray backwards from the camera, in a given direction
//=====================================================
PositionSph traceRayBack(PositionSph cameraPos, DirectionSph viewDir, DirectionSph rayDir)
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
    // The conversions are written in page 193 of my notes
    double n_theta = - std::sin( M_PI - rayDir.theta + viewDir.theta ) * std::sin( rayDir.phi - viewDir.phi );
    double n_phi   = - std::sin( M_PI - rayDir.theta + viewDir.theta ) * std::cos( rayDir.phi - viewDir.phi );
    // Convert these to the initial values of the 4-momentum of the ray
    double p_theta = cameraPos.r * n_theta;
    double p_phi   = cameraPos.r * std::sin(cameraPos.theta) * n_phi;
    // Get the initial value of p_r from the two above
    double p_r = std::sqrt( 1. - ( (n_theta * n_theta) + (n_phi * n_phi) ) / (1. - 2./cameraPos.r) );

    // Solve the differential equations backwards in time
    int n_steps = 2000;
    std::vector<std::vector<double>> solution = solveRungeKutta4th5eqCustom(&dr_dlambda, &dp_r_dlambda, &dtheta_dlambda, &dp_theta_dlambda, &dphi_dlambda,
                                                                           n_steps, 0., -200., cameraPos.r, p_r, cameraPos.theta, p_theta, cameraPos.phi, p_phi);

    // float thr = 0.001;
    // if (fabs(rayDir.phi) < thr && rayDir.theta < M_PI + thr && rayDir.theta > M_PI - thr)
    // {
    //     std::cout << "algo\n";
    //     // std::vector<std::vector<double>> output { solution[0], solution[1], solution[3], solution[5] };
    //     // outputDataToFile(output, "out.csv");
    //     outputDataToFile(solution, "out.csv");
    // }

    // Return the final value of the ray angle
    // The angles are the variables solution[2] and solution[4]
    return PositionSph( solution[1][n_steps], solution[3][n_steps], solution[5][n_steps] );
}

//=====================================================
//  Propagate all the rays backwards in a given frame
//=====================================================
void computeFrame(PositionSph cameraPos, DirectionSph viewDir, int width, int height, double thetaFov, ImageRGB background)
{
    // Compute z_0 (distance from the camera point to the screen) from the field of view.
    // The field of view angle is given in degrees.
    double z_0 = width / ( 2. * std::tan(thetaFov / 2. * M_PI / 180.) );

    // Vector to store the angles obtained as a result of the propagation of the rays
    std::vector<double> anglesaux ( width * height * 2, 0. );
    std::vector<double> raux ( width * height , 0. );
    std::vector <DirectionSph> angles (width * height);

    // Iterate through the entire frame
    for (int ix = 0; ix < width; ++ix)
    {
        for (int iy = 0; iy < height; ++iy)
        {
            // Compute the x and y coordinates of the pixel
            double x = -width / 2. + (double)ix; 
            double y = -height / 2. + (double)iy; 

            // Get the angles from these x, y and z coordinates
            DirectionSph rayDir = computeAnglesFromxyz( x, y, -z_0 );

            // Propagate the rays
            PositionSph resultPos = traceRayBack( cameraPos, viewDir, rayDir );
            // std::cout << resultPos << std::endl;

            // Store these in the auxiliary vectors
            raux[ width * iy + ix ] = resultPos.r;
            anglesaux[ 2 * (width * iy + ix) ] = resultPos.theta;
            anglesaux[ 2 * (width * iy + ix) + 1 ] = resultPos.phi;

            angles[ width * iy + ix ] = DirectionSph(resultPos.theta, resultPos.phi);
        }
        if (ix % 100 == 0) std::cout << ix << " / " << width << "\n";
    }

    // =========================================================================
    // Make an image from the background
    // =========================================================================
    ImageRGB outImage(width, height);
    for (int i = 0; i < width * height; ++i)
    {
        outImage.pixels[i] = mapDirectionToImage(background, angles[i]);
        // If the ray has gone inside the horizon, the pixels must be black
        if (raux[i] < 2.)
            // outImage.pixels[i] = PixelRGB(0, 255, 0);
            outImage.pixels[i] = PixelRGB(0, 0, 0);
    }
    // Write the image to a file
    saveToPpm("out.ppm", outImage);

    // =========================================================================
    // Make an image with the values of the angles
    // =========================================================================
    // Convert the angles to the appropriate range
    std::vector<int> colors ( width * height * 3, 0 );
    for (int i = 0; i < width * height; ++i)
    {
        if (raux[i] < 2.)
        {
            // std::cout << "dentro\n";
            colors[3 * i] = 255;
            colors[3 * i + 1] = 0;
        }
        else
        {
            // Map the angle phi to the range (0, 2 Pi)
            double phi = std::fmod( anglesaux[2*i+1], 2. * M_PI );
            if (phi < 0.)
                phi = 2. * M_PI + phi;
            if ( std::isnan(phi) )
                phi = 0.;
            colors[3 * i] = (int)( anglesaux[2 * i] / M_PI * 255. );
            colors[3 * i + 1] = (int)(( phi ) / ( 2. * M_PI ) * 255. );
        }
    }
    // Store these values in an image
    saveToPpm("outangles.ppm", colors, width, height, 3);
}
