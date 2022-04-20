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

// Function to get the color value of an image in a given angle
PixelRGB mapDirectionToImage( const ImageRGB& image, const DirectionSph& dir )
{
    // Map the angle phi to the range (0, 2 Pi)
    float phi = std::fmod( dir.phi, 2. * M_PI );
    if (phi < 0.)
    {
        phi = 2. * M_PI + phi;
    }
    if ( std::isnan(phi) )
    {
        phi = 0.;
    }

    // Map the angle theta to the range (0, Pi)
    // float theta = dir.theta > 0. ? dir.theta : 0.;
    // if (theta > M_PI) theta = 0.;
    float theta = std::fmod( dir.theta, 2. * M_PI );
    if (theta < 0.)
    {
        theta = 2. * M_PI + theta;
    }
    if (theta > M_PI)
    {
        theta -= M_PI;
        phi += M_PI;
        phi = std::fmod( phi, 2.*M_PI );
        // std::cout << "====== " << theta << ' ' << phi << std::endl;
    }

    // Map the horizontal direction (0, w) to (0, 2 Pi) in phi
    // int ix = (int)( dir.phi / ( 2. * M_PI ) * (image.width - 1) );
    int ix = (int)( phi / ( 2. * M_PI ) * (image.width - 1) );
    // Map the vertical direction (0, h) to (0, Pi) in theta
    // int iy = (int)( dir.theta / M_PI * (image.height - 1) );
    int iy = (int)( theta / M_PI * (image.height - 1) );
    //
    // std::cout << "---" << std::endl;
    // std::cout << phi << ' ' << dir.theta << std::endl;
    // std::cout << ix << ' ' << image.width << ' ' << iy << ' ' << image.height << std::endl;
    // std::cout << image.width * iy + ix << ' ' << image.pixels.size() << std::endl;

    return image.pixels[ image.width * iy + ix ];
    // return image.pixels[ image.width * iy + image.width - ix ];
    // return image.pixels[ image.width * (image.height - iy) + ix ];
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
    // float n_theta = std::cos( rayDir.theta - viewDir.theta );
    // float n_theta = - std::sin( rayDir.theta - viewDir.theta ) * std::sin( rayDir.phi - viewDir.phi );
    // float n_phi   = - std::sin( rayDir.theta - viewDir.theta ) * std::cos( rayDir.phi - viewDir.phi );
    // rayDir.theta = M_PI;
    float n_theta = - std::sin( M_PI - rayDir.theta + viewDir.theta ) * std::sin( rayDir.phi - viewDir.phi );
    float n_phi   = - std::sin( M_PI - rayDir.theta + viewDir.theta ) * std::cos( rayDir.phi - viewDir.phi );
    // Since I propagate the rays backwards, I change the sign of these
    // float p_theta = - cameraPos.r * n_theta;
    // float p_phi   = - cameraPos.r * std::sin(cameraPos.theta) * n_phi;
    float p_theta = cameraPos.r * n_theta;
    float p_phi   = cameraPos.r * std::sin(cameraPos.theta) * n_phi;

    // std::cout << p_theta << " --- " << p_phi << '\n';

    // Get the initial value of p_r from the two above
    float p_r = std::sqrt( 1. - ( (n_theta * n_theta) + (n_phi * n_phi) ) / (1. - 2./cameraPos.r) );
    // float p_r = std::sqrt( 1. - ( (n_theta * n_theta) + (n_phi * n_phi) ) / (1.) );
    // std::cout << "p_phi = ";
    // std::cout << p_phi << '\n';
    // std::cout << p_r << '\n';

    // Solve the differential equations backwards in time
    int n_steps = 20000;
    std::vector<std::vector<float>> solution = solveRungeKutta4th5eqCustom(&dr_dlambda, &dp_r_dlambda, &dtheta_dlambda, &dp_theta_dlambda, &dphi_dlambda,
                                                                           n_steps, 0., -60., cameraPos.r, p_r, cameraPos.theta, p_theta, cameraPos.phi, p_phi);

    // std::cout << solution[2][n_steps-1] << ' ' << solution[4][n_steps-1] << ' ' << std::fmod(solution[4][n_steps-1], 2*M_PI) << '\n';

    // std::vector<std::vector<float>> output { solution[0], solution[1], solution[3], solution[5] };
    // outputDataToFile(output, "out.csv");

    // Return the final value of the ray angle
    // The angles are the variables solution[2] and solution[4]
    return PositionSph( solution[1][n_steps], solution[3][n_steps], solution[5][n_steps] );

    // // return DirectionSph(0., 0.);
    // return DirectionSph( -n_theta, -std::sin(cameraPos.theta) * n_phi );
}

//=====================================================
//  Propagate all the rays backwards in a given frame
//=====================================================
void computeFrame(PositionSph cameraPos, DirectionSph viewDir, int width, int height, float thetaFov, ImageRGB background)
{
    // Compute z_0 (distance from the camera point to the screen) from the field of view.
    // The field of view angle is given in degrees.
    float z_0 = width / ( 2. * std::tan(thetaFov / 2. * M_PI / 180.) );

    // Vector to store the angles obtained as a result of the propagation of the rays
    std::vector<float> anglesaux ( width * height * 2, 0. );
    std::vector<float> raux ( width * height , 0. );
    std::vector <DirectionSph> angles (width * height);

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
            PositionSph resultPos = traceRayBack( cameraPos, viewDir, rayDir );
            // std::cout << resultPos << std::endl;

            // Store these in the auxiliary vectors
            raux[ width * iy + ix ] = resultPos.r;
            anglesaux[ 2 * (width * iy + ix) ] = resultPos.theta;
            anglesaux[ 2 * (width * iy + ix) + 1 ] = resultPos.phi;

            // std::cout << resultDir << '\n';
            angles[ width * iy + ix ] = DirectionSph(resultPos.theta, resultPos.phi);
        }
        if (ix % 100 == 0) std::cout << ix << " / " << width << "\n";
    }

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // Test output

    // Make an image from the background
    ImageRGB outImage(width, height);
    for (int i = 0; i < width * height; ++i)
    {
        // std::cout << angles[i] << '\n';
        // std::cout << raux[i] << '\n';
        outImage.pixels[i] = mapDirectionToImage(background, angles[i]);
        if (raux[i] < 2.)
        {
            outImage.pixels[i] = PixelRGB(0, 255, 0);
        }
    }
    // Write the image to a file
    saveToPpm("out.ppm", outImage);

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
            float phi = std::fmod( anglesaux[2*i+1], 2. * M_PI );
            if (phi < 0.)
            {
                phi = 2. * M_PI + phi;
            }
            if ( std::isnan(phi) )
            {
                std::cout << "phi is nan \n";
                phi = 0.;
            }
            // std::cout << anglesaux[2*i] << ' ' << phi << '\n';
            colors[3 * i] = (int)( anglesaux[2 * i] / M_PI * 255. );
            // colors[3 * i + 1] = (int)(( anglesaux[2 * i + 1]) / ( 2. * M_PI ) * 255. );
            colors[3 * i + 1] = (int)(( phi ) / ( 2. * M_PI ) * 255. );
        }
    }
    // Store these values in an image
    saveToPpm("outangles.ppm", colors, width, height, 3);
}
