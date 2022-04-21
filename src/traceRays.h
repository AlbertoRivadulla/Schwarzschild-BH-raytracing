#ifndef TRACE_RAYS_H
#define TRACE_RAYS_H

#include <cmath>
#include "rungeKutta.h"
#include "ioFiles.h"

struct DirectionSph
{
    double theta;
    double phi;

    DirectionSph(double thetaVal, double phiVal) :
        theta { thetaVal }, phi { phiVal }
    {}
    DirectionSph() : theta { 0 }, phi { 0 }
    {}
};

struct PositionSph
{
    double r;
    double theta;
    double phi;
    // DirectionSph direction;

    PositionSph(double rVal, double thetaVal, double phiVal) : 
        // r { rVal }, direction { DirectionSph(thetaVal, phiVal) }
        r { rVal }, theta { thetaVal }, phi { phiVal }
    {}
};

std::ostream& operator<<(std::ostream& os, const DirectionSph dir);
std::ostream& operator<<(std::ostream& os, const PositionSph pos);

// Function to compute the angles for a point in spherical coordinates
DirectionSph computeAnglesFromxyz( double x, double y, double z );

// Function to get the color value of an image in a given angle
PixelRGB mapDirectionToImage( const ImageRGB& image, const DirectionSph& dir );

//=====================================================
//  Propagate a ray backwards from the camera, in a given direction
//=====================================================
PositionSph traceRayBack(PositionSph cameraPos, DirectionSph viewDir, DirectionSph rayDir);
/* Input:
   - Position of the camera (only r).
   - View direction of the camera, relative to the angular direction.
   - Direction of the incoming ray relative to the view direction of the camera.
*/

//=====================================================
//  Propagate all the rays backwards in a given frame
//=====================================================
void computeFrame(PositionSph cameraPos, DirectionSph viewDir, int width, int height, double thetaFov, ImageRGB background);

//=====================================================
//  Equations that determine the evolution of a light ray
//=====================================================
// I always set M = 1, so all dimensional values are measured in units of M
/* 
   - Independent variable:
        - Affine parameter lambda
   - Dependent variables:
        - r
        - p_r
        - theta
        - p_theta
        - phi
        - p_phi (this is constant)
*/

// Equation for dr/dlambda
inline double dr_dlambda( double t, double r, double p_r, double theta, double p_theta, double phi, double p_phi )
{
    return ( 1. - 2. / r ) * p_r;
    // return p_r;
    // return - ( 1. - 2. / r ) * p_r;
    // return -1. * p_r;
}
// Equation for dp_r/dlambda
inline double dp_r_dlambda( double t, double r, double p_r, double theta, double p_theta, double phi, double p_phi )
{
    double sintheta = std::sin(theta);
    return ( p_theta * p_theta + p_phi * p_phi / (sintheta * sintheta) ) / (r * r * r) - 1. / ( r * r * (1-2./r) * (1-2./r) ) - p_r * p_r / (r * r);
    // return ( p_theta * p_theta + p_phi * p_phi / (sintheta * sintheta) ) / (r * r * r);
}
// Equation for dtheta/dlambda
inline double dtheta_dlambda( double t, double r, double p_r, double theta, double p_theta, double phi, double p_phi )
{
    return p_theta / (r * r);
    // return - p_theta / (r * r);
}
// Equation for dp_theta/dlambda
inline double dp_theta_dlambda( double t, double r, double p_r, double theta, double p_theta, double phi, double p_phi )
{
    double sintheta = std::sin(theta);
    return std::cos(theta) / (r * r * sintheta * sintheta * sintheta) * p_phi * p_phi;
    // return - std::cos(theta) / (r * r * sintheta * sintheta * sintheta) * p_phi * p_phi;
}
// Equation for dphi/dlambda
inline double dphi_dlambda( double t, double r, double p_r, double theta, double p_theta, double phi, double p_phi )
{
    double sintheta = std::sin(theta);
    return p_phi / ( r * r * sintheta * sintheta );
    // return - p_phi / ( r * r * sintheta * sintheta );
}
// Equation for dp_phi/dlambda
// This quantity is a constant of movement!
inline double dp_phi_dlambda( double t, double r, double p_r, double theta, double p_theta, double phi, double p_phi )
{
    return 0.;
}

#endif
