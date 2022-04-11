#ifndef TRACE_RAYS_H
#define TRACE_RAYS_H

#include <cmath>
#include "rungeKutta.h"
#include "ioFiles.h"

struct DirectionSph
{
    float theta;
    float phi;

    DirectionSph(float thetaVal, float phiVal) :
        theta { thetaVal }, phi { phiVal }
    {}
    DirectionSph() : theta { 0 }, phi { 0 }
    {}
};

struct PositionSph
{
    float r;
    float theta;
    float phi;
    // DirectionSph direction;

    PositionSph(float rVal, float thetaVal, float phiVal) : 
        // r { rVal }, direction { DirectionSph(thetaVal, phiVal) }
        r { rVal }, theta { thetaVal }, phi { phiVal }
    {}
};

std::ostream& operator<<(std::ostream& os, const DirectionSph dir);
std::ostream& operator<<(std::ostream& os, const PositionSph pos);

// Function to compute the angles for a point in spherical coordinates
DirectionSph computeAnglesFromxyz( float x, float y, float z );

// Function to get the color value of an image in a given angle
PixelRGB mapDirectionToImage( const ImageRGB& image, const DirectionSph& dir );

//=====================================================
//  Propagate a ray backwards from the camera, in a given direction
//=====================================================
DirectionSph traceRayBack(PositionSph cameraPos, DirectionSph viewDir, DirectionSph rayDir);
/* Input:
   - Position of the camera (only r).
   - View direction of the camera, relative to the angular direction.
   - Direction of the incoming ray relative to the view direction of the camera.
*/

//=====================================================
//  Propagate all the rays backwards in a given frame
//=====================================================
void computeFrame(PositionSph cameraPos, DirectionSph viewDir, int width, int height, float thetaFov, ImageRGB background);

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
        - p_phi
*/

// Equation for dr/dlambda
inline float dr_dlambda( float t, float r, float p_r, float theta, float p_theta, float phi, float p_phi )
{
    // return ( 1. - 2. / r ) * p_r;
    return - ( 1. - 2. / r ) * p_r;
}
// Equation for dp_r/dlambda
inline float dp_r_dlambda( float t, float r, float p_r, float theta, float p_theta, float phi, float p_phi )
{
    float sintheta = std::sin(theta);
    // return ( p_theta * p_theta + p_phi * p_phi / (sintheta * sintheta) ) / (r * r * r);
    return - ( p_theta * p_theta + p_phi * p_phi / (sintheta * sintheta) ) / (r * r * r);
}
// Equation for dtheta/dlambda
inline float dtheta_dlambda( float t, float r, float p_r, float theta, float p_theta, float phi, float p_phi )
{
    // return p_theta / (r * r);
    return - p_theta / (r * r);
}
// Equation for dp_theta/dlambda
inline float dp_theta_dlambda( float t, float r, float p_r, float theta, float p_theta, float phi, float p_phi )
{
    float sintheta = std::sin(theta);
    // return std::cos(theta) / (r * r * sintheta * sintheta * sintheta) * p_phi * p_phi;
    return - std::cos(theta) / (r * r * sintheta * sintheta * sintheta) * p_phi * p_phi;
}
// Equation for dphi/dlambda
inline float dphi_dlambda( float t, float r, float p_r, float theta, float p_theta, float phi, float p_phi )
{
    float sintheta = std::sin(theta);
    // return p_phi / ( r * r * sintheta * sintheta );
    return - p_phi / ( r * r * sintheta * sintheta );
}
// Equation for dp_phi/dlambda
// This is constant!
inline float dp_phi_dlambda( float t, float r, float p_r, float theta, float p_theta, float phi, float p_phi )
{
    return 0.;
}

#endif
