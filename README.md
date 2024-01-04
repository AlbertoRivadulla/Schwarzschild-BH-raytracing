# Schwarzschild black hole raytracer

A simple program written in C++ that computes the appearance of a Schwarzschild
black hole and its gravitational lens as seen by an external observer.
It computes the trajectory of the ray of light that reaches each pixel by
integrating backwards the geodesic equations of motion.

The background is given by an image whose pixel locations map linearly to the spherical coordinates $(\theta, \phi)$.
The image used should have the format `.ppm`, and its path must be given as an
argument to the function `loadFromPpm()` in `main.cpp`.
The output of the program is an image with the format `.ppm`.

## Example of the resulting image

![Example result](examples/example.png?raw=true)

