#ifndef IO_FILES_H
#define IO_FILES_H

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>

// Struct for pixels with colors
struct PixelRGB
{
    int R;
    int G;
    int B;

    PixelRGB() : 
        R { 0 } , G { 0 }, B { 0 } 
    {}

    PixelRGB( int r, int g, int b ) :
        R { r }, G { g }, B { b }
    {}
};

std::ostream& operator<<(std::ostream& os, const PixelRGB pixel);
std::ifstream& operator>>(std::ifstream &in, PixelRGB& pixel);

// Struct for all the pixels in an image, with its size
struct ImageRGB
{
    // Vector of pixels
    std::vector<PixelRGB> pixels;
    // Dimensions
    int width;
    int height;
    // Constructor
    ImageRGB() : 
        pixels { std::vector<PixelRGB> ( 1, PixelRGB() ) }, width { 1 }, height { 1 }
    {}
    ImageRGB( int w, int h ) :
        pixels { std::vector<PixelRGB> ( w * h, PixelRGB() ) } , width { w }, height { h }
    {}
    ImageRGB( std::vector<PixelRGB> pixelValues, int w, int h) :
        pixels { pixelValues }, width { w }, height { h }
    {}
};

// Save a the color data in a frame buffer to a ppm file. The data must be given as
// a vector of integers, from 0 to 255, with the given amount of chanels per pixel.
void saveToPpm(std::string fileName, const std::vector<int>& frameBuffer, int width, int height, int channels = 3);
void saveToPpm(std::string fileName, const ImageRGB& image);

void outputDataToFile(std::vector<std::vector<float>> data, std::string fileName);

// Function to load the data from a ppm file
ImageRGB loadFromPpm(std::string fileName);

#endif
