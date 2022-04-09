#ifndef SAVE_TO_FILES_H
#define SAVE_TO_FILES_H

#include <fstream>
#include <vector>
#include <iostream>

// Save a the color data in a frame buffer to a ppm file. The data must be given as
// a vector of integers, from 0 to 255, with the given amount of chanels per pixel.
void saveToPpm(std::string fileName, std::vector<int> frameBuffer, int width, int height, int channels = 3);

void outputDataToFile(std::vector<std::vector<float>> data, std::string fileName);

#endif
