#include "saveToFiles.h"

// Save a the color data in a frame buffer to a ppm file. The data must be given as
// a vector of integers, from 0 to 255, with the given amount of chanels per pixel.
void saveToPpm(std::string fileName, std::vector<int> frameBuffer, int width, int height, int channels)
{
    // Open a file
    std::ofstream fileStream;
    fileStream.open(fileName);

    // Write the file in binary format (taking 1.2 MB in this case)
    // Write the header of the file, as described in the Wikipedia article
    fileStream << "P6\n";                           // Type of the file (here, RGB color image in ASCII)
    fileStream << width << ' ' << height << '\n';   // Width and height in pixels
    fileStream << "255\n";                          // Maximum value for each color

    // Iterate over the different pixels
    for (int px = 0; px < width * height; ++px)
    {
        // Iterate over the three channels in each pixel
        for (int ch = 0; ch < channels; ++ch)
        {
            // Clip the values between 0 and 255
            fileStream << static_cast<char>(std::max(0, std::min(255, frameBuffer[3 * px + ch])));
        }
    }

    // Close the file
    fileStream.close();
}

void outputDataToFile(std::vector<std::vector<float>> data, std::string fileName)
{
    // Open a file to output the results
    std::ofstream fileStream;
    fileStream.open(fileName);

    // Save the values in each point to the file
    for (int i = 0; i < data[0].size(); ++i)
    {
        for (int j = 0; j < data.size() - 1; ++j)
        {
            fileStream << std::to_string(data[j][i]) << ';';
        }
        fileStream << std::to_string(data[data.size() - 1][i]) << '\n';
    }

    // close the file
    fileStream.close();
}
