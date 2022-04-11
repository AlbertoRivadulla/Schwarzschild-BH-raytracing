#include "ioFiles.h"

// Print method for the PixelRGB struct
std::ostream& operator<<(std::ostream& os, const PixelRGB pixel)
{
    os << "RGB = (" << pixel.R << ", " << pixel.G << ", " << pixel.B << ")";
    return os;
}
// Method for assigning values to a PixelRGB struct with the >> operator
std::ifstream& operator>>(std::ifstream &in, PixelRGB& pixel)
{
    in >> pixel.R;
    in >> pixel.G;
    in >> pixel.B;

    // in >> pixel.R >> pixel.G >> pixel.B;
    return in;
}

// Save a the color data in a frame buffer to a ppm file. The data must be given as
// a vector of integers, from 0 to 255, with the given amount of chanels per pixel.
void saveToPpm(std::string fileName, const std::vector<int>& frameBuffer, int width, int height, int channels)
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

// Save the data in a ImageRGB instance to a file
void saveToPpm(std::string fileName, const ImageRGB& image)
{
    // Convert the data to a vector
    std::vector<int> frameBuffer ( image.width * image.height * 3 , 0 );
    int index = 0;
    for (auto px : image.pixels)
    {
        frameBuffer[ index++ ] = px.R;
        frameBuffer[ index++ ] = px.G;
        frameBuffer[ index++ ] = px.B;
    }

    // Call the function inplemented for a vector of image data
    saveToPpm( fileName, frameBuffer, image.width, image.height, 3 );
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

    // Close the file
    fileStream.close();
}

// Function to load the data from a ppm file
ImageRGB loadFromPpm(std::string fileName)
{
    // Open the file
    // std::ifstream file(fileName.c_str());
    std::ifstream file;
    file.open(fileName.c_str(), std::ifstream::in);
    if (!file.is_open())
    {
        std::cout << "Failed to open the file " << fileName << '\n';
    }

    // Set this to not skip spaces, which also carry a value
    file >> std::noskipws;

    // Read the preamble
    std::string line;
    // The first line is discarded
    std::getline(file, line);
    // Get the width and height of the image from the second line
    int width, height;
    std::getline(file, line);
    std::stringstream ss;       // Convert the line to a stringstream
    ss << line;
    ss >> width >> height;      // Get the width and height from this
    // The third line is discarded
    std::getline(file, line);

    // Create a vector with the pixels
    std::vector<PixelRGB> pixels ( width * height );
    // Define some auxiliary variables for reading the pixels
    char c; // Variable to read each value, which ranges from 0 to 255 (4 bits, the size of a char)
    int red = -1, green = -1, blue = -1; // Values of the colors in the current pixel
    PixelRGB pixel ( -1, -1, -1 );
    int rgbCount = 0; // Counter, that will become 3 when one entire pixel has been read
    // Iterate over the different pixels
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            // Read the red value
            file >> c;
            red = (int)(unsigned char) c;
            // Read the green value
            file >> c;
            green = (int)(unsigned char) c;
            // Read the blue value
            file >> c;
            blue = (int)(unsigned char) c;

            // Store the values of the colors in the current pixel in the image
            pixels[ width * j + i ] = PixelRGB( red, green, blue );
        }
    }

    // Close the file
    file.close();

    // Construct the image struct and return it
    return ImageRGB(pixels, width, height);
}
