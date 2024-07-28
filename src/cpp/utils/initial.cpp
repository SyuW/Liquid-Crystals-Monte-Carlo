#include "../auxiliary/list"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <vector>


Matrix initializePositionsCircle(int numParticles, double majorAxis, double minorAxis, double radius)
{
    Matrix posArray(numParticles, 3);
    double x_min, x_max;
    double y_min, y_max;
    double bottom_edge_pos, left_edge_pos;
    int particleIndex;
    int xgrid_size, ygrid_size;
    bool arrayComplete {false};

    y_min = -sqrt(radius * radius - majorAxis * majorAxis);
    y_max = 0;

    // ygrid_size = static_cast<int>((y_max - y_min) / (2 * minorAxis));

    // figure out the y-positions by wrapping each ellipse in bounding box
    bottom_edge_pos = y_min;
    particleIndex = 0;
    while (bottom_edge_pos < y_max)
    {
        // find the x-positions for each particle
        x_min = -sqrt(radius * radius - bottom_edge_pos * bottom_edge_pos);
        x_max =  sqrt(radius * radius - bottom_edge_pos * bottom_edge_pos);

        left_edge_pos = x_min;
        xgrid_size = static_cast<int>((x_max - x_min) / (2 * majorAxis));

        for (int j = 0; j < xgrid_size; ++j)
        {

            if (particleIndex >= numParticles)
            {
                arrayComplete = true;
                break;
            }
            else
            {
                // transform from bottom-left corner of bounding box to ellipse center
                posArray(particleIndex, 0) = left_edge_pos + majorAxis;
                posArray(particleIndex, 1) = bottom_edge_pos + minorAxis;
                posArray(particleIndex, 2) = 0;
                particleIndex += 1;
            }

            if (particleIndex >= numParticles)
            {
                arrayComplete = true;
                break;
            }
            else
            {
                // by symmetry, can reflect start position about x-axis to get equally valid start
                posArray(particleIndex, 0) = left_edge_pos + majorAxis;
                posArray(particleIndex, 1) = -(bottom_edge_pos + minorAxis);
                posArray(particleIndex, 2) = 0;
                particleIndex += 1;
            }

            left_edge_pos += 2 * majorAxis;
        }

        if (arrayComplete)
        {
            break;
        }
        else
        {
            bottom_edge_pos += 2 * minorAxis;
        }
    }

    if (!arrayComplete)
    {
        std::cout << "Number of particles in requested initial config exceeds capacity of container,"
                  << "returning resized positions array with " << particleIndex+1 << " particles." << "\n";

        // return resized position array
        Matrix resized(particleIndex, 3);
        for (int i=0; i < particleIndex; ++i)
        {
            resized(i, 0) = posArray(i, 0);
            resized(i, 1) = posArray(i, 1);
            resized(i, 2) = posArray(i, 2);
        }

        return resized;
    }

    return posArray;
}

int main()
{
    int numParticlesToSimulate {100};
    double majorAxis {3};
    double minorAxis {3};
    double boundaryRadius {25};

    Matrix initialPositions { initializePositionsCircle(numParticlesToSimulate, majorAxis, minorAxis, boundaryRadius) };

    std::ofstream outFile { "initialPositions.txt" };

    if (!outFile)
    {
        // Print an error and exit
        std::cerr << "Could not open initialPositions.txt for writing.\n";
        return 1;
    }

    outFile << "# x y theta\n";

    for (int particleIndex=0; particleIndex < initialPositions.getNumberOfRows(); ++particleIndex)
    {
        outFile << initialPositions(particleIndex, 0) 
                << "," 
                << initialPositions(particleIndex, 1) 
                << "," 
                << initialPositions(particleIndex, 2) << "\n";
    }

    outFile.close();

    return 0;
}