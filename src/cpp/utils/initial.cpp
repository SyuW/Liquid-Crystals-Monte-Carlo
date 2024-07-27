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
    int particle_index;
    int xgrid_size, ygrid_size;
    bool arrayComplete {false};

    y_min = -sqrt(radius * radius - majorAxis * majorAxis);
    y_max =  sqrt(radius * radius - majorAxis * majorAxis);

    ygrid_size = static_cast<int>((y_max - y_min) / (2 * minorAxis));

    // figure out the y-positions by wrapping each ellipse in bounding box
    bottom_edge_pos = y_min;
    particle_index = 0;
    for (int i = 0; i < ygrid_size; ++i)
    {
        // find the x-positions for each particle
        x_min = -sqrt(radius * radius - bottom_edge_pos * bottom_edge_pos);
        x_max =  sqrt(radius * radius - bottom_edge_pos * bottom_edge_pos);
        left_edge_pos = x_min;
        xgrid_size = static_cast<int>((x_max - x_min) / (2 * majorAxis));

        for (int j = 0; j < xgrid_size; ++j)
        {
            // transform from bottom-left corner of bounding box to ellipse center
            posArray(particle_index, 0) = left_edge_pos + majorAxis;
            posArray(particle_index, 1) = bottom_edge_pos + minorAxis;
            posArray(particle_index, 2) = 0;

            std::cout << "(x,y)=" << "(" << left_edge_pos + majorAxis << ", " << bottom_edge_pos + minorAxis << ")\n";

            left_edge_pos += 2 * majorAxis;
            particle_index += 1;

            if (particle_index >= numParticles)
            {
                arrayComplete = true;
                break;
            }
        }

        if (arrayComplete)
        {
            break;
        }

        bottom_edge_pos += 2 * minorAxis;
    }

    return posArray;
}

int main()
{
    int numParticlesToSimulate {120};
    double majorAxis {3};
    double minorAxis {1};
    double boundaryRadius {25};

    Matrix initialPositions { initializePositionsCircle(numParticlesToSimulate, majorAxis, minorAxis, boundaryRadius) };

    std::ofstream outFile { "initialPositions.txt" };

    if (!outFile)
    {
        // Print an error and exit
        std::cerr << "Could not open initialPositions.txt for writing.\n";
        return 1;
    }

    outFile << 

    for (int particleIndex=0; particleIndex < numParticlesToSimulate; ++particleIndex)
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