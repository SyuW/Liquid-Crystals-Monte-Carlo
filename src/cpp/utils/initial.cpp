#include "../auxiliary/list"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>


Matrix initializePositionsCircle(int numParticles, double majorAxis, double minorAxis, double radius)
{
    Matrix posArray(numParticles, 3);
    double x_min, x_max;
    double y_min, y_max;
    double current_bottom_edge_pos, current_left_edge_pos;
    int particle_index, numParticlesPossible;
    int xgrid_size, ygrid_size;
    bool arrayComplete {false};

    y_min = -sqrt(radius * radius - majorAxis * majorAxis);
    y_max = sqrt(radius * radius - majorAxis * majorAxis);

    ygrid_size = static_cast<int>((y_max - y_min) / (2 * minorAxis));

    Vector yGrid(ygrid_size);

    // figure out the y-positions by wrapping each ellipse in bounding box
    current_bottom_edge_pos = y_min;
    for (int i = 0; i < ygrid_size; ++i)
    {
        yGrid[i] = current_bottom_edge_pos;
        current_bottom_edge_pos += 2 * minorAxis;
    }

    // now find the x-positions of each particle
    particle_index = 0;
    for (int i = 0; i < ygrid_size; ++i)
    {
        x_min = -sqrt(radius * radius - yGrid[i] * yGrid[i]);
        x_max = sqrt(radius * radius - yGrid[i] * yGrid[i]);
        current_left_edge_pos = x_min;
        xgrid_size = static_cast<int>((x_max - x_min) / (2 * majorAxis));

        // figure out the range of x positions based on current y position
        for (int j = 0; j < xgrid_size; ++j)
        {
            posArray(particle_index, 0) = current_left_edge_pos;
            posArray(particle_index, 1) = yGrid[i];
            posArray(particle_index, 2) = 0.0; 

            current_left_edge_pos += 2 * majorAxis;
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
    }

    // transform to ellipse centers
    for (int i = 0; i < numParticles; ++i)
    {
        posArray(i, 0) += majorAxis;
        posArray(i, 1) += minorAxis; 
    }

    posArray.print();

    return posArray;
}

int main()
{
    int numParticlesToSimulate {7};
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

    for (int particleIndex=0; particleIndex < numParticlesToSimulate; ++particleIndex)
    {
        outFile << initialPositions(particleIndex, 0) 
                << " " 
                << initialPositions(particleIndex, 1) 
                << " " 
                << initialPositions(particleIndex, 2) << "\n";
    }

    return 0;
}