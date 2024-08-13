#include "../auxiliary/list"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <vector>


Matrix initializePositionsBox(const int numParticles, const double majorAxis, const double minorAxis,
                              const double height, const double width)
{
    /*
     * numParticles
     * majorAxis
     * minorAxis
     * height
     * width
     */

    Matrix posArray(numParticles, 3);
    double x_min, x_max;
    double y_min, y_max;
    double bottom_edge_pos, left_edge_pos;
    int particleIndex;
    int xgrid_size;
    bool arrayComplete {false};

    y_min = 0;
    y_max = height - 2 * minorAxis;

    // use bounding box strategy
    bottom_edge_pos = y_min;

    particleIndex = 0;
    while (bottom_edge_pos < y_max)
    {
        x_min = 0;
        x_max = width;

        left_edge_pos = x_min;
        xgrid_size = static_cast<int>((x_max - x_min) / (2 * majorAxis));

        for (int j = 0; j < xgrid_size; ++j)
        {
            posArray(particleIndex, 0) = left_edge_pos + majorAxis;
            posArray(particleIndex, 1) = bottom_edge_pos + minorAxis;
            posArray(particleIndex, 2) = 0;
            particleIndex += 1;

            if (particleIndex >= numParticles)
            {
                arrayComplete = true;
                break;
            }
            // else, find position of next particle at same y-coord
            else
            {
                left_edge_pos += 2 * majorAxis;
            }
        }

        if (arrayComplete)
            break;
        else
            bottom_edge_pos += 2 * minorAxis;
    }

    if (!arrayComplete)
    {
        std::cout << "Number of particles in requested initial config exceeds capacity of container, "
                  << "returning resized positions array with " << particleIndex+1 << " particles." << "\n";

        // return resized positions array
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


Matrix initializePositionsCircle(const int numParticles, const double majorAxis, const double minorAxis,
                                 const double radius)
{
    /*
     * numParticles
     * majorAxis
     * minorAxis
     * radius
     *
     */

    Matrix posArray(numParticles, 3);
    double x_min, x_max;
    double y_min, y_max;
    double bottom_edge_pos, left_edge_pos;
    double intersect;
    int particleIndex;
    int xgrid_size;
    bool arrayComplete {false};

    y_min = -sqrt(radius * radius - majorAxis * majorAxis);
    y_max =  sqrt(radius * radius - majorAxis * majorAxis) - 2 * minorAxis;

    // figure out the y-positions by wrapping each ellipse in bounding box
    bottom_edge_pos = y_min;
    particleIndex = 0;
    while (bottom_edge_pos < y_max)
    {
        // find the x-positions for each particle
        if (bottom_edge_pos < 0)
            intersect = bottom_edge_pos;
        else
            intersect = bottom_edge_pos + 2 * minorAxis;

        x_min = -sqrt(radius * radius - intersect * intersect);
        x_max =  sqrt(radius * radius - intersect * intersect);

        left_edge_pos = x_min;
        xgrid_size = static_cast<int>((x_max - x_min) / (2 * majorAxis));

        for (int j = 0; j < xgrid_size; ++j)
        {
            posArray(particleIndex, 0) = left_edge_pos + majorAxis;
            posArray(particleIndex, 1) = bottom_edge_pos + minorAxis;
            posArray(particleIndex, 2) = 0;
            particleIndex += 1;

            if (particleIndex >= numParticles)
            {
                arrayComplete = true;
                break;
            }
            // else, find position of next particle at same y-coord
            else
            {
                left_edge_pos += 2 * majorAxis;
            }
        }

        if (arrayComplete)
            break;
        else
            bottom_edge_pos += 2 * minorAxis;
    }

    if (!arrayComplete)
    {
        std::cout << "Number of particles in requested initial config exceeds capacity of container,"
                  << "returning resized positions array with " << particleIndex << " particles." << "\n";

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