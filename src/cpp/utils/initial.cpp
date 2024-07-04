#include "../auxiliary/list"

#include <iostream>
#include <cmath>
#include <cassert>

Matrix initializePositionsCircle(int numParticles, double majorAxis, double minorAxis, double radius)
{
    Matrix posArray(numParticles, 3);
    double x_min, x_max;
    double y_min, y_max;

    double current_bottom_edge_pos, current_left_edge_pos;
    int xgrid_size, ygrid_size;
    int particle_index;

    y_min = -sqrt(radius * radius - majorAxis * majorAxis);
    y_max = sqrt(radius * radius - majorAxis * majorAxis);

    ygrid_size = static_cast<int>((y_max - y_min) / (2 * minorAxis));

    Vector yGrid(ygrid_size);

    // figure out the y-positions by wrapping each ellipse in bounding box
    current_bottom_edge_pos = y_min;
    for (int i = 0; i < ygrid_size; ++i)
    {
        std::cout << current_bottom_edge_pos << "\n";
        current_bottom_edge_pos += 2 * minorAxis;
        yGrid[i] = current_bottom_edge_pos;
    }

    // now allocate positions to each particle
    current_left_edge_pos = x_min;

    particle_index = 0;
    for (int i = 0; i < ygrid_size; ++i)
    {
        // figure out the range of x positions based on current y position
        x_min = -sqrt(radius * radius - yGrid[i] * yGrid[i]);
        x_max = sqrt(radius * radius - yGrid[i] * yGrid[i]);

        xgrid_size = static_cast<int>((x_max - x_min) / (2 * majorAxis));
        for (int j = 0; j < xgrid_size; ++j)
        {
            current_left_edge_pos += 2 * majorAxis;
        }
    
    } 

    return posArray;
}

int main()
{
    initializePositionsCircle(5, 3, 1, 25);
    return 0;
}