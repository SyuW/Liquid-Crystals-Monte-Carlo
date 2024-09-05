#ifndef CELLLISTHEADERDEF
#define CELLLISTHEADERDEF

#include <cmath>

class CellList
{
public:
    // constructors
    CellList(double cellWidth, double height, double width);            // box boundary condition
    CellList(double cellWidth, double boundaryRadius);                  // circle boundary condition
    CellList(double cellWidth, double innerRadius, double outerRadius); // annulus boundary condition

private:

};

#endif

