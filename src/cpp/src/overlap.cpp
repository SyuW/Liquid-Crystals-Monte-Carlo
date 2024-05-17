#include "types.h"
#include <cmath>

bool checkOverlap(REAL x1, REAL y1, REAL x2, REAL y2, REAL theta1, REAL theta2, REAL minorAxis, REAL majorAxis)
{
    REAL aspectRatio, xi, w, rVal, sigma2D;
    bool val;

    /* default initialize the array types here*/

    aspectRatio = majorAxis / minorAxis;
    // using square root function from the standard library
    xi = (pow(aspectRatio, 2) - 1) / (pow(aspectRatio, 2) + 1);
    w = 2 * minorAxis;

    // set rVal and sigma to arbitrary values for now
    rVal = 2;
    sigma2D = 3;

    if (rVal >= sigma2D) {
        val = true;
    } else {
        val = false;
    }

    return val;
}