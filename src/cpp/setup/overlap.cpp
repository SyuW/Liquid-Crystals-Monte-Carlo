#include "../auxiliary/list"

// standard library imports
#include <cmath>
#include <iostream>
#include <string>


// check overlap between two identical ellipses -- alternate method using Vieillard-Baron
const bool checkEllipseEllipseOverlap(const double x1, const double y1, const double x2, const double y2, 
                                      const double theta1, const double theta2,
                                      const double minorAxis, const double majorAxis, const bool debug=false)
{
    /*
     * function for checking whether two ellipses overlap, using the Viellard-Baron criterion
     * Arguments       : double x1
     *                   double y1
     *                   double x2
     *                   double y2
     *                   double theta1
     *                   double theta2
     *                   double minorAxis
     *                   double majorAxis
     *                   bool debug
     * return Type     : bool
     */

    double aspectRatio;
    double contactFunction;
    double f_1, f_2;
    double r_dot_u1a, r_dot_u2a, r_dot_u1b, r_dot_u2b;
    double G;
    bool overlap;
    Vector u1a(2), u1b(2), u2a(2), u2b(2), rVec(2);

    // unit vector pointed in direction of major axis of ellipse 1
    u1a[0] = cos(theta1);
    u1a[1] = sin(theta1);

    // unit vector pointed in direction of minor axis of ellipse 1 (modulo sign)
    u1b[0] = -sin(theta1);
    u1b[1] = cos(theta1);

    // unit vector pointed in direction of major axis of ellipse 2
    u2a[0] = cos(theta2);
    u2a[1] = sin(theta2);

    // unit vector pointed in direction of minor axis of ellipse 2 (modulo sign)
    u2b[0] = -sin(theta2);
    u2b[1] = cos(theta2);

    // unit vector pointing from ellipse 1 to ellipse 2
    rVec[0] = x2 - x1;
    rVec[1] = y2 - y1;

    // dot products
    r_dot_u1a = rVec.dot(u1a);
    r_dot_u1b = rVec.dot(u1b);
    r_dot_u2a = rVec.dot(u2a);
    r_dot_u2b = rVec.dot(u2b);

    // aspect ratio is always long axis over short axis
    aspectRatio = majorAxis / minorAxis;

    // alternate method using Vieillard-Baron
    G = 2 + pow((aspectRatio - 1/aspectRatio) * sin(theta1 - theta2), 2);
    f_1 = 1 + G - pow(r_dot_u1a / majorAxis, 2) - pow(r_dot_u1b / minorAxis, 2);
    f_2 = 1 + G - pow(r_dot_u2a / majorAxis, 2) - pow(r_dot_u2b / minorAxis, 2);
    contactFunction = 4 * (f_1 * f_1 - 3 * f_2) * (f_2 * f_2 - 3 * f_1) - pow(9 - f_1 * f_2, 2);
    
    if (contactFunction < 1e-8) // very small contact function indicates ellipses are tangent
        overlap = true;
    else
        overlap = !((contactFunction > 0) && ((f_1 < 0) || (f_2 < 0)));

    // debug information
    if (debug)
    {
        std::string_view msg {};
        if (overlap)
        {
            msg = "The two ellipses overlap";
        }
        else
        {
            msg = "The two ellipses do not overlap";
        }

        std::cout << "-------------------" << "\n";
        std::cout << "Debug option set. Printing out variable info..\n";
        std::cout << "\n";
        std::cout << "Inputs: \n";
        std::cout << "x1: " << x1 << "\n";
        std::cout << "y1: " << y1 << "\n";
        std::cout << "x2: " << x2 << "\n";
        std::cout << "y2: " << y2 << "\n";
        std::cout << "theta1: " << theta1 << "\n";
        std::cout << "theta2: " << theta2 << "\n";
        std::cout << "minorAxis: " << minorAxis << "\n";
        std::cout << "majorAxis: " << majorAxis << "\n";
        std::cout << "\n";
        std::cout << "Outputs: " << "\n";
        std::cout << "Contact function: " << contactFunction << "\n";
        std::cout << "f_1: " << f_1 << "\n";
        std::cout << "f_2: " << f_2 << "\n";
        std::cout << "G: " << G << "\n";
        std::cout << "\n";
        std::cout << msg << "\n";
        std::cout << "-------------------" << "\n";
    }
    
    return overlap;
}

const bool checkBoundaryOverlapCircle(const double R, const double minorAxis, const double majorAxis,
                                      const double xc, const double yc, const double theta, const bool debug=false)
{
    /*
     * function for checking ellipse overlap with a circle - for boundary conditions of container
     * Arguments       : double R
     *                   double minorAxis
     *                   double majorAxis
     *                   double xc
     *                   double yc
     *                   double theta
     *                   bool debug
     * return Type     : bool
     */

    bool overlap {false};
    double cosine {cos(theta)};
    double sine {sin(theta)};
    double sine2 {sin(2 * theta)};
    
    double A, B, C, D, E, P, R_d, D_d;
    double delta, delta_0;

    A = pow(majorAxis, 2) + pow(xc, 2) + pow(yc, 2) - 2 * majorAxis * (xc * cosine + yc * sine) - pow(R, 2);
    B = 4 * minorAxis * (yc * cosine - xc * sine);
    C = 4 * pow(minorAxis, 2) - 2 * pow(majorAxis, 2) + 2 * pow(xc, 2) + 2 * pow(yc, 2) - 2 * pow(R, 2);
    D = 4 * minorAxis * (yc * cosine - xc * sine);
    E = pow(majorAxis, 2) + pow(xc, 2) + pow(yc, 2) + 2 * majorAxis * (xc * cosine + yc * sine) - pow(R, 2);

    P = 8 * A * C - 3 * pow(B, 2);
    R_d = pow(B, 3) + 8 * D * pow(A, 2) - 4 * A * B * C;

    delta = 256 * pow(A * E, 3)
            - 192 * (B * D) * pow(A * E, 2)
            - 128 * pow(A * C * E, 2)
            + 144 * (C * E) * pow(A * D, 2)
            - 27 * pow(A, 2) * pow(D, 4)
            + 144 * (A * C) * pow(B * E, 2)
            - 6 * (A * E) * pow(B * D, 2)
            - 80 * (A * B * D * E) * pow(C, 2)
            + 18 * (A * B * C) * pow(D, 3)
            + 16 * (A * E) * pow(C, 4)
            - 4 * A * pow(C, 3) * pow(D, 2)
            - 27 * pow(E, 2) * pow(B, 4)
            + 18 * (C * D * E) * pow(B, 3)
            - 4 * pow(B * D, 3)
            - 4 * E * pow(B, 2) * pow(C, 3)
            + pow(B, 2) * pow(C, 2) * pow(D, 2);

    delta_0 = pow(C, 2) - 3 * B * D + 12 * A * E;

    D_d = 64 * pow(A, 3) * E - 16 * pow(A, 2) * pow(C, 2) + 16 * A * pow(B, 2) * C - 16 * pow(A, 2) * B * D - 3 * pow(B, 4);

    if (delta < 0)
        overlap = true;
    else if (delta > 0 && P < 0 && D_d < 0)
        overlap = true;
    else if (delta == 0)
    {
        if (P < 0 && D_d < 0 && delta_0 != 0)
            overlap = true;
        else if (D_d > 0 || (P > 0 && (D_d != 0 || R_d != 0)))
            overlap = false;
        else if (delta_0 == 0 && D_d != 0)
            overlap = true;
        else if (D_d == 0 && P < 0)
            overlap = true;
        else if (D_d == 0 && P > 0 && R_d == 0)
            overlap = false;
        else
            overlap = false;
    }

    if (debug)
    {
        std::cout << "-------------------" << "\n";
        std::cout << "Debug option set. Printing out variable info..\n";
        std::cout << "\n";
        std::cout << "Inputs:\n";
        std::cout << "Boundary radius: " << R << "\n";
        std::cout << "Ellipse center x-position: " << xc << "\n";
        std::cout << "Ellipse center y-position: " << yc << "\n";
        std::cout << "Ellipse rotation angle: " << theta << "\n";
        std::cout << "Major axis: " << majorAxis << "\n";
        std::cout << "Minor axis: " << minorAxis << "\n";
        std::cout << "\n";
        std::cout << "Outputs:" << "\n";
        std::cout << "A: " << A << "\n";
        std::cout << "B: " << B << "\n";
        std::cout << "C: " << C << "\n";
        std::cout << "D: " << D << "\n";
        std::cout << "E: " << E << "\n";
        std::cout << "P: " << P << "\n";
        std::cout << "R_d: " << R_d << "\n";
        std::cout << "D_d: " << D_d << "\n";
        std::cout << "delta: " << delta << "\n";
        std::cout << "delta_0: " << delta_0 << "\n";
        std::cout << "overlap: " << overlap << "\n";
        std::cout << "-------------------" << "\n";
    }

    return overlap;
}

const bool checkBoundaryOverlapLine(const double slope, const double intercept, const double minorAxis, const double majorAxis,
                                    const double xc, const double yc, const double theta, const bool debug=false)
{
    bool overlap { false };
    double cosine { cos(theta) };
    double sine { sin(theta) };

    double A, B, C;
    double delta;

    A = majorAxis * (slope * cosine - sine) + yc - slope * xc - intercept;
    B = 2 * minorAxis * (cosine + slope * sine);
    C = majorAxis * (sine - slope * cosine) + yc - slope * xc - intercept;

    delta = pow(B, 2) - 4 * A * C;

    if (delta > 0)
        overlap = true; 
    else
        overlap = false;

    if (debug)
    {
        std::cout << "------------------------" << "\n";
        std::cout << "Debug option set. Printing out variable info..\n";
        std::cout << "\n";
        std::cout << "Inputs:\n";
        std::cout << "Slope of line: " << slope << "\n";
        std::cout << "Intercept of line: " << intercept << "\n";
        std::cout << "Ellipse center x-position: " << xc << "\n";
        std::cout << "Ellipse center y-position: " << yc << "\n";
        std::cout << "Ellipse rotation angle: " << theta << "\n";
        std::cout << "Major axis: " << majorAxis << "\n";
        std::cout << "Minor axis: " << minorAxis << "\n";
        std::cout << "\n";
        std::cout << "Outputs:" << "\n";
        std::cout << "A: " << A << "\n";
        std::cout << "B: " << B << "\n";
        std::cout << "C: " << C << "\n";
        std::cout << "delta: " << delta << "\n";
        std::cout << "overlap: " << overlap << "\n";
        std::cout << "-------------------" << "\n";
    }

    return overlap;
}

// degenerate case for computing ellipse-line overlap when line is vertical
bool checkBoundaryOverlapVertical(const double xIntercept, const double minorAxis, const double majorAxis,
                                  const double xc, const double yc, const double theta, const bool debug=false)
{
    bool overlap { false };
    double cosine { cos(theta) };
    double sine { sin(theta) };
    double A, B, C;
    double delta;

    A = xc - xIntercept - majorAxis * cosine;
    B = -2 * minorAxis * sine;
    C = majorAxis * cosine + xc - xIntercept;

    delta = pow(B, 2) - 4 * A * C;

    if (delta > 0)
        overlap = true;
    else
        overlap = false;

    if (debug)
    {
        std::cout << "-------------------------" << "\n";
        std::cout << "Debug option set. Printing out variable info...\n";
        std::cout << "\n";
        std::cout << "Inputs:\n";
        std::cout << "x-intercept of vertical line: " << xIntercept << "\n";
        std::cout << "Ellipse center x-position: " << xc << "\n";
        std::cout << "Ellipse center y-position: " << yc << "\n";
        std::cout << "Ellipse rotation angle: " << theta << "\n";
        std::cout << "Major axis: " << majorAxis << "\n";
        std::cout << "Minor axis: " << minorAxis << "\n";
        std::cout << "\n";
        std::cout << "Outputs:" << "\n";
        std::cout << "A: " << A << "\n";
        std::cout << "B: " << B << "\n";
        std::cout << "C: " << C << "\n";
        std::cout << "delta: " << delta << "\n";
        std::cout << "overlap: " << overlap << "\n";
        std::cout << "-------------------" << "\n";
    }

    return overlap;
}

// check overlap between two identical ellipses
bool checkEllipseEllipseOverlap_old(double x1, double y1, double x2, double y2, double theta1, double theta2, double minorAxis, double majorAxis)
{
    double aspectRatio, xi, w, centerDistance, sigma2D;
    double r_dot_u1, r_dot_u2, u1_dot_u2;
    double denom;
    bool val;
    Vector u1(2), u2(2), rVec(2);

    // unit vector pointed in direction of major axis of ellipse 1
    u1[0] = cos(theta1);
    u1[1] = sin(theta1);
    // unit vector pointed in direction of major axis of ellipse 2
    u2[0] = cos(theta2);
    u2[1] = sin(theta2);

    // unit vector pointing from ellipse 1 to ellipse 2
    centerDistance = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)); // center to center distance between ellipses
    rVec[0] = (x2 - x1) / centerDistance;
    rVec[1] = (y2 - y1) / centerDistance;

    // dot products
    r_dot_u1 = rVec.dot(u1);
    r_dot_u2 = rVec.dot(u2);
    u1_dot_u2 = u1.dot(u2);

    aspectRatio = majorAxis / minorAxis; // aspect ratio is always long axis over short axis
    xi = (pow(aspectRatio, 2) - 1) / (pow(aspectRatio, 2) + 1);

    denom = pow(
                1 - 0.5 * xi * 
                (
                    (pow(r_dot_u1 + r_dot_u2, 2) / (1 + xi * u1_dot_u2))
                    + 
                    (pow(r_dot_u1 - r_dot_u2, 2) / (1 - xi * u1_dot_u2))
                ), 
                0.5);
    
    w = 2 * minorAxis;
    sigma2D = w / denom;

    if (centerDistance >= sigma2D) {
        val = false;
    } else {
        val = true;
    }

    // debug information

    return val;
}