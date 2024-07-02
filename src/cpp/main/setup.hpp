class Setup
{
public:
    Setup(void);
    ~Setup(void);

    // function prototypes
    bool checkEllipseEllipseOverlap(double x1, double y1, double x2, double y2, double theta1, double theta2, double minorAxis, double majorAxis, bool verbose);
    bool checkBoundaryOverlapCircle(double R, double minorAxis, double majorAxis, double xc, double yc, double theta);

private:

};

#include "../lib/overlap.cpp"