#ifndef OVERLAPHEADERDEF
#define OVERLAPHEADERDEF

const bool checkEllipseEllipseOverlap(const double x1, const double y1, const double x2, const double y2, 
                                      const double theta1, const double theta2,
                                      const double minorAxis, const double majorAxis, const bool debug=false);

const bool checkBoundaryOverlapCircle(const double R, const double minorAxis, const double majorAxis,
                                      const double xc, const double yc, const double theta, const bool debug=false);

const bool checkBoundaryOverlapLine(const double slope, const double intercept, const double minorAxis, const double majorAxis,
                                    const double xc, const double yc, const double theta, const bool debug=false);

bool checkBoundaryOverlapVertical(const double xIntercept, const double minorAxis, const double majorAxis,
                                  const double xc, const double yc, const double theta, const bool debug=false);

#endif