#ifndef OVERLAPHEADERDEF
#define OVERLAPHEADERDEF

// Helper for checking if two identical ellipses overlap, using Vieillard-Baron (1970) criterion.
const bool checkEllipseEllipseOverlap(const double x1, const double y1, const double x2, const double y2,
                                      const double theta1, const double theta2,
                                      const double minorAxis, const double majorAxis, const bool debug = false);

// Helper for checking ellipse overlap with a circle - for boundary conditions of circle container.
const bool checkBoundaryOverlapCircle(const double R, const double minorAxis, const double majorAxis,
                                      const double xc, const double yc, const double theta, const bool debug = false);

const bool checkBoundaryOverlapLine(const double slope, const double intercept, const double minorAxis, const double majorAxis,
                                    const double xc, const double yc, const double theta, const bool debug = false);

// For handling degenerate case of computing ellipse-line overlap when the line is vertical.
bool checkBoundaryOverlapVertical(const double xIntercept, const double minorAxis, const double majorAxis,
                                  const double xc, const double yc, const double theta, const bool debug = false);

#endif