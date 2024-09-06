#ifndef OVERLAPTESTHEADER
#define OVERLAPTESTHEADER

struct ee_TestCase;
struct ec_TestCase;

ee_TestCase makeEllipseEllipseTestCase(double x1, double y1, double x2, double y2, 
                                       double t1, double t2, 
                                       double minorAxis, double majorAxis, bool answer);

ec_TestCase makeEllipseCircleTestCase(double R, double minorAxis, double majorAxis,
                                      double xc, double yc, double theta,
                                      bool answer);

void printEETestCase(ee_TestCase test);
void printECTestCase(ec_TestCase test);
bool runEllipseEllipseTest(ee_TestCase test);
bool runEllipseCircleTest(ec_TestCase test);

#endif