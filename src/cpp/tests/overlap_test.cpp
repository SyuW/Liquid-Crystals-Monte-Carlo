#include "../auxiliary/list"
#include "../setup/overlap.cpp"
#include "../main/constants.hpp"

#include <vector>


// ellipse-ellipse overlap test case
struct ee_TestCase
{
    double x1;
    double y1;
    double x2;
    double y2;
    double t1;
    double t2;
    double minorAxis;
    double majorAxis;
    bool answer;
};


// ellipse-circle overlap test case
struct ec_TestCase
{
    double R;
    double minorAxis;
    double majorAxis;
    double xc;
    double yc;
    double theta;
    bool answer;
};


ee_TestCase makeEllipseEllipseTestCase(double x1, double y1, double x2, double y2, 
                                       double t1, double t2, 
                                       double minorAxis, double majorAxis, bool answer)
{
    ee_TestCase test;

    test.x1 = x1;
    test.y1 = y1;
    test.x2 = x2;
    test.y2 = y2;
    test.t1 = t1;
    test.t2 = t2;
    test.minorAxis = minorAxis;
    test.majorAxis = majorAxis;
    test.answer = answer;

    return test; 
}


ec_TestCase makeEllipseCircleTestCase(double R, double minorAxis, double majorAxis,
                                      double xc, double yc, double theta,
                                      bool answer)
{
    ec_TestCase test;

    test.R = R;
    test.minorAxis = minorAxis;
    test.majorAxis = majorAxis;
    test.xc = xc;
    test.yc = yc;
    test.theta = theta;
    test.answer = answer;

    return test;
}


void printEETestCase(ee_TestCase test)
{
    std::cout << "EETestCase(";
    std::cout << "params={";
    std::cout << test.x1 << ", ";
    std::cout << test.y1 << ", ";
    std::cout << test.x2 << ", ";
    std::cout << test.y2 << ", ";
    std::cout << test.t1 << ", ";
    std::cout << test.t2 << ", ";
    std::cout << test.minorAxis << ", ";    
    std::cout << test.majorAxis << "}, ";
    std::cout << "Answer=" << test.answer << ")";        
}


void printECTestCase(ec_TestCase test)
{
    std::cout << "ECTestCase(";
    std::cout << "params={";
    std::cout << test.R << ", ";
    std::cout << test.minorAxis << ", ";
    std::cout << test.majorAxis << ", ";
    std::cout << test.xc << ", ";
    std::cout << test.yc << ", ";
    std::cout << test.theta << ", ";
    std::cout << "Answer=" << test.answer << ")";
}


// feeling cute might template later
bool runEllipseEllipseTest(ee_TestCase test)
{
    bool result;
    result = checkEllipseEllipseOverlap(test.x1, test.y1, test.x2, test.y2,
                                        test.t1, test.t2,
                                        test.minorAxis, test.majorAxis, false);
    if (result == test.answer)
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool runEllipseCircleTest(ec_TestCase test)
{
    bool result;
    result = checkBoundaryOverlapCircle(test.R, test.minorAxis, test.majorAxis,
                                        test.xc, test.yc, test.theta, true);
    if (result == test.answer)
    {
        return true;
    }
    else
    {
        return false;
    }
}


int main()
{
    std::vector<ee_TestCase> ellipseEllipseTests;
    std::vector<ec_TestCase> ellipseCircleTests;

    // ellipse-ellipse overlap tests
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(-1, 0, 1, 0, PI/2, PI/2, 1, 3, true));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(-1-0.01, 0, 1, 0, PI/2, PI/2, 1, 3, false));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(0, 0, 8, 0, 0, PI/2, 1, 7, true));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(0, 0, 0, 0, 0, PI/2, 1, 1, true));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(0, 0, 0, 0, 0, PI/2, 7, 1, true));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(-5, 0, -6, 10, -0.9, 0.9, 1, 7, true));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(-5, 0, 1, 10, -0.9, 0.9, 1, 7, false));

    bool testOutcome;
    for (const auto& eetest: ellipseEllipseTests)
    {
        testOutcome = runEllipseEllipseTest(eetest);

        std::cout << "Running: ";
        printEETestCase(eetest);

        if (testOutcome)
            std::cout << ", passed.\n";
        else
            std::cout << ", failed.\n";
    }

    // ellipse-circle overlap tests
    ellipseCircleTests.push_back(makeEllipseCircleTestCase(4, 1, 3, 1, 1, 7, true));
    ellipseCircleTests.push_back(makeEllipseCircleTestCase(6.6, 1, 3, 4, 1, 4, true));
    ellipseCircleTests.push_back(makeEllipseCircleTestCase(3.5, 1, 3, 0, 0, 0, false));
    ellipseCircleTests.push_back(makeEllipseCircleTestCase(3.5, 1, 3, 0.5, 0, 0.8, false));
    ellipseCircleTests.push_back(makeEllipseCircleTestCase(3.5, 1, 3, 0.6, 0, 0, true));
    ellipseCircleTests.push_back(makeEllipseCircleTestCase(25, 1, 5, 22.8785, 7.73919, -4.33774, true));

    for (const auto& ectest: ellipseCircleTests)
    {
        testOutcome = runEllipseCircleTest(ectest);
        std::cout << "Running: ";
        printECTestCase(ectest);

        if (testOutcome)
            std::cout << ", passed.\n";
        else
            std::cout << ", failed.\n";
    }

}