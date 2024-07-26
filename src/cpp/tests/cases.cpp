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
class ec_TestCase
{
    int a {5};
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


bool runEllipseEllipseTest(ee_TestCase test)
{
    bool result;
    result = checkEllipseEllipseOverlap(test.x1, test.y1, test.x2, test.y2,
                                        test.t1, test.t2,
                                        test.minorAxis, test.majorAxis);
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

    // ellipse-ellipse overlap tests

    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(-1, 0, 1, 0, PI/2, PI/2, 1, 3, true));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(-1-0.01, 0, 1, 0, PI/2, PI/2, 1, 3, false));
    ellipseEllipseTests.push_back(makeEllipseEllipseTestCase(0, 0, 8, 0, 0, PI/2, 1, 7, true));

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

}