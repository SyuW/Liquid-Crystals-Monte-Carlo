class Setup
{
public:
    Setup(void);
    ~Setup(void);

    // function prototypes
    bool checkOverlap(double x1, double y1, double x2, double y2, double theta1, double theta2, double minorAxis, double majorAxis, bool verbose);
    
private:

};

#include "../lib/overlap.cpp"