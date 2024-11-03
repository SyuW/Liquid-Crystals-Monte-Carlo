typedef unsigned long int ULONG;
typedef unsigned short int USHORT;
typedef signed long int LONG;
typedef signed short int SHORT;
typedef double REAL;

struct circleSimData
{
    int numMonteCarloSteps {};
    int boundaryRadius {};
    double majorAxis {};
    double minorAxis {};
};

struct annulusSimData
{
    int numMonteCarloSteps {};
    int outerRadius {};
    int innerRadius {};
    double majorAxis {};
    double minorAxis {};
};