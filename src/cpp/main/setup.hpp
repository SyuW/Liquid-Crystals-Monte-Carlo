#include "../auxiliary/list"
#include "../setup/overlap.cpp"
#include "../setup/simulators.cpp"
#include "../setup/input.cpp"
#include "../setup/output.cpp"
#include "../utils/initial.cpp"


// class Simulation
// {
// public:
//     Simulation(void);
//     ~Simulation(void);

//     // overlap functions
//     bool checkEllipseEllipseOverlap(double x1, double y1, double x2, double y2, 
//                                     double theta1, double theta2,
//                                     double minorAxis, double majorAxis, bool verbose);
                                    
//     bool checkBoundaryOverlapCircle(double R, double minorAxis, double majorAxis,
//                                     double xc, double yc, double theta, bool debug);

//     // writing out
//     void writeOutPositionsCircle(const double majorAxis, const double minorAxis, const double boundaryRadius,
//                                  Matrix& posArray, const std::string fileName);

//     void writeOutPositionsBox(const double majorAxis, const double minorAxis,
//                               const double length, const double width, Matrix& posArray, const std::string fileName);
    
//     // simulation functions
//     void tuneAcceptanceRate(const double rate, double& stepXY, double& stepTh);
//     void circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
//                                       const double boundaryRadius, const double majorAxis, const double minorAxis);

    
//     // getting input

// private:

//     // system setup
//     int numParticles {100};
//     int numMonteCarloSteps {20000};
//     double boundaryRadius {25};
//     double majorAxis {7};
//     double minorAxis {1};
// };

// Simulation::Simulation(void)
// {

// }

// Simulation::~Simulation(void)
// {

// }