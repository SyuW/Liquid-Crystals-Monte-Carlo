#ifndef SIMULATORSHEADERDEF
#define SIMULATORSHEADERDEF

#include "src/cpp/auxiliary/list"

// Helper function for tuning acceptance rates throughout duration of a simulation.
// If acceptance rate is too high, make step sizes larger.
// If rate is too low, make them smaller.
void tuneAcceptanceRate(const double rate, double &stepXY, double &stepTh);

// Hard particle Monte Carlo with box hard boundary conditions.
Matrix boxHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                 const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                 Matrix posArray, const std::string outDir);

// Hard particle Monte Carlo with box periodic boundary conditions.
Matrix boxPeriodicBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                     const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                     Matrix posArray, const std::string outDir);

// Hard particle Monte Carlo with hard circle boundary conditions.
Matrix circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double boundaryRadius,
                                    const double majorAxis, const double minorAxis,
                                    Matrix posArray, const std::string outDir);

#endif