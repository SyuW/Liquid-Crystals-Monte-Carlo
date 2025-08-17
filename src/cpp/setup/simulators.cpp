#include <cmath>
#include <random>
#include <string>
#include <map>
#include <vector>

#include "src/cpp/auxiliary/list"
#include "src/cpp/setup/constants.hpp"
#include "src/cpp/setup/overlap.hpp"
#include "src/cpp/setup/output.hpp"

void tuneAcceptanceRate(const double rate, double &stepXY, double &stepTh)
{
    /* Helper function for tuning acceptance rates throughout the duration of a simulation
     *
     * rate     - current acceptance rate of Monte Carlo steps
     * stepXY   - maximum step sizes for translation x and y directions
     * stepTh   - maximum angle of rotation
     */

    // if acceptance rate is too high, make step sizes larger
    // if rate is too low, make them smaller.
    if (rate <= 0.07)
    {
        stepTh *= 0.2;
        stepXY *= 0.2;
    }
    else if (rate > 0.07 && rate <= 0.17)
    {
        stepTh *= 0.35;
        stepXY *= 0.35;
    }
    else if (rate > 0.17 && rate <= 0.27)
    {
        stepTh *= 0.6;
        stepXY *= 0.6;
    }
    else if (rate > 0.27 && rate <= 0.37)
    {
        stepTh *= 0.75;
        stepXY *= 0.75;
    }
    else if (rate > 0.37 && rate <= 0.47)
    {
        stepTh *= 0.9;
        stepXY *= 0.9;
    }
    else if (rate > 0.47 && rate <= 0.57)
    {
        stepTh *= 1;
        stepXY *= 1;
    }
    else if (rate > 0.57 && rate <= 0.67)
    {
        stepTh *= 1.1;
        stepXY *= 1.1;
    }
    else if (rate > 0.67 && rate <= 0.77)
    {
        stepTh *= 1.35;
        stepXY *= 1.35;
    }
    else if (rate > 0.77 && rate <= 0.87)
    {
        stepTh *= 1.5;
        stepXY *= 1.5;
    }
    else if (rate > 0.87 && rate <= 0.97)
    {
        stepTh *= 1.65;
        stepXY *= 1.65;
    }
    else if (rate > 0.97)
    {
        stepTh *= 1.8;
        stepXY *= 1.8;
    }
}

Matrix boxHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                 const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                 Matrix posArray, const std::string outDir)
{
    /*
     * Hard particle Monte Carlo with hard box boundary conditions
     *
     * numParticles         - number of particles to simulate
     * numMonteCarloSteps   - number of Monte Carlo steps to perform
     * boxHeight            - height of box
     * boxWidth             - width of box
     * majorAxis            - major axis of ellipse particle
     * minorAxis            - minor axis of ellipse particle
     * posArray             - array of particle positions
     *
     */

    // seed a Mersenne Twister random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(-1, 1);
    double stepXY{0.5 * boxWidth};
    double stepTh{PI / 2};
    double proposedX{0};
    double proposedY{0};
    double proposedTh{0};
    int acceptedMoves{0};
    int totalMoves{0};
    int tunePeriod{numMonteCarloSteps / 50};
    int writeOutPeriod{1000};
    // main simulation loop
    for (int stepNo{1}; stepNo <= numMonteCarloSteps; ++stepNo)
    {
        // tune displacements to maintain acceptance rate
        if (stepNo % tunePeriod == 0)
        {
            tuneAcceptanceRate(static_cast<double>(acceptedMoves) / totalMoves, stepXY, stepTh);
        }
        // write out configurations periodically
        if (stepNo % writeOutPeriod == 0)
        {
            std::string outName{"positionsStep" + std::to_string(stepNo)};
            writeOutPositions(posArray, outDir + outName);
        }
        // attempt a Monte Carlo move for each particle
        for (int particleIndex1{0}; particleIndex1 < numParticles; ++particleIndex1)
        {
            // generate a trial position
            proposedX = posArray(particleIndex1, 0) + stepXY * uniform_dist(rng);
            proposedY = posArray(particleIndex1, 1) + stepXY * uniform_dist(rng);
            proposedTh = posArray(particleIndex1, 2) + stepTh * uniform_dist(rng);
            totalMoves += 1;
            // check for boundary overlaps
            if ((proposedX < 0 || proposedX > boxWidth) || (proposedY < 0 || proposedY > boxHeight) ||            // ellipse is outside
                checkBoundaryOverlapLine(0, 0, minorAxis, majorAxis, proposedX, proposedY, proposedTh) ||         // bottom boundary
                checkBoundaryOverlapLine(0, boxHeight, minorAxis, majorAxis, proposedX, proposedY, proposedTh) || // top boundary
                checkBoundaryOverlapVertical(0, minorAxis, majorAxis, proposedX, proposedY, proposedTh) ||        // left boundary
                checkBoundaryOverlapVertical(boxWidth, minorAxis, majorAxis, proposedX, proposedY, proposedTh))   // right boundary
            {
                continue;
            }
            // if no boundary overlaps, check for overlaps with other ellipse particles
            bool overlapVar{false};
            for (int particleIndex2{0}; particleIndex2 < numParticles; ++particleIndex2)
            {
                if (particleIndex2 == particleIndex1)
                {
                    continue;
                }
                double x2{posArray(particleIndex2, 0)};
                double y2{posArray(particleIndex2, 1)};
                double t2{posArray(particleIndex2, 2)};
                if (checkEllipseEllipseOverlap(proposedX, proposedY, x2, y2, proposedTh, t2, minorAxis, majorAxis))
                {
                    overlapVar = true;
                    break;
                }
            }
            if (!overlapVar)
            {
                posArray(particleIndex1, 0) = proposedX;
                posArray(particleIndex1, 1) = proposedY;
                posArray(particleIndex1, 2) = proposedTh;
                acceptedMoves += 1; // move has been accepted, so increment counter by one
            }
        }
    }
    std::cout << "Done simulation.\n";
    std::cout << "Final acceptance rate: " << static_cast<double>(acceptedMoves) / totalMoves << "\n";
    return posArray;
}

Matrix boxPeriodicBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps,
                                     const double boxHeight, const double boxWidth, const double majorAxis, const double minorAxis,
                                     Matrix posArray, const std::string outDir)
{
    /*
     * Hard Particle Monte Carlo with box periodic boundary conditions
     *
     * numParticles             - number of particles to simulate
     * numMonteCarloSteps       - number of Monte Carlo steps to perform
     * boxHeight                - height of box
     * boxWidth                 - width of box
     * majorAxis                - major axis of ellipse particle
     * minorAxis                - minor axis of ellipse particle
     * posArray                 - array of particle positions
     *
     */

    // seed a random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(-1, 1);
    // periodic image translation directions
    int translationDirections[9][2] =
        {
            {-1, +1}, {+0, +1}, {+1, +1}, {-1, +0}, {+0, +0}, {+1, +0}, {-1, -1}, {+0, -1}, {+1, -1}};
    double stepXY{0.5 * boxWidth};
    double stepTh{PI / 2};
    double proposedX{};
    double proposedY{};
    double proposedTh{};
    int acceptedMoves{0};
    int totalMoves{0};
    int tunePeriod{numMonteCarloSteps / 50};
    int writeOutPeriod{1000};
    // main simulation loop
    for (int stepNo{1}; stepNo <= numMonteCarloSteps; ++stepNo)
    {
        if (stepNo % tunePeriod == 0)
        {
            tuneAcceptanceRate(static_cast<double>(acceptedMoves) / totalMoves, stepXY, stepTh);
        }
        if (stepNo % writeOutPeriod == 0)
        {
            std::string outName{"positionsStep" + std::to_string(stepNo)};
            writeOutPositions(posArray, outDir + outName);
        }
        for (int particleIndex1{0}; particleIndex1 < numParticles; ++particleIndex1)
        {
            // generate a trial position
            proposedX = posArray(particleIndex1, 0) + stepXY * uniform_dist(rng);
            proposedY = posArray(particleIndex1, 1) + stepXY * uniform_dist(rng);
            proposedTh = posArray(particleIndex1, 2) + stepTh * uniform_dist(rng);
            // need to check overlaps with other ellipses under periodic boundary conditions
            bool overlapVar{false};
            for (int particleIndex2{0}; particleIndex2 < numParticles; ++particleIndex2)
            {
                if (particleIndex2 == particleIndex1)
                {
                    continue;
                }
                double x2{posArray(particleIndex2, 0)};
                double y2{posArray(particleIndex2, 1)};
                double t2{posArray(particleIndex2, 2)};
                for (int latticeTransIndex{0}; latticeTransIndex < 9; ++latticeTransIndex)
                {
                    double lX{boxWidth * translationDirections[latticeTransIndex][0]};
                    double lY{boxHeight * translationDirections[latticeTransIndex][1]};
                    // ellipse-ellipse overlap detected: exit the loop and reject
                    if (checkEllipseEllipseOverlap(proposedX + lX, proposedY + lY, x2, y2,
                                                   proposedTh, t2, minorAxis, majorAxis))
                    {
                        overlapVar = true;
                        break;
                    }
                }
                if (overlapVar)
                {
                    break;
                }
            }
            totalMoves += 1;
            if (!overlapVar)
            {
                // with periodic boundary conditions, need to reduce modulo lattice translation
                posArray(particleIndex1, 0) = (proposedX < 0) ? boxWidth - std::abs(fmod(proposedX, boxWidth)) : fmod(proposedX, boxWidth);
                posArray(particleIndex1, 1) = (proposedY < 0) ? boxHeight - std::abs(fmod(proposedY, boxHeight)) : fmod(proposedY, boxHeight);
                posArray(particleIndex1, 2) = proposedTh;
                acceptedMoves += 1;
            }
        }
    }
    return posArray;
}

Matrix circleHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double boundaryRadius,
                                    const double majorAxis, const double minorAxis,
                                    Matrix posArray, const std::string outDir)
{

    /*
     * Hard particle Monte Carlo with hard circle boundary conditions
     *
     * numParticles         - number of particles to simulate
     * numMonteCarloSteps   - number of Monte Carlo steps to perform
     * boundaryRadius       - radius of circle boundary
     * majorAxis            - major axis of ellipse particle
     * minorAxis            - minor axis of ellipse particle
     * posArray             - array of particle positions
     *
     */

    // seed a random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(-1, 1);
    double stepXY{0.5 * boundaryRadius};
    double stepTh{PI / 2};
    double proposedX{};
    double proposedY{};
    double proposedTh{};
    int acceptedMoves{0};
    int totalMoves{0};
    int writeOutPeriod{1000};
    int tunePeriod{numMonteCarloSteps / 50};
    // each Monte Carlo step (or sweep) involves attempting a move for each particle
    for (int stepNo{1}; stepNo <= numMonteCarloSteps; ++stepNo)
    {
        // tune displacements to maintain acceptance rate
        if (stepNo % tunePeriod == 0)
        {
            tuneAcceptanceRate(static_cast<double>(acceptedMoves) / totalMoves, stepXY, stepTh);
        }
        // write out configurations periodically
        if (stepNo % writeOutPeriod == 0)
        {
            std::string outName{"positionsStep" + std::to_string(stepNo)};
            writeOutPositions(posArray, outDir + outName);
        }
        for (int particleIndex1{0}; particleIndex1 < numParticles; ++particleIndex1)
        {
            // generate trial position
            proposedX = posArray(particleIndex1, 0) + stepXY * uniform_dist(rng);
            proposedY = posArray(particleIndex1, 1) + stepXY * uniform_dist(rng);
            proposedTh = posArray(particleIndex1, 2) + stepTh * uniform_dist(rng);
            totalMoves += 1;
            // check boundary overlaps
            if (pow(proposedX, 2) + pow(proposedY, 2) > pow(boundaryRadius, 2) ||                                   // cannot be outside
                checkBoundaryOverlapCircle(boundaryRadius, minorAxis, majorAxis, proposedX, proposedY, proposedTh)) // cannot overlap boundary
            {
                continue;
            }
            bool overlapVar{false};
            for (int particleIndex2{0}; particleIndex2 < numParticles; ++particleIndex2)
            {
                if (particleIndex2 == particleIndex1)
                {
                    continue;
                }
                double x2{posArray(particleIndex2, 0)};
                double y2{posArray(particleIndex2, 1)};
                double t2{posArray(particleIndex2, 2)};
                // ellipse-ellipse overlap detected: exit the loop and reject
                if (checkEllipseEllipseOverlap(proposedX, proposedY, x2, y2, proposedTh, t2, minorAxis, majorAxis))
                {
                    overlapVar = true;
                    break;
                }
            }
            if (!overlapVar)
            {
                posArray(particleIndex1, 0) = proposedX;
                posArray(particleIndex1, 1) = proposedY;
                posArray(particleIndex1, 2) = proposedTh;
                acceptedMoves += 1;
            }
        }
    }
    std::cout << "Done simulation.\n";
    std::cout << "Final acceptance rate: " << static_cast<double>(acceptedMoves) / totalMoves << "\n";
    return posArray;
}

Matrix annulusHardBoundaryMonteCarlo(const int numParticles, const int numMonteCarloSteps, const double innerRadius,
                                     const double outerRadius, const double majorAxis, const double minorAxis, Matrix posArray)
{

    // seed a random number generator
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<> uniform_dist(-1, 1);

    return posArray;
}
