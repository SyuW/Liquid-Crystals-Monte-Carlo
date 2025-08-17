#include "src/cpp/auxiliary/list"
#include <cmath>

// pair correlation function
// can compare the pair correlation function with literature for hard disk case
// should see characteristic decaying oscillations ~ exp(-k*x) * cos(l * x)
Vector computePairCorrelationFunc(Matrix &posArray, int numBins, double rMin, double rMax)
{
    Vector pairCorrelationFunc(numBins);
    double numParticles{posArray.getNumberOfRows()};
    double dr{(rMax - rMin) / numBins};

    for (int refIndex{0}; refIndex < numParticles; ++refIndex)
    {
        double x1{posArray(refIndex, 0)};
        double y1{posArray(refIndex, 1)};
        for (int index2{0}; index2 < numParticles; ++index2)
        {
            if (index2 == refIndex)
            {
                continue;
            }
            double x2{posArray(index2, 0)};
            double y2{posArray(index2, 1)};
            double dist{sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2))};
            if (dist < rMin || dist > rMax)
            {
                continue;
            }
        }
    }
}

// structure factor: for detecting positional order: connected to pair correlation function
// via a Fourier transformation

// nematic order parameters, such as the Q-tensor: for detecting orientational order without reference
// to an existing nematic director

// hexatic order parameter...
// Hexactic order parameter is useful for hard disk system case: a = b = sigma: for orientational, but not spatial order

// other order parameters for detecting rotational order for confinements which are not simply
// connected, such as annulus.. in order to capture and classify discrete rotational symmetries present (dihedral group)