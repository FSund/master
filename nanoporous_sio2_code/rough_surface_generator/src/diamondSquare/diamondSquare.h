#ifndef DIAMONDSQUARE_H
#define DIAMONDSQUARE_H

#include <iostream>
#include <armadillo>
//#include "lib.h"

using namespace std;
using namespace arma;

class DiamondSquare {
public:
    DiamondSquare(const int power2, long idum, const int RNG = 1, const bool PBC = false);
    mat generate(const double H, const double minZValue, const double maxZValue);

private:
    void runDiamondSquare(mat &R, const double H);
    void square(const uint x, const uint y, const uint squareSize, const uint squareSizeHalf, const double RNGstddv, mat &R);
    void diamond(const uint x, const uint y, const uint stepLength, const uint halfStepLength, const double RNGstddv, mat &R);

    void bottomEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, mat &R);
    void rightEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, mat &R);

    double random();
    double gaussianDeviate(long *seed);

    uint power2, systemSize;
    long idum;
    mat R;
    const int RNG;
    const bool PBC;

//    uint xpos[4];
//    uint ypos[4];
//    double sum;
//    uint nPoints;
};

//inline double DiamondSquare::random(long *seed)
inline double DiamondSquare::random()
{
    /* returns random number with mean 0 */
    if (RNG == 0) {
        return 0.0;
    } else if (RNG == 1) {
//        return (ran2(seed) - 0.5);
        return (randu<double>() - 0.5); // uniform distribution in [-0.5,0.5]
    } else if (RNG == 2) {
//        return gaussianDeviate(seed);
        return randn<double>();
    } else {
        return NAN;
    }
}

//inline double DiamondSquare::gaussianDeviate(long *seed)
//{
//    /* returns normally distributed random number (mean 0, stddev 1) */
//    double R, randomNormal;
//    // Box-Muller transform
//    R = sqrt(-2.0*log(ran2(seed)));
//    randomNormal = R*cos(2.0*pi*ran2(seed));

//    return randomNormal;
//}

#endif // DIAMONDSQUARE_H
