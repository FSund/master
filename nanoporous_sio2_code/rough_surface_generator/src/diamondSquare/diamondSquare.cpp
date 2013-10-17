#include <src/diamondSquare/diamondSquare.h>

DiamondSquare::DiamondSquare(const int power2, long idum, const int RNG, const bool PBC):
    power2(power2),
    systemSize(pow(2.0, power2) + 1),
    idum(idum),
    R(zeros<mat>(systemSize,systemSize)),
    RNG(RNG),
    PBC(PBC) {
    if (RNG > 2) {
        cout << "RNG too large, have only implemented 3 random number generators" << endl;
        exit(1);
    }
}

mat DiamondSquare::generate(const double H, const double minZValue, const double maxZValue) {

//    R(0,0)                       = ran2(&idum) - 0.5;
//    R(0,systemSize-1)            = ran2(&idum) - 0.5;
//    R(systemSize-1,0)            = ran2(&idum) - 0.5;
//    R(systemSize-1,systemSize-1) = ran2(&idum) - 0.5;
    if (PBC) { // need the same value in the corners if we are using periodic boundaries
        R(0,0)                       = randu<double>() - 0.5;
        R(0,systemSize-1)            = R(0,0);
        R(systemSize-1,0)            = R(0,0);
        R(systemSize-1,systemSize-1) = R(0,0);
    } else {
        R(0,0)                       = randu<double>() - 0.5;
        R(0,systemSize-1)            = randu<double>() - 0.5;
        R(systemSize-1,0)            = randu<double>() - 0.5;
        R(systemSize-1,systemSize-1) = randu<double>() - 0.5;
    }

    runDiamondSquare(R, H);

    if (RNG == 0 && PBC == 0) {
        return R;
    }
    // normalize to minZ maxZ
//    for (uint i = 0; i < system)

    return R;
}

void DiamondSquare::runDiamondSquare(mat& R, const double H) {

    double RNGstddv;

    uint stepLength = systemSize-1;
    uint halfStepLength = stepLength/2;

    for (uint iteration = 0; iteration < 2*power2; iteration+=2) { // Doing two iterations inside the loop
        uint nSteps = (systemSize-1)/stepLength;

        // Squares
        RNGstddv = sqrt(1.0/pow(2.0, 2.0*H*iteration));
        for (uint x = 0; x < nSteps*stepLength; x += stepLength) {
            for (uint y = 0; y < nSteps*stepLength; y += stepLength) {
                square(x, y, stepLength, halfStepLength, RNGstddv, R);
            }
        }
//        cout << "R after squares" << endl << R << endl;

        // Diamonds
        RNGstddv = sqrt(1.0/pow(2.0, 2.0*H*(iteration+1)));
        for (uint x = 0; x < nSteps*stepLength; x += stepLength) {
            for (uint y = 0; y < nSteps*stepLength; y += stepLength) {
                diamond(x, y, stepLength, halfStepLength, RNGstddv, R);
            }
        }
//        cout << "R after diamonds" << endl << R << endl;

        if (PBC) {
            for (uint i = 0; i < nSteps; i++) {
                uint idx = halfStepLength + i*stepLength;
                R(idx, systemSize-1) = R(idx,0);
                R(systemSize-1, idx) = R(0,idx);
            }
//            cout << "R after PBC" << endl << R << endl;
        } else {
            // We have to do the bottom and right edge diamonds manually

            // Bottom edge diamonds
            for (uint y = halfStepLength; y < nSteps*stepLength; y += stepLength) {
                uint x = (nSteps-1)*stepLength + halfStepLength;
                bottomEdgeDiamonds(x, y, halfStepLength, RNGstddv, R);
            }
//            cout << "R after bottom diamonds" << endl << R << endl;

            // Right edge diamonds
            for (uint x = halfStepLength; x < nSteps*stepLength; x+= stepLength) {
                uint y = (nSteps-1)*stepLength + halfStepLength;
                rightEdgeDiamonds(x, y, halfStepLength, RNGstddv, R);
            }
//            cout << "R after right diamonds" << endl << R << endl;
        }

        stepLength /= 2;
        halfStepLength /= 2;
    }
}

void DiamondSquare::square(
        const uint x,
        const uint y,
        const uint stepLength,
        const uint halfStepLength,
        const double RNGstddv,
        mat &R) {

    double average = 0.25*(R(x, y) + R(x+stepLength, y) + R(x, y+stepLength) + R(x+stepLength, y+stepLength));
    R(x + halfStepLength, y + halfStepLength) = average + random()*RNGstddv; // change to = instead of +=
}

void DiamondSquare::diamond(
        const uint x,
        const uint y,
        const uint stepLength,
        const uint halfStepLength,
        const double RNGstddv,
        mat &R) {

    double average;

    // Point centered at left edge of square
    if (y == 0) { // At left edge of system
        if (PBC) {
            average = 0.25*(
                R(x, y) +
                R(x+stepLength, y) +
                R(x+halfStepLength, y+halfStepLength) +
                R(x+halfStepLength, systemSize - 1 - halfStepLength));
        } else {
            average = 0.33333333333333333333333*(
                R(x, y) +
                R(x+stepLength, y) +
                R(x+halfStepLength, y+halfStepLength));
        }
    } else { // Inside the system -- nothing to worry about
        average = 0.25*(
            R(x, y) +
            R(x+stepLength, y) +
            R(x+halfStepLength, y+halfStepLength) +
            R(x+halfStepLength, y-halfStepLength));
    }
    R(x + halfStepLength, y) = average + random()*RNGstddv;

    // Point centered at top edge of square
    if (x == 0) { // At top edge of system
        if (PBC) {
            average = 0.25*(
                R(x, y) +
                R(x, y+stepLength) +
                R(x+halfStepLength, y+halfStepLength) +
                R(systemSize - 1 - halfStepLength, y+halfStepLength));
        } else {
            average = 0.33333333333333333333333*(
                R(x, y) +
                R(x, y+stepLength) +
                R(x+halfStepLength, y+halfStepLength));
        }
    } else { // Inside the system -- nothing to worry about
        average = 0.25*(
            R(x, y) +
            R(x, y+stepLength) +
            R(x+halfStepLength, y+halfStepLength) +
            R(x-halfStepLength, y+halfStepLength));
    }
    R(x, y + halfStepLength) = average + random()*RNGstddv;
}

void DiamondSquare::bottomEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, mat& R) {

    double average = 0.33333333333333333333333*(
        R(x, y) +
        R(x+halfStepLength, y+halfStepLength) +
        R(x+halfStepLength, y-halfStepLength));
    R(x + halfStepLength, y) = average + random()*RNGstddv;
}

void DiamondSquare::rightEdgeDiamonds(const uint x, const uint y, const uint halfStepLength, const double RNGstddv, mat& R) {

    double average = 0.33333333333333333333333*(
        R(x, y) +
        R(x+halfStepLength, y+halfStepLength) +
        R(x-halfStepLength, y+halfStepLength));
    R(x, y + halfStepLength) = average + random()*RNGstddv;
}
