#include <iostream>
#include <cstdlib>      // atoi, atof, atol
#include <armadillo>
#include <src/diamondSquare/diamondSquare.h>
#include <src/mesher/mesher.h>

using namespace std;
using namespace arma;

void normalize(mat &R, double position);

int main(int nArgs, const char *argv[]) {
    string filename;
    double positionOfBottomSurface, positionOfTopSurface;

    int power2;
    double H;
    vec corners = {0.0};
    double sigma;
    bool addition = true;
    bool PBC = true;
    int RNG = 2;
    int seed;

    if (nArgs < 3) {
        cout << "Usage: ./fracture_generator  power2  filename  optional:(H  position_of_bottom_surface  position_of_top_surface  initial_RNG_stddv  seed[unsigned int])" << endl;
        exit(1);
    }

    // arguments that are needed
    power2   = atoi(argv[1]);
    filename = argv[2];

    // argument that have default values
    H                       = nArgs > 3 ? atof(argv[3]) : 0.75;
    positionOfBottomSurface = nArgs > 4 ? atof(argv[4]) : 0.2;
    positionOfTopSurface    = nArgs > 5 ? atof(argv[5]) : 0.8;
    sigma                   = nArgs > 6 ? atof(argv[6]) : 0.3;
    seed                    = nArgs > 7 ? atol(argv[7]) : 1;

    cout << "--- Diamond-square settings --------------------" << endl;
    cout << "power2 = " << power2  << endl;
    cout << "filename = " << filename << endl;
    cout << "H (Hurst exponent) = " << H << endl;
    cout << "position_of_bottom_surface = " << positionOfBottomSurface << endl;
    cout << "position_of_top_surface    = " << positionOfTopSurface << endl;
    cout << "initial_RNG_stddv          = " << sigma << endl;
//    cout << "PBC  = " << std::boolalpha << PBC << std::noboolalpha << endl;
    cout << "seed = " << seed << endl;
//    cout << "RNG  = " << RNG << " (0 == no RNG, 1 == uniform, 2 == standard normal distribution)" << endl;
    cout << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
    cout << "------------------------------------------------" << endl;

    srand(seed); // setting the seed of the RNG for both C++'s rand()/srand() and Armadillo's randu()/randn()

    DiamondSquare generator;

    seed = rand();
    mat bottomHeighmap = generator.generate(power2, H, corners, seed, sigma, addition, PBC, RNG);
    normalize(bottomHeighmap, positionOfBottomSurface);

    seed = rand();
    mat topHeightmap = generator.generate(power2, H, corners, seed, sigma, addition, PBC, RNG);
    normalize(topHeightmap, positionOfTopSurface);

    Mesher mesher;
    mesher.mesh(topHeightmap, bottomHeighmap);
    mesher.printToMsh(filename);

    return 0;
}

void normalize(mat &R, double position) {
//    // normalize to range [minZ maxZ]
//    zMin = positionOfBottomSurface - surfaceDeltaZ/2.0;
//    zMax = positionOfBottomSurface + surfaceDeltaZ/2.0;
//    minValue = min(min(bottomHeighmap));
//    normFactor = 1.0/(max(max(bottomHeighmap)) - minValue);
//    for (uint i = 0; i < bottomHeighmap.n_rows; i++) {
//        for (uint j = 0; j < bottomHeighmap.n_cols; j++) {
//            bottomHeighmap(i,j) = (bottomHeighmap(i,j) - minValue)*normFactor*(zMax - zMin) + zMin;
//        }
//    }

    // move
    double zAverage = mean(mean(R));
    for (uint i = 0; i < R.n_rows; i++) {
        for (uint j = 0; j < R.n_cols; j++) {
            R(i,j) = R(i,j) - zAverage + position;
        }
    }
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
