#include <iostream>
#include <cstdlib>      // atoi, atof, atol
#include <armadillo>
#include <src/diamondSquare/diamondSquare.h>
#include <src/mesher/mesher.h>

using namespace std;
using namespace arma;

void parseArgs(int nArgs, const char* argv[], int &power2, string &filename, double &H, double &a, double& b, double &surfaceDeltaZ, bool &PBC, long &idum, int &RNG);

int main(int nArgs, const char *argv[]) {
    int power2;
    string filename;
    double H;
    double positionOfBottomSurface, positionOfTopSurface, surfaceDeltaZ;
    bool PBC;
    long idum;
    int RNG;

    parseArgs(nArgs, argv, power2, filename, H, positionOfBottomSurface, positionOfTopSurface, surfaceDeltaZ, PBC, idum, RNG);

    cout << "--- Diamond-square settings --------------------" << endl;
    cout << "power2 = " << power2  << endl;
    cout << "filename = " << filename << endl;
    cout << "H (Hurst exponent) = " << H << endl;
    cout << "position_of_bottom_surface = " << positionOfBottomSurface << endl;
    cout << "position_of_top_surface    = " << positionOfTopSurface << endl;
    cout << "surface_delta_z            = " << surfaceDeltaZ << endl;
    cout << "PBC = " << std::boolalpha << PBC << std::noboolalpha << endl;
    cout << "idum = " << idum << endl;
    cout << "RNG = " << RNG << " (0 == no RNG, 1 == uniform, 2 == standard normal distribution)" << endl;
    cout << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
    cout << "------------------------------------------------" << endl;

    srand(idum); // setting the seed of the RNG for both C++'s rand()/srand() and Armadillo's randu()/randn()

    double zMin, zMax;

    DiamondSquare generator;
    zMin = positionOfBottomSurface - surfaceDeltaZ/2.0;
    zMax = positionOfBottomSurface + surfaceDeltaZ/2.0;
    mat bottomHeighmap = generator.generate(power2, H, zMin, zMax, PBC, RNG);
    zMin = positionOfTopSurface - surfaceDeltaZ/2.0;
    zMax = positionOfTopSurface + surfaceDeltaZ/2.0;
    mat topHeightmap = generator.generate(power2, H, zMin, zMax, PBC, RNG);

    Mesher mesher;
    mesher.mesh(topHeightmap, bottomHeighmap);
    mesher.printToMsh(topHeightmap, bottomHeighmap, filename);

    return 0;
}

void parseArgs(
        int nArgs,
        const char *argv[],
        int &power2,
        string &filename,
        double &H,
        double &positionOfBottomSurface,
        double &positionOfTopSurface,
        double &surfaceDeltaZ,
        bool &PBC,
        long &idum,
        int &RNG) {

    if (nArgs < 3) {
        cout << "Usage: ./diamondSquare  power2  filename  optional:(H[1,2]  position_of_bottom_surface  position_of_top_surface  surface_delta_z  PBC[0|1]  idum[unsigned int]  RNG[0|1|2])" << endl;
        exit(1);
    }

    // arguments that are needed
    power2   = atoi(argv[1]);
    filename = argv[2];

    // argument that have default values
    H         = nArgs > 3 ? atof(argv[3]) : 1.5;
    positionOfBottomSurface = nArgs > 4 ? atof(argv[4]) : 0.2;
    positionOfTopSurface    = nArgs > 5 ? atof(argv[5]) : 0.8;
    surfaceDeltaZ           = nArgs > 6 ? atof(argv[6]) : 0.3;
    PBC       = nArgs > 7 ? atoi(argv[7]) : true;
    idum      = nArgs > 8 ? atol(argv[8]) : 1;
    RNG       = nArgs > 9 ? atoi(argv[9]) : 2;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
