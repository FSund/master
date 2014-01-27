#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>      // atoi, atof, atol
#include <cmath>        // sqrt, pow

#include <armadillo>

#include <src/diamondSquare/diamondSquare.h>
#include <src/mesher/mesher.h>

using namespace std;
using namespace arma;

template<typename T> mat convertVectorToArma(const vector<vector<T> > vector);
void displace(mat &R, double position);
void normalize(mat &R, double position, double range);

int main(int nArgs, const char *argv[]) {
    string filename;
    double positionOfBottomSurface, positionOfTopSurface, surfaceDelta;

    int power2;
    double H;
    vector<double> corners; // Empty vector, so DiamondSquare generates random corners
    double sigma;
    double randomFactor;
    bool addition = true;
    bool PBC = true;
    int RNG = 2;
    int seed = 1;

    if (nArgs < 3) {
        cout << "Usage: ./fracture_generator  power2  filename  optional:(H  position_of_bottom_surface  position_of_top_surface  surface_delta  initial_RNG_stddv  RNG_range_reduction_factor)" << endl;
        exit(1);
    }

    // arguments that are needed
    power2   = atoi(argv[1]);
    filename = argv[2];

    // argument that have default values
    H                       = nArgs > 3 ? atof(argv[3]) : 0.75;
    positionOfBottomSurface = nArgs > 4 ? atof(argv[4]) : 0.2;
    positionOfTopSurface    = nArgs > 5 ? atof(argv[5]) : 0.8;
    surfaceDelta            = nArgs > 6 ? atof(argv[6]) : 0.2;
    sigma                   = nArgs > 7 ? atof(argv[7]) : 1.0;
    randomFactor            = nArgs > 8 ? atof(argv[8]) : 1.0/sqrt(2);
//    seed                    = nArgs > 9 ? atol(argv[9]) : 1;

    cout << "--- Diamond-square settings --------------------------------------------" << endl;
    cout << "power2                     = " << power2  << endl;
    cout << "filename                   = " << filename << endl;
    cout << "H (Hurst exponent)         = " << H << endl;
    cout << "position_of_bottom_surface = " << positionOfBottomSurface << endl;
    cout << "position_of_top_surface    = " << positionOfTopSurface << endl;
    cout << "surface_delta              = " << surfaceDelta << endl;
    cout << "initial_RNG_stddv          = " << sigma << endl;
    cout << "random factor              = " << randomFactor << endl;
    cout << "PBC                        = " << std::boolalpha << PBC << std::noboolalpha << endl;
    cout << "seed                       = " << seed << endl;
    cout << "RNG                        = " << RNG << " (0 == no RNG, 1 == uniform, 2 == standard normal distribution)" << endl;
    cout << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
    cout << "------------------------------------------------------------------------" << endl;

    // Printing settings to <filename>-fracture_generator-settings.txt
    stringstream fNameGen;
    fNameGen << filename << "-fracture_generator-settings.txt";
    ofstream ofile(fNameGen.str());

    ofile << "--- Diamond-square settings --------------------" << endl;
    ofile << "Input arguments: " << endl;
    ofile << "  power2                      = " << power2  << endl;
    ofile << "  filename                    = " << filename << endl;
    ofile << "  H (Hurst exponent)          = " << H << endl;
    ofile << "  position_of_bottom_surface  = " << positionOfBottomSurface << endl;
    ofile << "  position_of_top_surface     = " << positionOfTopSurface << endl;
    ofile << "  surface_delta               = " << surfaceDelta << endl;
    ofile << "  initial_RNG_stddv           = " << sigma << endl;
    ofile << "  RNG_range_reduction_factor  = " << randomFactor << endl;
    ofile << endl << "Constant settings:" << endl;
    if (corners.size() > 0) {
        ofile << "  corners                     = ";
        for (uint i = 0; i < corners.size(); i++) {
            ofile << corners[i] << " ";
        }
        ofile << endl;
    } else {
        ofile << "  Using random corners" << endl;
    }
    ofile << "  addition                    = " << addition << endl;
    ofile << "  PBC                         = " << PBC << endl;
    ofile << "  RNG                         = " << RNG << endl;
    ofile << "  seed                        = " << seed << endl;
    ofile << endl;
    ofile << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
    ofile << "------------------------------------------------" << endl;

    DiamondSquare generator(power2, RNG, seed);

    mat bottomHeightmap = convertVectorToArma(generator.generate(H, corners, sigma, randomFactor, addition, PBC));
    normalize(bottomHeightmap, positionOfBottomSurface, surfaceDelta);

    mat topHeightmap = convertVectorToArma(generator.generate(H, corners, sigma, randomFactor, addition, PBC));
    normalize(topHeightmap, positionOfTopSurface, surfaceDelta);

    Mesher mesher;
    mesher.mesh(topHeightmap, bottomHeightmap);
    mesher.printToMsh(filename);

    cout << "max = " << max(max(topHeightmap)) << " " << max(max(bottomHeightmap)) << endl;
    cout << "min = " << min(min(topHeightmap)) << " " << min(min(bottomHeightmap)) << endl;

    return 0;
}

template<typename T>
mat convertVectorToArma(const vector<vector<T> > vector) {
    mat matrix(vector.size(), vector[0].size());
    for (uint i = 0; i < vector.size(); i++) {
        for (uint j = 0; j < vector[i].size(); j++) {
            matrix(i,j) = vector[i][j];
        }
    }
    return matrix;
}

void normalize(mat &R, double position, double range) {
    // normalize to range
    double zMin = position - range/2.0;
    double zMax = position + range/2.0;

    double minValue = min(min(R));
    double maxValue = max(max(R));
    double normFactor = 1.0/(maxValue - minValue);
    for (uint i = 0; i < R.n_rows; i++) {
        for (uint j = 0; j < R.n_cols; j++) {
            R(i,j) = (R(i,j) - minValue)*normFactor*(zMax - zMin) + zMin;
        }
    }
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
