#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>      // atoi, atof, atol
#include <cmath>        // sqrt, pow

#include <armadillo>

#include <src/diamondSquare/diamondSquare.h>
#include <src/mesher/mesher.h>
//#include <src/fileIO.h>

using namespace std;
using namespace arma;

template<typename T> mat convertVectorToArma(const vector<vector<T> > vector);
void displace(mat &R, double position);
void normalize(mat &R, double position, double range);

int main(int nArgs, const char *argv[]) {

    // Arguments for both modes
    string mode;
    bool generate_mode;
    string filename;
    double positionOfBottomSurface, positionOfTopSurface, surfaceDelta;

    // Arguments for generate mode
    int power2;
    double H;
    vector<double> corners; // Empty vector, so DiamondSquare generates random corners
    double sigma;
    double randomFactor;
    bool addition;
    bool PBC;
    int RNG;
    int seed;

    // Arguments for load mode
    string bottomSurfaceFilename, topSurfaceFilename;

    if (nArgs < 2) {
        cout << "Usage: ./fracture_generator "
             << "mode[\"generate\"|\"load\"] " << endl
             << "    arguments for \"generate\" mode: " << endl
             << "        power2  output_filename " << endl
             << "        optional: H  position_of_bottom_surface  position_of_top_surface " << endl
             << "                  surface_delta  initial_RNG_stddv  RNG_range_reduction_factor " << endl
             << "    arguments for \"load\" mode: " << endl
             << "        bottom_surface_filename  top_surface_filename  output_filename " << endl
             << "        optional: position_of_bottom_surface  position_of_top_surface  surface_delta "
             << endl;
        exit(1);
    }

    int arg_counter = 1;
    mode = argv[arg_counter++];
    if (mode == "generate") {
        if (nArgs < 4) {
            cout << "Error: Generate mode needs at least \"power2\" and \"output_filename\" arguments. Aborting.!" << endl;
            exit(1);
        } else if (nArgs > 10) {
            cout << "Error: Too many arguments to generate mode. Aborting!" << endl;
            exit(1);
        }
        generate_mode = true;

        // Arguments that are needed
        power2   = atoi(argv[arg_counter++]);
        filename = argv[arg_counter++];

        // argument that have default values
        H                       = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.75;
        positionOfBottomSurface = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.2;
        positionOfTopSurface    = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.8;
        surfaceDelta            = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.2;
        sigma                   = nArgs > arg_counter ? atof(argv[arg_counter++]) : 1.0;
        randomFactor            = nArgs > arg_counter ? atof(argv[arg_counter++]) : 1.0/sqrt(2);
        addition = true;
        PBC = true;
        RNG = 2;
        seed = 1;

        cout << "--- Diamond-square settings --------------------------------------------" << endl;
        cout << "power2                     = " << power2  << endl;
        cout << "filename                   = " << filename << endl;
        cout << "H (Hurst exponent)         = " << H << endl;
        cout << "position_of_bottom_surface = " << positionOfBottomSurface << endl;
        cout << "position_of_top_surface    = " << positionOfTopSurface << endl;
        cout << "surface_delta              = " << surfaceDelta << endl;
        cout << "initial_RNG_stddv          = " << sigma << endl;
        cout << "random factor              = " << randomFactor << endl;
        cout << "addition                   = " << std::boolalpha << addition << std::noboolalpha << endl;
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
            ofile << "  random corners              = true" << endl;
        }
        ofile << "  addition                    = " << addition << endl;
        ofile << "  PBC                         = " << PBC << endl;
        ofile << "  RNG                         = " << RNG << endl;
        ofile << "  seed                        = " << seed << endl;
        ofile << endl;
        ofile << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
        ofile << "------------------------------------------------" << endl;

    } else if (mode == "load") {
        if (nArgs < 4) {
            cout << "Error: Load mode needs at least \"bottom_surface_filename\", \"top_surface_filename\" and \"output_filename\" arguments. Aborting.!" << endl;
            exit(1);
        }
        generate_mode = false;

        // Arguments that are needed
        bottomSurfaceFilename = argv[arg_counter++];
        topSurfaceFilename = argv[arg_counter++];
        filename = argv[arg_counter++];

        // argument that have default values
        positionOfBottomSurface = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.2;
        positionOfTopSurface    = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.8;
        surfaceDelta            = nArgs > arg_counter ? atof(argv[arg_counter++]) : 0.2;
    } else {
        cout << "Error: Mode \"" << mode << "\" not supported. Aborting!" << endl;
        exit(1);
    }

    mat topHeightmap, bottomHeightmap;

    if (generate_mode == true)  {
        DiamondSquare generator(power2, RNG, seed);

        bottomHeightmap = convertVectorToArma(generator.generate(H, corners, sigma, randomFactor, addition, PBC));
        normalize(bottomHeightmap, positionOfBottomSurface, surfaceDelta);

        topHeightmap = convertVectorToArma(generator.generate(H, corners, sigma, randomFactor, addition, PBC));
        normalize(topHeightmap, positionOfTopSurface, surfaceDelta);
    } else {
        if (!bottomHeightmap.load(bottomSurfaceFilename, raw_ascii)) {
            cout << "Error loading file \"" << bottomSurfaceFilename << "\", please check the file. Aborting!" << endl;
            exit(1);
        }
        normalize(bottomHeightmap, positionOfBottomSurface, surfaceDelta);

        if (!topHeightmap.load(topSurfaceFilename, raw_ascii)) {
            cout << "Error loading file \"" << topSurfaceFilename << "\", please check the file. Aborting!" << endl;
            exit(1);
        }
        normalize(topHeightmap, positionOfTopSurface, surfaceDelta);
    }

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
