#include <iostream>
#include <cstdlib>      // atoi, atof, atol
#include <armadillo>
#include <src/diamondSquare/diamondSquare.h>
#include <src/heightmapMesher/heightmapMesher.h>

using namespace std;
using namespace arma;

void parseArgs(int nArgs, const char* argv[], int &power2, string &filename, double &H, double &corners, double& maxZValue, bool &PBC, long &idum, int &RNG);

int main(int nArgs, const char *argv[]) {
    int power2;
    string filename;
    double H;
    double minZValue;
    double maxZValue;
    bool PBC;
    long idum;
    int RNG;

    parseArgs(nArgs, argv, power2, filename, H, minZValue, maxZValue, PBC, idum, RNG);

    cout << "--- Diamond-square settings --------------------" << endl;
    cout << "power2 = " << power2  << endl;
    cout << "filename = " << filename << endl;
    cout << "H (Hurst exponent) = " << H << endl;
    cout << "minZValue = " << minZValue << endl;
    cout << "maxZValue = " << maxZValue << endl;
    cout << "PBC = " << std::boolalpha << PBC << std::noboolalpha << endl;
    cout << "idum = " << idum << endl;
    cout << "RNG = " << RNG << " (0 == no RNG, 1 == uniform, 2 == standard normal distribution)" << endl;
    cout << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
    cout << "------------------------------------------------" << endl;

    srand(idum); // setting the seed of the RNG for both C++'s rand()/srand() and Armadillo's randu()/randn()

    DiamondSquare generator(power2, idum, RNG, PBC);
    mat R = generator.generate(H, minZValue, maxZValue);

//    cout << "R = " << endl << R << endl;

    HeightmapMesher mesher;
    mesher.mesh(R, filename);

//    R.save("rmat.dat", raw_ascii);

    return 0;
}

void parseArgs(
        int nArgs,
        const char *argv[],
        int &power2,
        string &filename,
        double &H,
        double &minZValue,
        double &maxZValue,
        bool &PBC,
        long &idum,
        int &RNG) {

    if (nArgs < 3) {
        cout << "Usage: ./diamondSquare  power2  filename  optional:(H[1,2]  min_value  max_value  PBC[0|1]  idum[unsigned int]  RNG[0|1|2])" << endl;
        exit(1);
    }

    // arguments that are needed
    power2   = atoi(argv[1]);
    filename = argv[2];

    // argument that have default values
    H         = nArgs > 3 ? atof(argv[3]) : 1.5;
    minZValue = nArgs > 4 ? atof(argv[4]) : 0.0;
    maxZValue = nArgs > 5 ? atof(argv[5]) : 1.0;
    PBC       = nArgs > 6 ? atoi(argv[6]) : true;
    idum      = nArgs > 7 ? atol(argv[7]) : 1;
    RNG       = nArgs > 8 ? atoi(argv[8]) : 2;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
