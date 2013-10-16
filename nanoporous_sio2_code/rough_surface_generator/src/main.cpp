#include <iostream>
#include <cstdlib>      // atoi, atof, atol
#include <armadillo>
#include <src/diamondSquare/diamondSquare.h>
#include <src/heightmapMesher/heightmapMesher.h>

using namespace std;
using namespace arma;

void parseArgs(int nArgs, const char* argv[], int &power2, double &H, double &corners, bool &PBC, long &idum, int &RNG);

int main(int nArgs, const char *argv[])
{
    int power2;
    double H;
    double corners;
    bool PBC;
    long idum;
    int RNG;

    parseArgs(nArgs, argv, power2, H, corners, PBC, idum, RNG);

    cout << "--- Diamond-square settings --------------------" << endl;
    cout << "power2 = " << power2  << endl;
    cout << "H = " << H << endl;
    cout << "corners = " << corners << endl;
    cout << "PBC = " << std::boolalpha << PBC << std::noboolalpha << endl;
    cout << "idum = " << idum << endl;
    cout << "RNG = " << RNG << endl;
    cout << "------------------------------------------------" << endl;

    srand(idum); // setting the seed of the RNG for both C++'s rand()/srand() and Armadillo's randu()/randn()

    DiamondSquare generator(power2, idum, RNG, PBC);
    mat R = generator.generate(H, corners);

//    cout << "R = " << endl << R << endl;

    HeightmapMesher mesher;
    string filename = "test";

    mesher.mesh(R, filename);

//    R.save("rmat.dat", raw_ascii);

    return 0;
}



//inline void convert_linear_index_to_3d_indices(int &index, int &nx, int &ny, int &nz, int* subscript) {
//    /* converts linear index to 3d subscripts (indexes) */

//    subscript[0] = index/(ny*nz);   // Node id in x-direction
//    subscript[1] = (index/nz)%ny;   // Node id in y-direction
//    subscript[2] = index%nz;        // Node id in z-direction
//}

void parseArgs(
        int nArgs,
        const char *argv[],
        int &power2,
        double &H,
        double &corners,
        bool &PBC,
        long &idum,
        int &RNG) {

    if (nArgs == 2) {
        power2  = atoi(argv[1]);
        H       = 0.8;
        corners = 0.5;
        PBC     = 1;
        idum    = 1;
        RNG     = 1;

//        cout << "Using default values" << endl;
//        cout << "H       = " << H << endl;
//        cout << "corners = " << corners << endl;
//        cout << "PBC     = " << PBC << endl;
//        cout << "idum    = " << idum << endl;
//        cout << "RNG     = " << RNG << endl;
//        cout << endl;
    } else if (nArgs == 3) {
        power2  = atoi(argv[1]);
        H       = atof(argv[2]);
        corners = 0.5;
        PBC     = 1;
        idum    = 1;
        RNG     = 1;

//        cout << "Using default values" << endl;
//        cout << "corners = " << corners << endl;
//        cout << "PBC     = " << PBC << endl;
//        cout << "idum    = " << idum << endl;
//        cout << "RNG     = " << RNG << endl;
//        cout << endl;
    } else if (nArgs == 4) {
        power2  = atoi(argv[1]);
        H       = atof(argv[2]);
        corners = atof(argv[3]);
        PBC     = 1;
        idum    = 1;
        RNG     = 1;

//        cout << "Using default values" << endl;
//        cout << "PBC  = " << PBC << endl;
//        cout << "idum = " << idum << endl;
//        cout << "RNG  = " << RNG << endl;
//        cout << endl;
    } else if (nArgs == 5) {
        power2  = atoi(argv[1]);
        H       = atof(argv[2]);
        corners = atof(argv[3]);
        PBC     = atoi(argv[4]);
        idum    = 1;
        RNG     = 1;

//        cout << "Using default values" << endl;
//        cout << "idum = " << idum << endl;
//        cout << "RNG  = " << RNG << endl;
//        cout << endl;
    } else if (nArgs == 6) {
        power2  = atoi(argv[1]);
        H       = atof(argv[2]);
        corners = atof(argv[3]);
        PBC     = atoi(argv[4]);
        idum    = atol(argv[5]);
        RNG     = 1;

//        cout << "Using default values" << endl;
//        cout << "RNG  = " << RNG << endl;
    } else if (nArgs == 7) {
        power2  = atoi(argv[1]);
        H       = atof(argv[2]);
        corners = atof(argv[3]);
        PBC     = atoi(argv[4]);
        idum    = atol(argv[5]);
        RNG     = atoi(argv[6]);
    } else {
//        power2  = 3;
//        H       = 0.8;
//        corners = 0.5;
//        PBC     = 1;
//        idum    = 1;
//        RNG     = 1;

        cout << "Usage: ./diamondSquare  power2  optional:(H  corners  PBC[0|1]  idum[unsigned int]  RNG[0|1|2])" << endl;
        exit(1);
//        cout << "Using default values" << endl;
//        cout << "power2  = " << power2 << endl;
//        cout << "H       = " << H << endl;
//        cout << "corners = " << corners << endl;
//        cout << "PBC     = " << PBC << endl;
//        cout << "idum    = " << idum << endl;
//        cout << "RNG     = " << RNG << endl;
//        cout << endl;
    }
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
