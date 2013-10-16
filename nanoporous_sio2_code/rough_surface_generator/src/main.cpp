#include <iostream>
#include <cstdlib>      // atoi, atof, atol
#include <armadillo>
#include <src/diamondSquare/diamondSquare.h>
#include <src/heightmapMesher/heightmapMesher.h>

using namespace std;
using namespace arma;

void parseArgs(int nArgs, const char* argv[], int &power2, string &filename, double &H, double &corners, bool &PBC, long &idum, int &RNG);

int main(int nArgs, const char *argv[])
{
    int power2;
    string filename;
    double H;
    double corners;
    bool PBC;
    long idum;
    int RNG;

    parseArgs(nArgs, argv, power2, filename, H, corners, PBC, idum, RNG);

    cout << "--- Diamond-square settings --------------------" << endl;
    cout << "power2 = " << power2  << endl;
    cout << "filename = " << filename << endl;
    cout << "H = " << H << endl;
    cout << "corners = " << corners << endl;
    cout << "PBC = " << std::boolalpha << PBC << std::noboolalpha << endl;
    cout << "idum = " << idum << endl;
    cout << "RNG = " << RNG << endl;
    cout << "total number of points in grid = " << pow(pow(2, power2)+1, 2) << endl;
    cout << "------------------------------------------------" << endl;

    srand(idum); // setting the seed of the RNG for both C++'s rand()/srand() and Armadillo's randu()/randn()

    DiamondSquare generator(power2, idum, RNG, PBC);
    mat R = generator.generate(H, corners);

//    cout << "R = " << endl << R << endl;

    HeightmapMesher mesher;

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
        string &filename,
        double &H,
        double &corners,
        bool &PBC,
        long &idum,
        int &RNG) {

    double default_H = 0.8;
    double default_corners = 0.5;
    bool default_PBC = true;
    long default_idum = 1;
    int default_RNG = 1;

    if (nArgs == 3) {
        power2   = atoi(argv[1]);
        filename = argv[2];
        H        = default_H;
        corners  = default_corners;
        PBC      = default_PBC;
        idum     = default_idum;
        RNG      = default_RNG;
    } else if (nArgs == 4) {
        power2   = atoi(argv[1]);
        filename = argv[2];
        H        = atof(argv[3]);
        corners  = default_corners;
        PBC      = default_PBC;
        idum     = default_idum;
        RNG      = default_RNG;
    } else if (nArgs == 5) {
        power2   = atoi(argv[1]);
        filename = argv[2];
        H        = atof(argv[3]);
        corners  = atof(argv[4]);
        PBC      = default_PBC;
        idum     = default_idum;
        RNG      = default_RNG;
    } else if (nArgs == 6) {
        power2   = atoi(argv[1]);
        filename = argv[2];
        H        = atof(argv[3]);
        corners  = atof(argv[4]);
        PBC      = atoi(argv[5]);
        idum     = default_idum;
        RNG      = default_RNG;
    } else if (nArgs == 7) {
        power2   = atoi(argv[1]);
        filename = argv[2];
        H        = atof(argv[3]);
        corners  = atof(argv[4]);
        PBC      = atoi(argv[5]);
        idum     = atol(argv[6]);
        RNG      = default_RNG;
    } else if (nArgs == 8) {
        power2   = atoi(argv[1]);
        filename = argv[2];
        H        = atof(argv[3]);
        corners  = atof(argv[4]);
        PBC      = atoi(argv[5]);
        idum     = atol(argv[6]);
        RNG      = atoi(argv[7]);
    } else {
        cout << "Usage: ./diamondSquare  power2  filename  optional:(H  corners  PBC[0|1]  idum[unsigned int]  RNG[0|1|2])" << endl;
        exit(1);
    }
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
