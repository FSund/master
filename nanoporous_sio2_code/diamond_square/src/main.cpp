#include <iostream>
#include <cstdlib>      // atoi, atof, atol
#include <armadillo>
#include "diamondSquare.h"

typedef unsigned int uint;

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

    DiamondSquare m(power2, idum, RNG, PBC);
    mat R = m.generate(H, corners);

    cout << "R = " << endl << R << endl;

//    R.save("rmat.dat", raw_ascii);

    return 0;
}

void createTetrahedra(const mat &terrainNodes) {
    // 6 unique tetrahedra

    // regular tetrahedron in top left corner
//    terrainNodes(x,y);
//    terrainNodes(x+1,y);
//    gridNodes(x,y);
//    gridNodes(x,y+1);
    vector<uint> tetrahedronNodes = {
        convert_2d_indices_to_linear_index(x,y,n),
        convert_2d_indices_to_linear_index(x+1,y,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y+1,n)
    };

    // regular tetrahedron in bottom right corner
//    terrainNodes(x+1,y);
//    terrainNodes(x+1,y+1);
//    gridNodes(x,y+1);
//    gridNodes(x+1,y+1);
    vector<uint> tetrahedronNodes = {
        convert_2d_indices_to_linear_index(x+1,y,n),
        convert_2d_indices_to_linear_index(x+1,y+1,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y+1,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,n)
    };

    // two irregular tetrahedron with three points in the grid, in bottom left corner
//    gridNodes(x,y);
//    gridNodes(x+1,y);
//    gridNodes(x+1,y+1);
//    terrainNodes(x+1,y);
    vector<uint> tetrahedronNodes = {
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,n),
        convert_2d_indices_to_linear_index(x+1,y)
    };

//    gridNodes(x+1,y);
//    gridNodes(x,y+1);
//    gridNodes(x+1,y+1);
//    terrainNodes(x+1,y);
    vector<uint> tetrahedronNodes = {
        nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y+1,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,n),
        convert_2d_indices_to_linear_index(x+1,y)
    };

    // two irregular tetrahedron with three points in the terrain, in top right corner
//    terrainNodes(x,y);
//    terrainNodes(x+1,y);
//    terrainNodes(x,y+1);
//    gridNodes(x,y+1);
    vector<uint> tetrahedronNodes = {
        convert_2d_indices_to_linear_index(x,y,n),
        convert_2d_indices_to_linear_index(x+1,y,n),
        convert_2d_indices_to_linear_index(x,y+1,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y+1,n)
    };

//    terrainNodes(x+1,y+1);
//    terrainNodes(x+1,y);
//    terrainNodes(x,y+1);
//    gridNodes(x,y+1);
    vector<uint> tetrahedronNodes = {
        convert_2d_indices_to_linear_index(x+1,y+1,n),
        convert_2d_indices_to_linear_index(x+1,y,n),
        convert_2d_indices_to_linear_index(x,y+1,n),
        nTerrainNodes + convert_2d_indices_to_linear_index(x,y+1,n)
    };
}

vector<vector<uint> > createElements(const mat &terrainNodes, const mat &gridNodes) {

}

inline int convert_2d_indices_to_linear_index(int subx, int suby, int &nx, int &ny) {
    /* converts 3d subscripts (indexes) to a linear index */
    return subx*ny + suby;
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
