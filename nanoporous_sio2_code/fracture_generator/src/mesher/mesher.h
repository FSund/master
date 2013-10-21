#ifndef MESHER_H
#define MESHER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <armadillo>

using namespace std;
using namespace arma;

typedef unsigned int uint;

class Mesher {
public:
    void mesh(const mat &topSurface, const mat &bottomSurface);
    const vector<vector<uint> > &getTetrahedrons() const {
        return tetrahedrons;
    }
    void printToMsh(const mat &topSurface, const mat &bottomSurface, string filename);
private:
    uint gridSize;
    uint nGridPoints;
    vector<vector<uint> > tetrahedrons;
    vector<uint> makeTetrahedron(const uint tetrahedronType, const uint x, const uint y);

    inline uint convert_2d_indices_to_linear_index(const uint x, const uint y);
};

#endif // MESHER_H
