#ifndef HEIGHTMAPMESHER_H
#define HEIGHTMAPMESHER_H

#include <fstream>
#include <sstream>
#include <vector>
#include <armadillo>

typedef unsigned int uint;

using namespace std;
using namespace arma;

class HeightmapMesher{
public:
    void mesh(const mat &heightMap, string filename);

private:
    uint nTerrainNodes;
    uint gridSize;
    vector<uint> makeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes);
    vector<uint> makeLeftEdgeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes);
    vector<uint> makeBottomEdgeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes);

    void printToEleNode(const mat& heightmap, const vector<vector<uint> >& tetrahedrons, string filename);
    void printToMsh(const mat& heightmap, const vector<vector<uint> >& tetrahedrons, string filename);
};

inline uint convert_2d_indices_to_linear_index(int subx, int suby, int nx, int ny) {
    /* converts 3d subscripts (indexes) to a linear index */
    return subx*ny + suby;
}

#endif // HEIGHTMAPMESHER_H
