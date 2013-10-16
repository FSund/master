#include <src/heightmapMesher/heightmapMesher.h>

void HeightmapMesher::mesh(const mat &heightmap, string filename) {

    if (heightmap.n_rows != heightmap.n_cols) {
        cout << "Error in HeightmapMesher::heighmapToEleNode, terrainNodes.n_rows != terrainNodes.n_cols" << endl;
        exit(1);
    }

    nTerrainNodes = heightmap.n_elem;
    gridSize = heightmap.n_rows;

//    cout << "nTerrainNodes = " << nTerrainNodes << endl;

    vector<vector<uint> > tetrahedrons; // TODO: figure out how big this will be beforehand (should be easy), possibly convert to umat?

    // loop over all all points except right and bottom edge
    for (uint x = 0; x < gridSize-1; x++) { // not including the bottom edge
        for (uint y = 0; y < gridSize-1; y++) { // not including right edge
            for (uint tetrahedronType = 0; tetrahedronType < 5; tetrahedronType++) {
                tetrahedrons.push_back(makeTetrahedron(tetrahedronType, x, y, gridSize, nTerrainNodes));
            }
        }
    }

//    printToEleNode(heightmap, tetrahedrons, filename);
    printToMsh(heightmap, tetrahedrons, filename);
}

vector<uint> HeightmapMesher::makeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes) {

    // only need 5 tetrahedra to fill a cube
    switch (tetrahedronType) {
    case 0:
    {
        // regular tetrahedron in the center of the cube, using terrain nodes (x+1,y) and (x,y+1), and zero height nodes (x,y) and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    case 1:
    {
        // irregular tetrahedron on top left corner, with triangle in terrain nodes. Using terrain nodes (x,y), (x+1,y) and (x,y+1), and zero heigh node (x,y)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    case 2:
    {
        // irregular tetrahedron on bottom left corner, with triangle in zero height nodes. Using terrain node (x+1,y), and zero heigh nodes (x+1,y), (x,y), and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    case 3:
    {
        // irregular tetrahedron on bottom right corner, with triangle in terrain nodes. Using terrain nodes (x+1,y+1), (x+1,y) and (x,y+1), and zero heigh node (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    case 4:
    {
        // irregular tetrahedron on top right corner, with triangle in zero height nodes. Using terrain node (x,y+1), and zero heigh nodes (x,y+1), (x,y), and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    default:
        cout << "Error in HeightmapMesher::makeTetrahedron" << endl;
        exit(1);
    }
}

void HeightmapMesher::printToMsh(const mat& heightmap, const vector<vector<uint> >& tetrahedrons, string filename) {

    uint nNodes = 2*nTerrainNodes;

    stringstream filenameMaker;
    filenameMaker << filename << ".msh";

    ofstream ofile;
    ofile.open(filenameMaker.str().c_str());
//    double norm = 1.0/gridSize;
    double norm = 1.0;

    ofile << "$MeshFormat" << endl;
    ofile << "2.2 0 8" << endl;
    ofile << "$EndMeshFormat" << endl;
    ofile << "$Comments" << endl;
    // optional comments
    ofile << "$EndComments" << endl;

    // Print nodes (points) to file //
        // formatting:
        // number-of-nodes
        // node-number  x1  y1  z1
    ofile << "$Nodes" << endl;
    ofile << nNodes << endl;

    uint counter = 1; // counting starts at 1 in .msh files
    // Printing the heightmap to file
    for (uint i = 0; i < gridSize; i++) {
        double x = i*norm;
        for (uint j = 0; j < gridSize; j++) {
            double y = j*norm;
            double z = heightmap(i,j);
            ofile << counter << " " << x << " " << y << " " << z << endl;
            counter++;
        }
    }
    // Printing the zero height grid to file
    for (uint i = 0; i < gridSize; i++) {
        double x = i*norm;
        for (uint j = 0; j < gridSize; j++) {
            double y = j*norm;
            double z = 0.0;
            ofile << counter << " " << x << " " << y << " " << z << endl;
            counter++;
        }
    }
    ofile << "$EndNodes" << endl;

    // Print elements (tetrahedrons) to file //
        // formatting:
        // number-of-elements
        // element-number  element-type(4 for 4-node tetrahedron)  number-of-tags(set to 2)  tags: "0 0"  node-numbers(4 of them for tetrahedron)
        // example: 4588 4 2 0 30 96 463 97 1677
    ofile << "$Elements" << endl;
    ofile << tetrahedrons.size() << endl;

    counter = 1;
    for (vector<vector<uint> >::const_iterator elements = tetrahedrons.begin(); elements != tetrahedrons.end(); ++elements) {
        ofile << counter << " 4 2 0 1 ";
        for (vector<uint>::const_iterator element = elements->begin(); element != elements->end(); ++element) {
            ofile << (*element) + 1 << " "; // counting starts at 1 in .msh files
        }
        ofile << endl;
        counter++;
    }

    ofile << "$EndElements" << endl;

    ofile.close();
}
