#include <src/mesher/mesher.h>

void Mesher::mesh(const mat &topSurface, const mat &bottomSurface) {

    if (topSurface.n_rows != topSurface.n_cols || bottomSurface.n_rows != bottomSurface.n_cols || topSurface.n_rows != bottomSurface.n_rows || topSurface.n_cols != bottomSurface.n_cols) {
        cout << "Error in Mesher::mesh() !" << endl;
        exit(1);
    }

    gridSize = topSurface.n_rows;
    nGridPoints = topSurface.n_elem;

    // we use 5 tetrahedrons to make up each "cube"
    uint nTetrahedrons = (gridSize-1)*(gridSize-1)*5;
    tetrahedrons = vector<vector<uint> >(nTetrahedrons, vector<uint>(4));

    // loop over all all points except right and bottom edge
    uint counter = 0;
    for (uint x = 0; x < gridSize-1; x++) { // not including the bottom edge
        for (uint y = 0; y < gridSize-1; y++) { // not including right edge
            for (uint tetrahedronType = 0; tetrahedronType < 5; tetrahedronType++) {
                tetrahedrons[counter] = makeTetrahedron(tetrahedronType, x, y);
                counter++;
            }
        }
    }

//    printToMsh(heightmap, tetrahedrons, filename, baseMeshHeight);

    cout << "Mesher::mesh() made " << nTetrahedrons << " tetrahedrons." << endl;
}

vector<uint> Mesher::makeTetrahedron(const uint tetrahedronType, const uint x, const uint y) {

    // only need 5 tetrahedra to fill a cube
    switch (tetrahedronType) {
    case 0:
    {
        // regular tetrahedron in the center of the cube, using top nodes (x+1,y) and (x,y+1), and bottom nodes (x,y) and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y),
            convert_2d_indices_to_linear_index(x,y+1),
            nGridPoints + convert_2d_indices_to_linear_index(x,y),
            nGridPoints + convert_2d_indices_to_linear_index(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    case 1:
    {
        // irregular tetrahedron on top left corner, with triangle in top nodes. Using top nodes (x,y), (x+1,y) and (x,y+1), and bottom node (x,y)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y),
            convert_2d_indices_to_linear_index(x+1,y),
            convert_2d_indices_to_linear_index(x,y+1),
            nGridPoints + convert_2d_indices_to_linear_index(x,y)
        };
        return tetrahedronNodes;
    }
    case 2:
    {
        // irregular tetrahedron on bottom left corner, with triangle in bottom nodes. Using top node (x+1,y), and bottom nodes (x+1,y), (x,y), and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y),
            nGridPoints + convert_2d_indices_to_linear_index(x+1,y),
            nGridPoints + convert_2d_indices_to_linear_index(x,y),
            nGridPoints + convert_2d_indices_to_linear_index(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    case 3:
    {
        // irregular tetrahedron on bottom right corner, with triangle in top nodes. Using top nodes (x+1,y+1), (x+1,y) and (x,y+1), and bottom node (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y+1),
            convert_2d_indices_to_linear_index(x+1,y),
            convert_2d_indices_to_linear_index(x,y+1),
            nGridPoints + convert_2d_indices_to_linear_index(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    case 4:
    {
        // irregular tetrahedron on top right corner, with triangle in bottom nodes. Using top node (x,y+1), and bottom nodes (x,y+1), (x,y), and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y+1),
            nGridPoints + convert_2d_indices_to_linear_index(x,y+1),
            nGridPoints + convert_2d_indices_to_linear_index(x,y),
            nGridPoints + convert_2d_indices_to_linear_index(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    default:
        cout << "Error in HeightmapMesher::makeTetrahedron" << endl;
        exit(1);
    }
}

inline uint Mesher::convert_2d_indices_to_linear_index(const uint x, const uint y) {
    /* converts 2d indexes to linear index */
    return x*gridSize + y;
}


void Mesher::printToMsh(const mat &topSurface, const mat &bottomSurface, std::string filename) {

//    uint nNodes = 2*nTerrainNodes;

    stringstream filenameMaker;
    filenameMaker << filename << ".msh";

    ofstream ofile;
    ofile.open(filenameMaker.str().c_str());

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
    ofile << 2*nGridPoints << endl;


    double gridNorm = 1.0/(gridSize-1); // scaling the grid so x and y is in [0,1]
    uint counter = 1; // counting starts at 1 in .msh files

    // Printing the the top surface to file
    for (uint i = 0; i < gridSize; i++) {
        double x = i*gridNorm;
        for (uint j = 0; j < gridSize; j++) {
            double y = j*gridNorm;
            double z = topSurface(i,j);
            ofile << counter << " " << x << " " << y << " " << z << endl;
            counter++;

        }
    }
    // Printing the the bottom surface to file
    for (uint i = 0; i < gridSize; i++) {
        double x = i*gridNorm;
        for (uint j = 0; j < gridSize; j++) {
            double y = j*gridNorm;
            double z = bottomSurface(i,j);
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
    for (auto elements = tetrahedrons.begin(); elements != tetrahedrons.end(); ++elements) {
        ofile << counter << " 4 2 0 1 ";
        for (auto element = elements->begin(); element != elements->end(); ++element) {
            ofile << (*element) + 1 << " "; // counting starts at 1 in .msh files
        }
        ofile << endl;
        counter++;
    }

    ofile << "$EndElements" << endl;

    ofile.close();
}
