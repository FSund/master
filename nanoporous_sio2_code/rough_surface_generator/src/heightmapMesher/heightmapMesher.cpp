#include <src/heightmapMesher/heightmapMesher.h>

void HeightmapMesher::mesh(const mat &heightmap, string filename) {

    if (heightmap.n_rows != heightmap.n_cols) {
        cout << "Error in HeightmapMesher::heighmapToEleNode, terrainNodes.n_rows != terrainNodes.n_cols" << endl;
        exit(1);
    }

    nTerrainNodes = heightmap.n_elem;
    gridSize = heightmap.n_rows;

    cout << "nTerrainNodes = " << nTerrainNodes << endl;

    vector<vector<uint> > tetrahedrons; // TODO: figure out how big this will be beforehand (should be easy), possibly convert to umat?

    // loop over inner points (not including the column at the left and right edge, and the two rows at the bottom edge)
    int a = 0;
    for (uint x = 0; x < gridSize-2; x++) { // not including the two bottom rows
        for (uint y = 1; y < gridSize-1; y++) { // not including first and last column
            cout << "default" << endl;
            for (uint tetrahedronType = 0; tetrahedronType < 6; tetrahedronType++) {
                tetrahedrons.push_back(makeTetrahedron(tetrahedronType, x, y, gridSize, nTerrainNodes));
                a++;
                cout << "a = " << a << ", x = " << x << ", y = " << y
                     << ", tetrahedron = " << tetrahedrons.back()[0] << "," << tetrahedrons.back()[1] << "," << tetrahedrons.back()[2] << "," << tetrahedrons.back()[3] << endl;
            }
        }
    }
    cout << "n tetrahedrons made = " << a << endl;

    // loop over left edge (not including endpoint (0,gridSize-1) or (0,gridSize-2))
    uint x;
    uint y = 0; // left edge
    for (x = 0; x < gridSize-2; x++) { // not including the two last points
        cout << "left edge" << endl;
        for (uint tetrahedronType = 0; tetrahedronType < 2; tetrahedronType++) {
            tetrahedrons.push_back(makeLeftEdgeTetrahedron(tetrahedronType, x, y, gridSize, nTerrainNodes));
        }
    }
    // manual irregular tetrahedron with (gridSize-1, 0) at top left corner of triangle
    x = gridSize-2; // next to the last point
    y = 0;
    uint tetrahedronType = 0;
    tetrahedrons.push_back(makeLeftEdgeTetrahedron(tetrahedronType, x, y, gridSize, nTerrainNodes));

    // loop over bottom edge (not including the startpoint (gridSize-1, 0) and the endpoint (gridSize-1,gridSize-1))
    x = gridSize-1; // bottom edge
    for (y = 1; y < gridSize-1; y++) { // not including the start and endpoints
        cout << "bottom edge" << endl;
        for (uint tetrahedronType = 0; tetrahedronType < 2; tetrahedronType++) {
            tetrahedrons.push_back(makeBottomEdgeTetrahedron(tetrahedronType, x, y, gridSize, nTerrainNodes));
        }
    }
    // manual irregular tetrahedron with (gridSize-1, gridSize-1) at bottom right corner of triangle
    x = gridSize-1;
    y = gridSize-1;
    tetrahedronType = 0;
    tetrahedrons.push_back(makeBottomEdgeTetrahedron(tetrahedronType, x, y, gridSize, nTerrainNodes));

    printToEleNode(heightmap, tetrahedrons, filename);
    printToMsh(heightmap, tetrahedrons, filename);
}

vector<uint> HeightmapMesher::makeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes) {

    // 6 unique tetrahedra
    switch (tetrahedronType) {
    case 0:
    {
        // regular tetrahedron in top left corner, with terrain node (x,y) at top and zero height node (x,y-1) at left
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y-1,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    case 1:
    {
        // regular tetrahedron in bottom right corner, with terrain node (x+1,y) at left and zero height node (x,y) at top
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    case 2:
    {
        // irregular tetrahedron with three points in the zero height grid, in bottom left corner, with zero heigh node (x,y-1) in top left corner of the triangle in the zero height grid
        vector<uint> tetrahedronNodes = {
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y-1,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y-1,gridSize-1,gridSize-1),
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    case 3:
    {
        // irregular tetrahedron with three points in the zero height grid, in bottom left corner, with zero height node (x+1,y) in bottom right corner of the triangle in the zero height grid
        vector<uint> tetrahedronNodes = {
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y-1,gridSize-1,gridSize-1),
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize)
        };
        return tetrahedronNodes;
    }
    case 4:
    {
        // irregular tetrahedron with three points in the terrain, in top right corner, with terrain node (x+1,y) in bottom left corner of the triangle in the terrain nodes
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    case 5:
    {
        // irregular tetrahedron with three points in the terrain, in top right corner, with terrain node (x,y+1) in top right corner of the triangle in the terrain nodes
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    default:
        cout << "Error in HeightmapMesher::makeTetrahedron" << endl;
        exit(1);
    }
}

vector<uint> HeightmapMesher::makeLeftEdgeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes) {
    // for left edge //
    switch (tetrahedronType) {
    case 0:
    {
        // irregular tetrahedron using terrain node x,y as top left corner of triangle
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    case 1:
    {
        // regular tetrahedron with terrain node (x+1,y) at left and zero node (x,y) at top
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x+1,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x+1,y+1,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x,y,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x+1,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    default:
        cout << "Error in HeightmapMesher::makeLeftEdgeTetrahedron" << endl;
        exit(1);
    }
}

vector<uint> HeightmapMesher::makeBottomEdgeTetrahedron(uint tetrahedronType, uint x, uint y, uint gridSize, uint nTerrainNodes) {
    // for bottom edge //
    switch (tetrahedronType) {
    case 0:
    {
        // irregular tetrahedron using x,y as bottom right corner of triangle
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x,y-1,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x-1,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x-1,y-1,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    case 1:
    {
        // regular tetrahedron with terrain node (x,y) at bottom and zero node (x,y) at the left
        vector<uint> tetrahedronNodes = {
            convert_2d_indices_to_linear_index(x,y,gridSize,gridSize),
            convert_2d_indices_to_linear_index(x-1,y,gridSize,gridSize),
            nTerrainNodes + convert_2d_indices_to_linear_index(x-1,y-1,gridSize-1,gridSize-1),
            nTerrainNodes + convert_2d_indices_to_linear_index(x-1,y,gridSize-1,gridSize-1)
        };
        return tetrahedronNodes;
    }
    default:
        cout << "Error in HeightmapMesher::makeBottomEdgeTetrahedron" << endl;
        exit(1);
    }
}

void HeightmapMesher::printToEleNode(const mat &heightmap, const vector<vector<uint> > &tetrahedrons, string filename) {
    // Print nodes (points) to file
    stringstream filenameMaker;
    filenameMaker << filename << ".node";

    ofstream ofile;
    ofile.open(filenameMaker.str().c_str());
    ofile << nTerrainNodes*2 << " 3 0 0" << endl;
    uint counter = 1; // counting starts at 1 in .node and .ele files
    double norm = 1.0/gridSize;

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

    // Printing the zero grid to file
    for (uint i = 0; i < gridSize; i++) {
        double x = i*norm;
        for (uint j = 0; j < gridSize; j++) {
            double y = j*norm;
            double z = 0.0;
            ofile << counter << " " << x << " " << y << " " << z << endl;
            counter++;
        }
    }
    ofile.close();

    filenameMaker.str(""); // clear the string stream
    filenameMaker << filename << ".ele";

    ofile.open(filenameMaker.str().c_str());
    ofile << tetrahedrons.size() << " 4 0" << endl;
    counter = 1;
    for (vector<vector<uint> >::const_iterator elements = tetrahedrons.begin(); elements != tetrahedrons.end(); ++elements) {
        ofile << counter << " ";
        for (vector<uint>::const_iterator element = elements->begin(); element != elements->end(); ++element) {
            ofile << (*element) + 1 << " "; // counting starts at 1 in .node and .ele files
        }
        ofile << endl;
        counter++;
    }
}

void HeightmapMesher::printToMsh(const mat& heightmap, const vector<vector<uint> >& tetrahedrons, string filename) {

    uint nNodes = nTerrainNodes + (gridSize-1)*(gridSize-1);

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
    // Printing the zero grid to file
    for (uint i = 0; i < gridSize-1; i++) {
        double x = (i + 0.5)*norm;
        for (uint j = 0; j < gridSize-1; j++) {
            double y = (j + 0.5)*norm;
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

