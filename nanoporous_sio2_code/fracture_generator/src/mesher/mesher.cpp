#include <src/mesher/mesher.h>

void Mesher::mesh(const mat &topSurface, const mat &bottomSurface) {

    if (topSurface.n_rows != topSurface.n_cols || bottomSurface.n_rows != bottomSurface.n_cols || topSurface.n_rows != bottomSurface.n_rows || topSurface.n_cols != bottomSurface.n_cols) {
        cout << "Error in Mesher::mesh() !" << endl;
        exit(1);
    }

    gridSize = topSurface.n_rows;
    nGridPoints = topSurface.n_elem;

    // Populating the list of nodes/points
    nodes = mat(2*nGridPoints, 3);
    for (uint i = 0; i < gridSize; i++) {
        for (uint j = 0; j < gridSize; j++) {
            nodes(linearIndexInTopSurface(i,j), 0) = i;               // x-coordinate for top surface
            nodes(linearIndexInTopSurface(i,j), 1) = j;               // y-coordinate for top surface
            nodes(linearIndexInTopSurface(i,j), 2) = topSurface(i,j); // z-coordinate for top surface
        }
    }
    for (uint i = 0; i < gridSize; i++) {
        for (uint j = 0; j < gridSize; j++) {
            nodes(linearIndexInBottomSurface(i,j), 0) = i;                  // x-coordinate for bottom surface
            nodes(linearIndexInBottomSurface(i,j), 1) = j;                  // y-coordinate for bottom surface
            nodes(linearIndexInBottomSurface(i,j), 2) = bottomSurface(i,j); // z-coordinate for bottom surface
        }
    }

    // Loop over all all points except right and bottom edge
    for (uint x = 0; x < gridSize-1; x++) { // not including the bottom edge
        for (uint y = 0; y < gridSize-1; y++) { // not including right edge
            for (uint tetrahedronType = 0; tetrahedronType < 5; tetrahedronType++) {
                vector<uint> tetrahedron = makeTetrahedron(tetrahedronType, x, y);

                // Check if we have a "negative" tetrahedron, ie. the top surface is below the bottom surface in any of
                // the points this tetrahedron consists of
                bool includeThisTetrahedron = true;
                for (auto it = tetrahedron.begin(); it != tetrahedron.end(); ++it) {
                    if (nodes(linearIndexInTopSurface(x,y), 2) <= nodes(linearIndexInBottomSurface(x,y), 2)) {
                        includeThisTetrahedron = false;
                        break;
                    }
                }
                if (includeThisTetrahedron) {
                    tetrahedrons.push_back(tetrahedron);
                }
            }
        }
    }

//    printToMsh(heightmap, tetrahedrons, filename, baseMeshHeight);

    cout << "Mesher::mesh() made " << tetrahedrons.size() << " tetrahedrons." << endl;
}

vector<uint> Mesher::makeTetrahedron(const uint tetrahedronType, const uint x, const uint y) {

    // only need 5 tetrahedra to fill a cube
    switch (tetrahedronType) {
    case 0:
    {
        // regular tetrahedron in the center of the cube, using top nodes (x+1,y) and (x,y+1), and bottom nodes (x,y) and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            linearIndexInTopSurface(x+1,y),
            linearIndexInTopSurface(x,y+1),
            linearIndexInBottomSurface(x,y),
            linearIndexInBottomSurface(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    case 1:
    {
        // irregular tetrahedron on top left corner, with triangle in top nodes. Using top nodes (x,y), (x+1,y) and (x,y+1), and bottom node (x,y)
        vector<uint> tetrahedronNodes = {
            linearIndexInTopSurface(x,y),
            linearIndexInTopSurface(x+1,y),
            linearIndexInTopSurface(x,y+1),
            linearIndexInBottomSurface(x,y)
        };
        return tetrahedronNodes;
    }
    case 2:
    {
        // irregular tetrahedron on bottom left corner, with triangle in bottom nodes. Using top node (x+1,y), and bottom nodes (x+1,y), (x,y), and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            linearIndexInTopSurface(x+1,y),
            linearIndexInBottomSurface(x+1,y),
            linearIndexInBottomSurface(x,y),
            linearIndexInBottomSurface(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    case 3:
    {
        // irregular tetrahedron on bottom right corner, with triangle in top nodes. Using top nodes (x+1,y+1), (x+1,y) and (x,y+1), and bottom node (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            linearIndexInTopSurface(x+1,y+1),
            linearIndexInTopSurface(x+1,y),
            linearIndexInTopSurface(x,y+1),
            linearIndexInBottomSurface(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    case 4:
    {
        // irregular tetrahedron on top right corner, with triangle in bottom nodes. Using top node (x,y+1), and bottom nodes (x,y+1), (x,y), and (x+1,y+1)
        vector<uint> tetrahedronNodes = {
            linearIndexInTopSurface(x,y+1),
            linearIndexInBottomSurface(x,y+1),
            linearIndexInBottomSurface(x,y),
            linearIndexInBottomSurface(x+1,y+1)
        };
        return tetrahedronNodes;
    }
    default:
        cout << "Error in HeightmapMesher::makeTetrahedron" << endl;
        exit(1);
    }
}

inline uint Mesher::linearIndexInTopSurface(const uint x, const uint y) {
    return x*gridSize + y;
}

inline uint Mesher::linearIndexInBottomSurface(const uint x, const uint y) {
    return nGridPoints + linearIndexInTopSurface(x,y);
}

void Mesher::printToMsh(std::string filename) {

//    uint nNodes = 2*nTerrainNodes;

    ofstream ofile;
    ofile.open(filename);

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

    // Printing nodes to file
    double gridNorm = 1.0/(gridSize-1); // scaling the grid so x and y is in [0,1]
    for (uint i = 0; i < 2*nGridPoints; i++) {
        ofile << i+1 << " ";                 // node number (numbering starts at 1 in .msh files)
        ofile << nodes(i,0)*gridNorm << " "; // x-coordinate
        ofile << nodes(i,1)*gridNorm << " "; // y-coordinate
        ofile << nodes(i,2);                 // z-coordinate
        ofile << endl;
    }
    ofile << "$EndNodes" << endl;

    // Print elements (tetrahedrons) to file //
        // formatting:
        // number-of-elements
        // element-number  element-type(4 for 4-node tetrahedron)  number-of-tags(set to 2)  tags: "0 0"  node-numbers(4 of them for tetrahedron)
        // example: 4588 4 2 0 30 96 463 97 1677
    ofile << "$Elements" << endl;
    ofile << tetrahedrons.size() << endl;

    uint counter = 1;
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

void Mesher::printToPLY(std::string filename) {

    ofstream ofile;
    ofile.open(filename);

    ofile << "ply" << endl;
    ofile << "format ascii 1.0" << endl;
    ofile << "comment made by Filip" << endl;

    // Define "vertex" element, and how many
    ofile << "element vertex " << 2*nGridPoints << endl;
    ofile << "property float x" << endl;
    ofile << "property float y" << endl;
    ofile << "property float z" << endl;

    // Define "face" element, and how many
    ofile << "element face " << 4*tetrahedrons.size() << endl;
    ofile << "property list uchar int vertex_index" << endl;

    ofile << "end_header" << endl;

    // Print nodes/vertices/points to file //
        // formatting:
        // x1  y1  z1
    for (uint i = 0; i < 2*nGridPoints; i++) {
        ofile << nodes(i,0) << " "; // x-coordinate
        ofile << nodes(i,1) << " "; // y-coordinate
        ofile << nodes(i,2);        // z-coordinate
        ofile << endl;
    }

    // Print faces to file //
        // formatting:
        // number_of_vertices vertex0 vertex1 vertex2
    vector<vector<uint> > face_indices(4);
    face_indices[0] = {0, 1, 2};
    face_indices[1] = {1, 2, 3};
    face_indices[2] = {0, 2, 3};
    face_indices[3] = {0, 1, 3};

    for (uint i = 0; i < tetrahedrons.size(); i++) {
        for (uint j = 0; j < face_indices.size(); j++) {
            ofile << "3 "; // number of vertices
            for (auto k = face_indices[j].begin(); k != face_indices[j].end(); ++k) {
                ofile << tetrahedrons[i][*k] << " ";
            }
            ofile << endl;
        }

//        for (auto element = elements->begin(); element != elements->end(); ++element) {
//            ofile << (*element) << " "; // vertex number
//        }
//        ofile << endl;
    }

    ofile.close();
}

//    ply
//    format ascii 1.0           { ascii/binary, format version number }
//    comment made by Greg Turk  { comments keyword specified, like all lines }
//    comment this file is a cube
//    element vertex 8           { define "vertex" element, 8 of them in file }
//    property float x           { vertex contains float "x" coordinate }
//    property float y           { y coordinate is also a vertex property }
//    property float z           { z coordinate, too }
//    element face 6             { there are 6 "face" elements in the file }
//    property list uchar int vertex_index { "vertex_indices" is a list of ints }
//    end_header                 { delimits the end of the header }
//    0 0 0                      { start of vertex list }
//    0 0 1
//    0 1 1
//    0 1 0
//    1 0 0
//    1 0 1
//    1 1 1
//    1 1 0
//    4 0 1 2 3                  { start of face list }
//    4 7 6 5 4
//    4 0 4 5 1
//    4 1 5 6 2
//    4 2 6 7 3
//    4 3 7 4 0

//void Mesher::printToSTL(std::string filename) {

////    solid name
////        facet normal ni nj nk
////            outer loop
////                vertex v1x v1y v1z
////                vertex v2x v2y v2z
////                vertex v3x v3y v3z
////            endloop
////        endfacet
////    endsolid name

//    ofstream ofile;
//    ofile.open(filename);

//    ofile << "$MeshFormat" << endl;
//    ofile << "2.2 0 8" << endl;
//    ofile << "$EndMeshFormat" << endl;
//    ofile << "$Comments" << endl;
//    // optional comments
//    ofile << "$EndComments" << endl;

//    // Print nodes (points) to file //
//        // formatting:
//        // number-of-nodes
//        // node-number  x1  y1  z1
//    ofile << "$Nodes" << endl;
//    ofile << 2*nGridPoints << endl;

//    // Printing nodes to file
//    double gridNorm = 1.0/(gridSize-1); // scaling the grid so x and y is in [0,1]
//    for (uint i = 0; i < 2*nGridPoints; i++) {
//        ofile << i+1 << " ";                 // node number (numbering starts at 1 in .msh files)
//        ofile << nodes(i,0)*gridNorm << " "; // x-coordinate
//        ofile << nodes(i,1)*gridNorm << " "; // y-coordinate
//        ofile << nodes(i,2);                 // z-coordinate
//        ofile << endl;
//    }
//    ofile << "$EndNodes" << endl;

//    // Print elements (tetrahedrons) to file //
//        // formatting:
//        // number-of-elements
//        // element-number  element-type(4 for 4-node tetrahedron)  number-of-tags(set to 2)  tags: "0 0"  node-numbers(4 of them for tetrahedron)
//        // example: 4588 4 2 0 30 96 463 97 1677
//    ofile << "$Elements" << endl;
//    ofile << tetrahedrons.size() << endl;

//    uint counter = 1;
//    for (auto elements = tetrahedrons.begin(); elements != tetrahedrons.end(); ++elements) {
//        ofile << counter << " 4 2 0 1 ";
//        for (auto element = elements->begin(); element != elements->end(); ++element) {
//            ofile << (*element) + 1 << " "; // counting starts at 1 in .msh files
//        }
//        ofile << endl;
//        counter++;
//    }

//    ofile << "$EndElements" << endl;

//    ofile.close();
//}
