#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "bit_vector.h"

using namespace std;

typedef unsigned int uint;

vector<double> calculate_vector_difference_using_minimum_image_convention(
    const vector<double> &u, 
    const vector<double> &v, 
    const vector<double> &system_size, 
    const vector<double> &half_system_size){

    vector<double> dr_vector(3);
    dr_vector[0] = u[0] - v[0];
    dr_vector[1] = u[1] - v[1];
    dr_vector[2] = u[2] - v[2];

    // minimum image convention
    if      (dr_vector[0] >  half_system_size[0]) dr_vector[0] -= system_size[0];
    else if (dr_vector[0] < -half_system_size[0]) dr_vector[0] += system_size[0];
    if      (dr_vector[1] >  half_system_size[1]) dr_vector[1] -= system_size[1];
    else if (dr_vector[1] < -half_system_size[1]) dr_vector[1] += system_size[1];
    if      (dr_vector[2] >  half_system_size[2]) dr_vector[2] -= system_size[2];
    else if (dr_vector[2] < -half_system_size[2]) dr_vector[2] += system_size[2];

    return dr_vector;
}

double calculate_distance_squared_using_minimum_image_convention(
    const vector<double> &u, 
    const vector<double> &v, 
    const vector<double> &system_size, 
    const vector<double> &half_system_size) {

    vector<double> dr_vector = calculate_vector_difference_using_minimum_image_convention(u, v, system_size, half_system_size);
    return dr_vector[0]*dr_vector[0] + dr_vector[1]*dr_vector[1] + dr_vector[2]*dr_vector[2];
}

void print_bitmatrix_to_vtk(const BitMatrix& matrix, const vector<double> voxel_size_angstrom, string output_file) {
    cout << "Printing to .vtk-format." << endl;

    ofstream ofile;
    ofile.open(output_file);

    uint nx = matrix.n_rows;
    uint ny = matrix.n_cols;
    uint nz = matrix.n_slices;

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << "structured point" << endl;
    ofile << "ASCII" << endl;
    ofile << endl;
    ofile << "DATASET STRUCTURED_POINTS" << endl;
    ofile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    // ofile << "ORIGIN 0.0 0.0 0.0" << endl;
    ofile << "ORIGIN " << setprecision(16) << voxel_size_angstrom[0]/2.0 << " " << voxel_size_angstrom[1]/2.0 << " " << voxel_size_angstrom[2]/2.0 << endl;
    // ofile << "SPACING " << 1.0/double(nx) << " " << 1.0/double(ny) << " " << 1.0/double(nz) << endl;
    ofile << "SPACING " << setprecision(16) << voxel_size_angstrom[0] << " " << voxel_size_angstrom[0] << " " << voxel_size_angstrom[0] << endl;
    ofile << "POINT_DATA " << nx*ny*nz << endl;
    ofile << "SCALARS occupation double" << endl;
    ofile << "LOOKUP_TABLE default" << endl;
    ofile << endl;

    // column major ordering
    for (uint k = 0; k < nz; k++){
        for (uint j = 0; j < ny; j++){
            for (uint i = 0; i < nx; i++){
                ofile << matrix(i,j,k) << endl;
            }
        }
    }

    // // row major ordering
    // for (int i = 0; i < matrix_dim_x; i++) {
    //     for (int j = 0; j < matrix_dim_y; j++) {
    //         for (int k = 0; k < matrix_dim_z; k++) {
    //             ofile << matrix[i][j][k] << endl;
    //         }
    //     }
    // }

    ofile.close();

    cout << "Finished printing to .vtk-format." << endl;
}

vector<vector<uint> > create_radial_neighbor_list(Mts0_io* mts0_io, const double neighbor_list_radius) {
    /*
     *Strategy: create neighbor list with method from Mts0_io. The method from Mts0 says that all atoms in neighboring boxes are neighbors.
     *          This method runs through the neighbor list from Mts0_io and checks if the atoms really are neighbors by calculating
     *          the distance between them and checking if it is smaller than neighbor_list_radius.
     *          The indexes of all neighboring atoms that have a distance smaller than neighbor_list_radius are placed in radial_neighbor_list.
     */

    vector<vector<int> > neighbor_list = mts0_io->create_neighbor_list(neighbor_list_radius);
    vector<vector<uint> > radial_neighbor_list(neighbor_list.size());

    vector<vector<double> > &positions = mts0_io->positions;
    vector<double> system_size = mts0_io->get_lx_ly_lz();
    vector<double> half_system_size = mts0_io->get_lx_ly_lz();
    half_system_size[0] *= 0.5;
    half_system_size[1] *= 0.5;
    half_system_size[2] *= 0.5;
    double neighbor_list_radius_squared = neighbor_list_radius*neighbor_list_radius;

    // Run through neighbor_list, and add the atoms that are within neighbor_list_radius to radial_neighbor_list
    for (uint i = 0; i < neighbor_list.size(); i++) {
        uint number_of_neighbors = neighbor_list[i].size();
        vector<double> main_atom_position = positions[i];
        vector<uint> my_neighbors;
        for (uint j = 0; j < number_of_neighbors; j++) {
            vector<double> neighbor_position = positions[neighbor_list[i][j]];

            // Find the distance between the atoms using minimium image convention
            double distance_between_main_atom_neighbor_squared = calculate_distance_squared_using_minimum_image_convention(neighbor_position, main_atom_position, system_size, half_system_size);
            if (distance_between_main_atom_neighbor_squared < neighbor_list_radius_squared) {
                my_neighbors.push_back(neighbor_list[i][j]);
            }
        }
        radial_neighbor_list[i] = my_neighbors;
    }
    return radial_neighbor_list;
}
