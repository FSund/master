// #pragma once // enabled by default -- and maybe not necessary in main file, only in headers/includes?
#include <iostream>
#include <iomanip>
// #include <bitset>       // std::bitset
// #include <cmath>        // ceil, floor
#include <sstream>
#include <fstream>
#include <algorithm>        // std::find
#include <cmath>            // ceil, floor

#include <mts0_io.h>
// #include "matrix.h"
#include "lib.h"
#include "bit_vector.h"
// #include "voro++.hh"

using namespace std;
// using namespace voro;

typedef unsigned int uint;

#define SI_TYPE 1
#define A_TYPE 2
#define H_TYPE 3
#define O_TYPE 4
#define NA_TYPE 5
#define CL_TYPE 6
#define X_TYPE 7

// inline bool check_if_atom_of_type_exists_within_distance_in_neighbor_voxels(
//     const Mts0_io* mts0_io, 
//     const vector<vector<vector<vector<uint> > > >& atoms_in_voxels, 
//     const uint atom_index, 
//     const vector<uint>& n_voxels,
//     const uint i,
//     const uint j,
//     const uint k, 
//     const double distance_squared,
//     const vector<double>& system_size_angstrom,
//     const vector<double>& half_system_size_angstrom);

int main(int n_args, char* arg_vec[]) {
    
    if (n_args < 9) {
        // cout << "Arguments: input_mts0_folder  nx  ny  nz  n_voxels_x  n_voxels_y  n_voxels_z  atom_radius" << endl;
        cout << "Arguments: input_mts0_folder  nx  ny  nz  r_start  r_step  n_steps  output_file" << endl;
        return 1;
    }

    string input_mts0_folder = arg_vec[1];
    int nx = atoi(arg_vec[2]);
    int ny = atoi(arg_vec[3]);
    int nz = atoi(arg_vec[4]);
    double r_start = atof(arg_vec[5]);
    double r_step = atof(arg_vec[6]);
    double N = atoi(arg_vec[7]);
    string output_file = arg_vec[8];
    // uint n_voxels_x = atoi(arg_vec[5]);
    // uint n_voxels_y = atoi(arg_vec[6]);
    // uint n_voxels_z = atoi(arg_vec[7]);
    // double distance = atof(arg_vec[5]);

    cout << "Loading data from mts0-files in folder '" << input_mts0_folder << "' ...";
    Mts0_io mts0_io(nx, ny, nz);
    mts0_io.load_atoms(input_mts0_folder);
    cout << "DONE" << endl;

    // Remove hydrogen atoms
    // This might be a bad idea, since some of the hydrogen atoms are there to passivate Si/N, and is really part of the
    // Si structure. BUT we do it for now anyways.
    for (uint i = 0; i < mts0_io.positions.size(); i++) {
        if (mts0_io.atom_types[i] == H_TYPE) {
            mts0_io.remove_atom(i);
        }
    }
    mts0_io.rearrange_vectors_by_moving_atoms_from_end_to_locations_where_atoms_have_been_removed(); // To update vectors

    int pore_atom_type = O_TYPE;
    int matrix_atom_type = SI_TYPE;
    int tag = 10;

    // If we don't start at 0.0 we need to tag the atoms from 0.0 to r_start with a different tag
    if (r_start > 0.0) {
        tag_atom_type_within_distance_from_other_atom_type(&mts0_io, matrix_atom_type, pore_atom_type, r_start, tag);
        tag++;
    }

    vector<vector<double> > results;
    for (uint i = 1; i <= N; i++) {
        double distance = r_start + r_step*i;
        uint n_tagged_atoms = tag_atom_type_within_distance_from_other_atom_type(&mts0_io, matrix_atom_type, pore_atom_type, distance, tag);

        double number_density = 0.0;
        if (n_tagged_atoms > 0) {
            number_density = find_number_density_of_atom_type(&mts0_io, tag);
        }

        cout << endl;
        cout << "Distance from matrix = [" << r_start + r_step*(i-1) << ", " << distance << ")" << endl;
        cout << "Number of tagged atoms = " << n_tagged_atoms << endl;
        cout << "Number density [atoms/Angstrom^3] = " << number_density << endl;

        double M = 0.0180158;       // kg/mol
        double NA = 6.0221413e+23;  // Avogadro's number
        double density_water = number_density*M/NA*(1.0e30);
        cout << "Density [kg/m^3] = " << density_water << endl;

        vector<double> temp_results = {distance, double(n_tagged_atoms), number_density, density_water};
        results.push_back(temp_results);
        
        tag++;
    }

    ofstream fout(output_file);
    fout << "# start = " << r_start << ", step = " << r_step << ", n_steps = " << N << ", removing hydrogen atoms, input mts0 folder \"" << input_mts0_folder << "\"" << endl;
    fout << "# r  n_atoms  number_density[atoms/Angstrom^3]  density[kg/m^3]" << endl;
    for (uint i = 0; i < results.size(); i++) {
        for (uint j = 0; j < results[i].size(); j++) {
            fout << results[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();

    // mts0_io.write_to_lammps("test_tagged.lmp");

    return 0;
}

// uint tag_atom_type_within_distance_from_other_atom_type(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag) {

//     vector<vector<uint> > radial_neighbor_list = create_radial_neighbor_list(mts0_io, distance);

//     uint n_tagged_atoms = 0;
//     vector<double> system_size = mts0_io->get_lx_ly_lz();
//     vector<double> half_system_size = {system_size[0]/2.0, system_size[1]/2.0, system_size[2]/2.0};
//     for (int atom_index = 0; atom_index < mts0_io->get_number_of_atoms(); atom_index++) {
//         if (mts0_io->atom_types[atom_index] == atom_type_to_tag) {
//             for (auto it = radial_neighbor_list[atom_index].begin(); it != radial_neighbor_list[atom_index].end(); ++it) {
//                 if (mts0_io->atom_types[*it] == atom_type_to_check_against) {
//                     mts0_io->atom_types[atom_index] = tag;
//                     n_tagged_atoms++;
//                     continue;
//                 }
//             }
//         }
//     }

//     // if (n_tagged_atoms == 0) {
//     //     cout << "Warning: tag_atom_type_within_distance_from_other_atom_type tagged 0 atoms." << endl;
//     // }

//     return n_tagged_atoms;
// }
