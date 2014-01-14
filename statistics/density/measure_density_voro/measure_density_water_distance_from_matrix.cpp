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

// uint tag_atom_type_within_distance_from_other_atom_type(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag);
uint tag_atom_type_within_distance_from_other_atom_type(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag);
inline bool check_if_atom_of_type_exists_within_distance_in_neighbor_voxels(
    const Mts0_io* mts0_io, 
    const vector<vector<vector<vector<uint> > > >& atoms_in_voxels, 
    const uint atom_index, 
    const vector<uint>& n_voxels,
    const uint i,
    const uint j,
    const uint k, 
    const double distance_squared,
    const vector<double>& system_size_angstrom,
    const vector<double>& half_system_size_angstrom);

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
    int tag = X_TYPE;

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

uint tag_atom_type_within_distance_from_other_atom_type(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag) {

    vector<double> system_size_angstrom = mts0_io->get_lx_ly_lz();
    vector<double> half_system_size_angstrom = {
        system_size_angstrom[0]/2.0,
        system_size_angstrom[1]/2.0,
        system_size_angstrom[2]/2.0,
    };

    vector<uint> n_voxels = {
        uint(ceil(system_size_angstrom[0]/distance)),
        uint(ceil(system_size_angstrom[1]/distance)),
        uint(ceil(system_size_angstrom[2]/distance))
    };

    // If we use a very small distance (which we usually only do when starting at 0.0 and using a small step in this program)
    // we sometimes get very small voxels in the first couple of cycles. This uses a lot of memory, without any real need
    // for it, and the program takes a lot longer to finish. So we limit the max number of voxels.
    uint max_n_voxels = 256; // 512 uses ~3.1 GiB memory, 256 uses ~400 MiB (and 256 is a lot faster than 512)
    while (n_voxels[0] > max_n_voxels || n_voxels[1] > max_n_voxels || n_voxels[2] > max_n_voxels) {
        for (uint i = 0; i < 3; i++) {
            n_voxels[i] = ceil(double(n_voxels[i])/2.0);
        }
    }
    // cout << "n_voxels = " << n_voxels[0] << ", " << n_voxels[1] << ", " << n_voxels[2] << endl; 

    vector<double> voxel_size_angstrom = {
        system_size_angstrom[0]/n_voxels[0],
        system_size_angstrom[1]/n_voxels[1],
        system_size_angstrom[2]/n_voxels[2]
    };

    vector<vector<vector<vector<uint> > > > atoms_in_voxels(n_voxels[0], 
        vector<vector<vector<uint> > >(n_voxels[1], 
            vector<vector<uint> >(n_voxels[2], 
                vector<uint>()
            )
        )
    );

    // // DEBUG //
    // cout << "hei" << endl;

    // // vector<vector<uint> > testvec;
    // vector<vector<uint> > testvec(1, vector<uint>());
    // cout << "max size = " << testvec.max_size() << endl;
    // cout << "wanted size = " << n_voxels[0]*n_voxels[1]*n_voxels[2] << endl;

    // testvec.resize(1);
    // cout << "max size = " << testvec.max_size() << endl;

    // testvec.resize(n_voxels[0]*n_voxels[1]*n_voxels[2]);
    // cout << "size = " << testvec.size() << endl;
    // cout << "hei" << endl;

    // // vector<vector<uint> > atom_in_voxels(n_voxels[0]*n_voxels[1]*n_voxels[2]);
    // // cout << "atom_in_voxels.size() = " << atom_in_voxels.size() << endl;
    // cout << "hei" << endl;
    // uint n_tagged_atoms = 0;
    // // DEBUG //

    // Filling atoms_in_voxels
    for (uint atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        // Only putting atoms of the type we want to check against in the list, so we don't have to check for atom type 
        // when looping through the atoms in each voxel later on
        if (mts0_io->atom_types[atom_index] == atom_type_to_check_against) {
            uint i = mts0_io->positions[atom_index][0]/voxel_size_angstrom[0];
            uint j = mts0_io->positions[atom_index][1]/voxel_size_angstrom[1];
            uint k = mts0_io->positions[atom_index][2]/voxel_size_angstrom[2];
            atoms_in_voxels[i][j][k].push_back(atom_index);
        }
    }

    // Loop over all atoms
    double distance_squared = distance*distance;
    uint n_tagged_atoms = 0;
    for (uint atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        if (mts0_io->atom_types[atom_index] == atom_type_to_tag) {
            uint i = mts0_io->positions[atom_index][0]/voxel_size_angstrom[0];
            uint j = mts0_io->positions[atom_index][1]/voxel_size_angstrom[1];
            uint k = mts0_io->positions[atom_index][2]/voxel_size_angstrom[2];

            if (check_if_atom_of_type_exists_within_distance_in_neighbor_voxels(
                mts0_io, 
                atoms_in_voxels, 
                atom_index, 
                n_voxels, 
                i,j,k, 
                distance_squared, 
                system_size_angstrom, 
                half_system_size_angstrom)) {

                mts0_io->atom_types[atom_index] = tag;
                n_tagged_atoms++;
            }

        }
    }

    // if (n_tagged_atoms == 0) {
    //     cout << "Warning: tag_atom_type_within_distance_from_other_atom_type tagged 0 atoms." << endl;
    // }

    return n_tagged_atoms;
}

inline bool check_if_atom_of_type_exists_within_distance_in_neighbor_voxels(
    const Mts0_io* mts0_io, 
    const vector<vector<vector<vector<uint> > > >& atoms_in_voxels, 
    const uint atom_index, 
    const vector<uint>& n_voxels,
    const uint i,
    const uint j,
    const uint k, 
    const double distance_squared,
    const vector<double>& system_size_angstrom,
    const vector<double>& half_system_size_angstrom) {

    // Loop over 27 nearest neighbors (including self (maybe not necessary? but probably faster to not check if di,dj,dk==0))
    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {
                // Periodic boundary conditions
                uint ii = (i + di + n_voxels[0]) % n_voxels[0];
                uint jj = (j + dj + n_voxels[1]) % n_voxels[1];
                uint kk = (k + dk + n_voxels[2]) % n_voxels[2];
                for (auto it = atoms_in_voxels[ii][jj][kk].begin(); it != atoms_in_voxels[ii][jj][kk].end(); ++it) {
                        vector<double> main_atom_position = mts0_io->positions[atom_index];
                        vector<double> position_to_check_against = mts0_io->positions[*it];
                        double dr_squared = calculate_distance_squared_using_minimum_image_convention(main_atom_position, position_to_check_against, system_size_angstrom, half_system_size_angstrom);
                        if (dr_squared < distance_squared) {
                            return true;
                        }
                }
            }
        }
    }

    return false;
}