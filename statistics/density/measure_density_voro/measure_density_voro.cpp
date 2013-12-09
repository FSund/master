// #pragma once // enabled by default
#include <iostream>
#include <iomanip>
// #include <bitset>       // std::bitset
// #include <cmath>        // ceil, floor
#include <sstream>
#include <algorithm>        // std::find
#include <cmath>            // ceil, floor

#include <mts0_io.h>
// #include "matrix.h"
#include "lib.h"
#include "bit_vector.h"
#include "voro++.hh"

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

uint tag_atom_type_within_distance_from_other_atom_type(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag);
uint tag_atom_type_within_distance_from_other_atom_type_ram(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag);
bool check_if_atom_of_type_exists_within_distance_in_neighbor_voxels(
    const Mts0_io* mts0_io, 
    const vector<vector<vector<vector<uint> > > >& atoms_in_voxels, 
    const uint atom_index, 
    const vector<uint>& n_voxels,
    const uint i,
    const uint j,
    const uint k, 
    const double distance_squared,
    const vector<double> system_size_angstrom,
    const vector<double> half_system_size_angstrom);
double find_number_density_of_atom_type(Mts0_io *mts0_io, int atom_type);

int main(int n_args, char* arg_vec[]) {
    
    if (n_args < 8) {
        // cout << "Arguments: input_mts0_folder  nx  ny  nz  n_voxels_x  n_voxels_y  n_voxels_z  atom_radius" << endl;
        cout << "Arguments: input_mts0_folder  nx  ny  nz  r_start  r_step  n_steps" << endl;
        return 1;
    }

    string input_mts0_folder = arg_vec[1];
    int nx = atoi(arg_vec[2]);
    int ny = atoi(arg_vec[3]);
    int nz = atoi(arg_vec[4]);
    double r_start = atof(arg_vec[5]);
    double r_step = atof(arg_vec[6]);
    double N = atoi(arg_vec[7]);
    // uint n_voxels_x = atoi(arg_vec[5]);
    // uint n_voxels_y = atoi(arg_vec[6]);
    // uint n_voxels_z = atoi(arg_vec[7]);
    // double distance = atof(arg_vec[5]);

    Mts0_io mts0_io(nx, ny, nz);
    mts0_io.load_atoms(input_mts0_folder);


    int pore_atom_type = O_TYPE;
    int matrix_atom_type = SI_TYPE;
    int tag = X_TYPE;

    vector<vector<double> > results(2, vector<double>(N));
    for (uint i = 0; i < N; i++) {
        double distance = r_start + r_step*i;
        uint n_tagged_atoms = tag_atom_type_within_distance_from_other_atom_type_ram(&mts0_io, matrix_atom_type, pore_atom_type, distance, tag);

        double number_density = 0.0;
        if (n_tagged_atoms > 0) {
            number_density = find_number_density_of_atom_type(&mts0_io, tag);
        }

        cout << "r = " << distance << endl;
        cout << "Number of tagged atoms = " << n_tagged_atoms << endl;
        cout << "Number density = " << number_density << endl;

        results[0][i] = distance;
        results[1][i] = number_density;
    }

    // mts0_io.write_to_lammps("test_tagged.lmp");

    return 0;
}

uint tag_atom_type_within_distance_from_other_atom_type(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag) {

    vector<vector<uint> > radial_neighbor_list = create_radial_neighbor_list(mts0_io, distance);

    uint n_tagged_atoms = 0;
    vector<double> system_size = mts0_io->get_lx_ly_lz();
    vector<double> half_system_size = {system_size[0]/2.0, system_size[1]/2.0, system_size[2]/2.0};
    for (int atom_index = 0; atom_index < mts0_io->get_number_of_atoms(); atom_index++) {
        if (mts0_io->atom_types[atom_index] == atom_type_to_tag) {
            for (auto it = radial_neighbor_list[atom_index].begin(); it != radial_neighbor_list[atom_index].end(); ++it) {
                if (mts0_io->atom_types[*it] == atom_type_to_check_against) {
                    mts0_io->atom_types[atom_index] = tag;
                    n_tagged_atoms++;
                    continue;
                }
            }
        }
    }

    if (n_tagged_atoms == 0) {
        cout << "Warning: tag_atom_type_within_distance_from_other_atom_type tagged 0 atoms." << endl;
    }

    return n_tagged_atoms;
}

double find_number_density_of_atom_type(Mts0_io *mts0_io, int atom_type) {

    vector<double> system_size_angstrom = mts0_io->get_lx_ly_lz();

    // Set up constants for the container geometry
    double x_min = 0, x_max = system_size_angstrom[0],
           y_min = 0, y_max = system_size_angstrom[1],
           z_min = 0, z_max = system_size_angstrom[2];

    // Set up the number of blocks that the container is divided into
    int n_x = 6, n_y = 6, n_z = 6;

    // Create a container with the geometry given above, and make it periodic in each of the three coordinates. Allocate 
    // space for "initial_n_particles_in_each_block" particles within each computational block
    int initial_n_particles_in_each_block = 10;
    // voro::container_periodic con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, initial_n_particles_in_each_block); // Couldn't get this to work...
    voro::container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, true, true, true, initial_n_particles_in_each_block);

    // Add atoms into the container, putting the atoms with atom_type equal to a type in "atom_types" into the
    // particle_order po so we can loop through those later
    voro::particle_order po;
    uint n_atoms_of_wanted_type = 0;
    for (int i = 0; i < mts0_io->get_number_of_atoms(); i++) {
        double x = mts0_io->positions[i][0], y = mts0_io->positions[i][1], z = mts0_io->positions[i][2];

        // Check if the atom with index i is the type we're looking for, if so add the atom to the particle_order "po",
        // else just add it to the container (need all atoms to get the correct voronoi cells, but we could probably 
        // get away with just adding the surrounding atoms).
        if (mts0_io->atom_types[i] == atom_type) {
            con.put(po, i, x, y, z);
            n_atoms_of_wanted_type++;
        } else {
            con.put(i, x, y, z);
        }
    }
    // cout << "Total number of atoms of wanted type = " << n_atoms_of_wanted_type << endl;

    // Loop over all points in the particle_order po
    /* 
        From example: http://math.lbl.gov/voro++/examples/loops/
        "
        On line 45, a c_loop_order class is created that is associated with the particle_order. A call is then made to 
        the start function, which sets the loop class to consider the first particle. Usually this returns true, but it 
        may return false in the case when there are no particles to consider.

        If a true result is returned, then a do/while loop is entered. The compute_cell routine is called to calculate 
        the Voronoi cell for the current particle in question. This routine returns true if the Voronoi cell was 
        computed, and false if the cell could not be computed for any reason (such as being completely removed by a 
        wall). Assuming a true result, the routine then queries the position of the current particle (on line 49), and 
        then draws the Voronoi cell as a POV-Ray mesh, and as POV-Ray collection of spheres and cylinders. The inc 
        routine is then called to advance the loop to the next particle to be considered. When the inc routine returns 
        false, there are no more particles to be considered and the routine exits.

        The loop therefore calculates Voronoi cells for all the particles that were previously stored on the 
        particle_order class
        "
    */
    voro::c_loop_order clo(con, po);
    voro::voronoicell cell;
    double total_volume_of_voronoi_cells = 0.0;
    if (clo.start()) { // May return false if there are no particles in particle_order po
        do if (con.compute_cell(cell, clo)) { // con.compute_cell() may return false if the cell couldn't be computed (ex: being completely removed by a wall)
            total_volume_of_voronoi_cells += cell.volume();
        } while (clo.inc()); // Increase c_loop_order to next particle, returns false if we're at the last particle
    }
    // cout << "Total volume of voronoi cells of atom with wanted type = " << total_volume_of_voronoi_cells << endl;

    // // Sum up the volumes, and check that this matches the container volume
    // double system_volume_angstrom = (x_max - x_min)*(y_max - y_min)*(x_max - x_min);
    // double vvol = con.sum_cell_volumes();
    // printf("Container volume : %g\n"
    // "Voronoi volume   : %g\n"
    // "Difference       : %g\n", system_volume_angstrom, vvol, vvol - system_volume_angstrom);

    double number_density = n_atoms_of_wanted_type/total_volume_of_voronoi_cells;
    return number_density;
}

uint tag_atom_type_within_distance_from_other_atom_type_ram(Mts0_io* mts0_io, int atom_type_to_check_against, int atom_type_to_tag, double distance, int tag) {

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

    // Filling atoms_in_voxels
    for (uint atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        // Only putting atoms of the type we want to check against in the list, so we don't have to check when looping
        // through the atoms in each voxel later on
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

    if (n_tagged_atoms == 0) {
        cout << "Warning: tag_atom_type_within_distance_from_other_atom_type tagged 0 atoms." << endl;
    }

    return n_tagged_atoms;
}

bool check_if_atom_of_type_exists_within_distance_in_neighbor_voxels(
    const Mts0_io* mts0_io, 
    const vector<vector<vector<vector<uint> > > >& atoms_in_voxels, 
    const uint atom_index, 
    const vector<uint>& n_voxels,
    const uint i,
    const uint j,
    const uint k, 
    const double distance_squared,
    const vector<double> system_size_angstrom,
    const vector<double> half_system_size_angstrom) {

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
