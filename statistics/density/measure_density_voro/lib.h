#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <voro++.hh>
#include <mts0_io.h>
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
    double total_volume_of_voronoi_cells_of_wanted_atom_type = 0.0;
    if (clo.start()) { // May return false if there are no particles in particle_order po
        do if (con.compute_cell(cell, clo)) { // con.compute_cell() may return false if the cell couldn't be computed (ex: being completely removed by a wall)
            total_volume_of_voronoi_cells_of_wanted_atom_type += cell.volume();
        } while (clo.inc()); // Increase c_loop_order to next particle, returns false if we're at the last particle
    }
    // cout << "Total volume of voronoi cells of atom with wanted type = " << total_volume_of_voronoi_cells_of_wanted_atom_type << endl;

    // // Sum up the volumes, and check that this matches the container volume
    // double system_volume_angstrom = (x_max - x_min)*(y_max - y_min)*(x_max - x_min);
    // double vvol = con.sum_cell_volumes();
    // printf("Container volume : %g\n"
    // "Voronoi volume   : %g\n"
    // "Difference       : %g\n", system_volume_angstrom, vvol, vvol - system_volume_angstrom);

    double number_density;
    if (n_atoms_of_wanted_type == 0) {
        number_density = 0.0;
    } else {
        number_density = n_atoms_of_wanted_type/total_volume_of_voronoi_cells_of_wanted_atom_type;
    }
    return number_density;
}
