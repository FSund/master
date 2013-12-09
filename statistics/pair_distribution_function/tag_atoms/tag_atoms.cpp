#include <iostream>
#include <sstream>

#include <mts0_io.h>

#define SI_TYPE 1
#define A_TYPE 2
#define H_TYPE 3
#define O_TYPE 4
#define NA_TYPE 5
#define CL_TYPE 6
#define X_TYPE 7

using namespace std;

vector<vector<uint> > create_radial_neighbor_list(Mts0_io* mts0_io, const double neighbor_list_radius);
double calculate_distance_squared_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &limits, const vector<double> &half_limits);
vector<double> calculate_vector_difference_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &limits, const vector<double> &half_limits);

int main(int n_args, char* argument_vector[]) {

    if (n_args < 7) {
        cout << "Usage: input_mts0_folder  nx  ny  nz  output_mts0_folder  radius" << endl;
        return 1;
    }

    string input_mts0_folder = argument_vector[1];
    int nx = atoi(argument_vector[2]);
    int ny = atoi(argument_vector[3]);
    int nz = atoi(argument_vector[4]);
    string output_mts0_folder = argument_vector[5];
    double radius = atof(argument_vector[6]);

    Mts0_io mts0_io(nx, ny
    , nz);
    mts0_io.load_atoms(input_mts0_folder);
    vector<vector<uint> > radial_neighbor_list = create_radial_neighbor_list(&mts0_io, radius);
    // double radius_squared = radius*radius;

    vector<double> system_size = mts0_io.get_lx_ly_lz();
    vector<double> half_system_size = {system_size[0]/2.0, system_size[1]/2.0, system_size[2]/2.0};
    for (int atom_index = 0; atom_index < mts0_io.get_number_of_atoms(); atom_index++) {
        if (mts0_io.atom_types[atom_index] == O_TYPE) {

            // Loop through neighbors to see if we are close enough to the surface (within "radius")
            for (auto it = radial_neighbor_list[atom_index].begin(); it != radial_neighbor_list[atom_index].end(); ++it) {
                // double dr_squared = calculate_distance_squared_using_minimum_image_convention(mts0_io.positions[atom_index], mts0_io.positions[*it], system_size, half_system_size);

                // If the atom with index atom_index has an atom of type A_TYPE in the list of neighbors
                if (mts0_io.atom_types[*it] == A_TYPE) {
                    mts0_io.atom_types[atom_index] = X_TYPE;
                    continue;
                }
            }
        }
    }

    mts0_io.save_atoms(output_mts0_folder);
    stringstream filename_gen;
    filename_gen << output_mts0_folder << ".lmp";
    mts0_io.write_to_lammps(filename_gen.str());

    return 0;
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

double calculate_distance_squared_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &limits, const vector<double> &half_limits) {
    vector<double> dr_vector = calculate_vector_difference_using_minimum_image_convention(u, v, limits, half_limits);

    return dr_vector[0]*dr_vector[0] + dr_vector[1]*dr_vector[1] + dr_vector[2]*dr_vector[2];
}

vector<double> calculate_vector_difference_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &limits, const vector<double> &half_limits){
    vector<double> dr_vector(3);
    dr_vector[0] = u[0] - v[0];
    dr_vector[1] = u[1] - v[1];
    dr_vector[2] = u[2] - v[2];

    // minimum image convention
    if      (dr_vector[0] >  half_limits[0]) dr_vector[0] -= limits[0];
    else if (dr_vector[0] < -half_limits[0]) dr_vector[0] += limits[0];
    if      (dr_vector[1] >  half_limits[1]) dr_vector[1] -= limits[1];
    else if (dr_vector[1] < -half_limits[1]) dr_vector[1] += limits[1];
    if      (dr_vector[2] >  half_limits[2]) dr_vector[2] -= limits[2];
    else if (dr_vector[2] < -half_limits[2]) dr_vector[2] += limits[2];

    return dr_vector;
}
