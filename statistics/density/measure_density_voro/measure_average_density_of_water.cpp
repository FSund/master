#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

#include "lib.h"

using namespace std;

typedef unsigned int uint;

#define SI_TYPE 1
#define A_TYPE 2
#define H_TYPE 3
#define O_TYPE 4
#define NA_TYPE 5
#define CL_TYPE 6
#define X_TYPE 7

int main(int n_args, char* arg_vec[]) {
    int min_n_args = 6;
    if (n_args < min_n_args) {
        cout << "Usage: executable.x  input_mts0_folder  nx  ny  nz  min_distance_from_matrix" << endl;
        exit(1);
    }

    string input_mts0_folder = arg_vec[1];
    int nx = atoi(arg_vec[2]);
    int ny = atoi(arg_vec[3]);
    int nz = atoi(arg_vec[4]);
    double min_distance_from_matrix = atof(arg_vec[5]);

    Mts0_io mts0_io(nx, ny, nz);
    mts0_io.load_atoms(input_mts0_folder);

    // Remove hydrogen atoms (might be a bad idea)
    for (uint i = 0; i < mts0_io.positions.size(); i++) {
        if (mts0_io.atom_types[i] == H_TYPE) {
            mts0_io.remove_atom(i);
        }
    }
    mts0_io.rearrange_vectors_by_moving_atoms_from_end_to_locations_where_atoms_have_been_removed(); // to update mts0_io

    int pore_atom_type = O_TYPE;
    int matrix_atom_type = SI_TYPE;
    int tag = 10;

    // Tag O-type atoms too close to matrix with another tag
    tag_atom_type_within_distance_from_other_atom_type(&mts0_io, matrix_atom_type, pore_atom_type, min_distance_from_matrix, tag);

    double number_density = find_number_density_of_atom_type(&mts0_io, pore_atom_type);
    cout << "Number density of water [atoms/Angstrom^3] = " << number_density << endl;

    double M = 0.0180158;       // kg/mol
    double NA = 6.0221413e+23;  // Avogadro's number
    cout << "Density [kg/m^3] = " << number_density*M/NA*(1.0e30) << endl;

    return 0;
}
