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
    int min_n_args = 5;
    if (n_args < min_n_args) {
        cout << "Usage: executable.x  input_mts0_folder  nx  ny  nz  atom_type_1  atom_type_2  ..." << endl;
        exit(1);
    }

    string input_mts0_folder = arg_vec[1];
    int nx = atoi(arg_vec[2]);
    int ny = atoi(arg_vec[3]);
    int nz = atoi(arg_vec[4]);

    vector<int> atom_types;
    if (n_args > min_n_args) {
        for (int i = min_n_args; i < n_args; i++) {
            atom_types.push_back(atoi(arg_vec[i]));
        }
        // cout << atom_types.size() << " atom types." << endl;
    } else {
        atom_types = {1, 2, 3, 4, 5, 6, 7};
        cout << "Doing calculations on all atom types." << endl;
    }

    Mts0_io mts0_io(nx, ny, nz);
    mts0_io.load_atoms(input_mts0_folder);

    map<string,int> tag_to_int_map;
    tag_to_int_map["Si"] = SI_TYPE;
    tag_to_int_map["A "]  = A_TYPE;
    tag_to_int_map["H "]  = H_TYPE;
    tag_to_int_map["O "]  = O_TYPE;
    tag_to_int_map["Na"] = NA_TYPE;
    tag_to_int_map["Cl"] = CL_TYPE;
    tag_to_int_map["X "]  = X_TYPE;

    map<int,string> int_to_tag_map;
    int_to_tag_map[SI_TYPE] = "Si";
    int_to_tag_map[A_TYPE]  = "A";
    int_to_tag_map[H_TYPE]  = "H";
    int_to_tag_map[O_TYPE]  = "O";
    int_to_tag_map[NA_TYPE] = "Na";
    int_to_tag_map[CL_TYPE] = "Cl";
    int_to_tag_map[X_TYPE]  = "X";

    vector<double> number_densities(atom_types.size());
    for (uint i = 0; i < atom_types.size(); i++) {
        number_densities[i] = find_number_density_of_atom_type(&mts0_io, atom_types[i]);
        cout << "rho(" << int_to_tag_map[atom_types[i]] << ") = " << number_densities[i] << endl;
    }

    return 0;
}
