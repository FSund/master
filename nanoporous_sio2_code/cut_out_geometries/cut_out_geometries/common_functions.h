#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib> // atoi, atof

#include <mts0_io.h>

using namespace std;

int check_files_and_folders(char *argv[]);
int save_to_xyz_and_mts0(char *argv[], Mts0_io* mts0_io);
vector<double> calculate_vector_difference_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &system_size, const vector<double> &half_system_size);
inline double calculate_squared_vector_difference_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &system_size, const vector<double> &half_system_size);

int check_files_and_folders(char *argv[]) {
    string mts0_directory_out = argv[5];
    string xyz_output_filename = argv[6];
    bool print_mts0_files = atoi(argv[7]);
    bool print_xyz_file = atoi(argv[8]);

    bool tests_went_ok = true;

    // Testing that we have write permissions, and that the needed output folder exists
    if (print_mts0_files) {
        ostringstream filename;
        filename << mts0_directory_out << "/mt0000";

        ofstream file(filename.str().c_str());
        if (!file.is_open()) {
            cout << "! Error: couldn't open file \"" << filename.str() << "\" for writing, please chech that the mts0_directory_out \"" << mts0_directory_out << "\" exists." << endl;
            tests_went_ok = false;
        }
        file.close();
    }
    if (print_xyz_file) {
        ofstream file(xyz_output_filename.c_str());
        if (!file.is_open()) {
            cout << "! Error: couldn't open file \"" << xyz_output_filename << "\" for writing, please chech that you have write permissions." << endl;
            tests_went_ok = false;
        }
        file.close();
    }

    return tests_went_ok;
}

int save_to_xyz_and_mts0(char *argv[], Mts0_io* mts0_io) {
    string mts0_directory_out = argv[5];
    string xyz_output_filename = argv[6];
    bool print_mts0_files = atoi(argv[7]);
    bool print_xyz_file = atoi(argv[8]);

    if (print_mts0_files) {
        cout << "Saving mts0 files ..." << endl;
        mts0_io->save_atoms(mts0_directory_out);
    }

    if (print_xyz_file) {
        cout << "Saving to .xyz-file ..." << endl;
        mts0_io->write_to_xyz_file(xyz_output_filename);
        cout << "Done saving to .xyz-file." << endl;
    }
}

inline double calculate_squared_vector_difference_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &system_size, const vector<double> &half_system_size){
    vector<double> dr = calculate_vector_difference_using_minimum_image_convention(u, v, system_size, half_system_size);
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}

vector<double> calculate_vector_difference_using_minimum_image_convention(const vector<double> &u, const vector<double> &v, const vector<double> &system_size, const vector<double> &half_system_size){
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


#endif // COMMON_FUNCTIONS_H