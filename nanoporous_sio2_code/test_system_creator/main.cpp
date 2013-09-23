#include <iostream>
#include <fstream>
#include <sstream>
#include "creator/creator.h"

using namespace std;

int main(int args, char *argv[]) {
    if(args != 9) {
        cout << "Run program with " << endl;
        cout << "./executable.x  nx  ny  nz  mts0_directory_in  mts0_directory_out  xyz_output_filename  print_mts0_files[0|1]  print_xyz_file[0|1]" << endl;
        return 0;
    }

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int nz = atoi(argv[3]);
    string mts0_directory_in = argv[4];
    string mts0_directory_out = argv[5];
    string xyz_output_filename = argv[6];
    bool print_mts0_files = atoi(argv[7]);
    bool print_xyz_file = atoi(argv[8]);

    // Testing that we have write permissions, and that the needed output folder exists
    if (print_mts0_files) {
        ostringstream filename;
        filename << mts0_directory_out << "/mt0000" << endl;

        ofstream file(filename.str().c_str());
        if (!file.is_open()) {
            cout << "! Error: couldn't open file " << filename.str() << " writing, please chech that the mts0_directory_out " << mts0_directory_out << " exists. Aborting!" << endl;
            exit(1);
        }
        file.close();
    }
    if (print_xyz_file) {
        ofstream file(xyz_output_filename.c_str());
        if (!file.is_open()) {
            cout << "! Error: couldn't open file " << xyz_output_filename << " writing, please chech that you have write permissions. Aborting!" << endl;
            exit(1);
        }
        file.close();
    }

    Creator creator(mts0_directory_in, nx, ny, nz);

    creator.remove_all_atoms();
//    creator.remove_all_atoms_in_half_system();
//    creator.cut_system_along_space_diagonal();

    creator.add_tetrahedra();
//    creator.test_mathildes_function();

    if (print_mts0_files) {
        cout << "Saving mts0 files ..." << endl;
        creator.mts0_io->save_atoms(mts0_directory_out);
    }

    if (print_xyz_file) {
        cout << "Saving to .xyz-file ..." << endl;
        creator.mts0_io->write_to_xyz_file(xyz_output_filename);
        cout << "Done saving to .xyz-file." << endl;
    }

    return 0;
}
