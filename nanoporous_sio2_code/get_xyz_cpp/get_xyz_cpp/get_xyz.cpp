#include <iostream>
#include <sstream>

#include <mts0_io.h>

int main(int args, char *argv[]) {
    if (args < 5) {
        cout << "Run program with " << endl;
        cout << "./executable.x  nx  ny  nz  mts0_directory_in  [xyz_output_filename](will put .xyz file in mts0 dir)" << endl;
        return 1;
    }

    int nx = atoi(argv[1]); 
    int ny = atoi(argv[2]); 
    int nz = atoi(argv[3]);
    string mts0_directory_in = argv[4];
    string xyz_output_filename;
    if (args == 5) {
        ostringstream filename;
        filename << mts0_directory_in << "/state.xyz" << endl;
        xyz_output_filename = filename.str();
        cout << "Will save .xyz file to \"" << xyz_output_filename << "\"." << endl;
    } else if (args == 6) {
        xyz_output_filename = argv[5];
    }

    Mts0_io* mts0_io = new Mts0_io(nx, ny, nz);

    cout << "Loading atoms from mts0-files in folder \"" << mts0_directory_in << "\" ..." << endl;
    mts0_io->load_atoms(mts0_directory_in);
    cout << "Done loading atoms." << endl;

    vector<double> system_size = mts0_io->get_lx_ly_lz();
    
    cout << endl;
    cout << "System size = " << system_size[0] << " " << system_size[1] << " " << system_size[2] << endl;
    cout << "Number of atoms = " << mts0_io->get_number_of_atoms() << endl;
    cout << endl;

    cout << "Saving to .xyz-file ..." << endl;
    mts0_io->write_to_xyz_file(xyz_output_filename);
    cout << "Done saving to .xyz-file." << endl;
}
