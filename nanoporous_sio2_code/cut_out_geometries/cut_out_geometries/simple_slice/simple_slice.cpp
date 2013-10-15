#include <iostream>
#include <iomanip> // setprecision, setw
#include <algorithm> // std::min

#include <common_functions.h> // calculate_squared_vector_difference_using_minimum_image_convention
#include <mts0_io.h>

int main(int args, char *argv[]) {
    if (args != 11) {
        cout << "Run program with " << endl;
        cout << "./executable.x  nx  ny  nz  mts0_directory_in  mts0_directory_out  xyz_output_filename  print_mts0_files[0|1]  print_xyz_file[0|1]  ";
        cout << "axis_to_cut_on[0|1|2]  where_on_axis_to_cut" << endl;
        return 1;
    }

    if (!check_files_and_folders(argv)) {
        exit(1);
    }

    int nx = atoi(argv[1]); 
    int ny = atoi(argv[2]); 
    int nz = atoi(argv[3]);
    string mts0_directory_in = argv[4];
    string mts0_directory_out = argv[5];
    string xyz_output_filename = argv[6];
    bool print_mts0_files = atoi(argv[7]);
    bool print_xyz_file = atoi(argv[8]);
    int axis_to_cut_on = atoi(argv[9]);
    double where_on_axis_to_cut = atof(argv[10]);

    if (axis_to_cut_on > 2) {
        cout << "Invalid axis to cut on: " << axis_to_cut_on << ", aborting!" << endl;
        exit(1);
    }

    Mts0_io* mts0_io = new Mts0_io(nx, ny, nz);
    mts0_io->load_atoms(mts0_directory_in);
    int n_atoms = mts0_io->get_number_of_atoms();

    vector<double> system_size = mts0_io->get_lx_ly_lz();
    where_on_axis_to_cut = where_on_axis_to_cut*system_size[axis_to_cut_on];

    char axes[] = {'x', 'y', 'z'};
    cout << "System size   = " << system_size[0] << " " << system_size[1] << " " << system_size[2] << endl;
    cout << "Cutting at " << where_on_axis_to_cut << " on the " << axes[axis_to_cut_on] << "-axis." << endl;

    int n_removed_atoms = 0;
    int n_removed_si = 0;
    int n_removed_o = 0;
    for (int i = 0; i < mts0_io->positions.size(); i++) {
        while (n_removed_o % (2*n_removed_si) != 0) {
            if (mts0_io->positions[i][axis_to_cut_on] > where_on_axis_to_cut) {
                mts0_io->remove_atom(i);
                n_removed_atoms++;
            }
        }
    }

    cout << "Removed " << n_removed_atoms << " atoms (" << setprecision(2) << 100*n_removed_atoms/double(n_atoms) << "%)." << endl;

    save_to_xyz_and_mts0(argv, mts0_io);
}


