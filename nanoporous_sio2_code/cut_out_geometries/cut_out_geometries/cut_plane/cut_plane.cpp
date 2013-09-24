#include <iostream>
#include <iomanip> // setprecision, setw
#include <algorithm> // std::min

#include <common_functions.h> // calculate_squared_vector_difference_using_minimum_image_convention
#include <mts0_io.h>

int main(int args, char *argv[]) {
    if (args != 16) {
        cout << "Run program with " << endl;
        cout << "./executable.x  nx  ny  nz  mts0_directory_in  mts0_directory_out  xyz_output_filename  print_mts0_files[0|1]  print_xyz_file[0|1]  ";
        cout << "x0  y0  z0(point on plane, [0,1])  normal_x  normal_y  normal_z(normal to plane) ";
        cout << "keep_atoms_on_which_side_of_plane[-1|+1]" << endl;
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
    double x0 = atof(argv[9]);
    double y0 = atof(argv[10]);
    double z0 = atof(argv[11]);
    double normal_x = atof(argv[12]);
    double normal_y = atof(argv[13]);
    double normal_z = atof(argv[14]);
    int sign_of_atoms_to_keep = atoi(argv[15]);

    Mts0_io* mts0_io = new Mts0_io(nx, ny, nz);
    mts0_io->load_atoms(mts0_directory_in);
    int n_atoms = mts0_io->get_number_of_atoms();

    vector<double> system_size = mts0_io->get_lx_ly_lz();
    x0 *= system_size[0];
    y0 *= system_size[1];
    z0 *= system_size[2];

    cout << "System size = " << system_size[0] << " " << system_size[1] << " " << system_size[2] << endl;
    cout << "Cutting plane defined by point at " << x0 << " " << y0 << " " << z0 << ", and normal " << normal_x << " " << normal_y << " " << normal_z << "." << endl;
    if (sign_of_atoms_to_keep == 1) cout << "Keeping atoms on the side of the plane that the normal points towards." << endl;
    else cout << "Keeping atoms on the opposite side of the plane that the normal points towards." << endl;

    int n_removed_atoms = 0;
    for (int i = 0; i < mts0_io->positions.size(); i++) {
        double x = mts0_io->positions[i][0];
        double y = mts0_io->positions[i][1];
        double z = mts0_io->positions[i][2];
        int point_is_on_right_side_of_plane = ( normal_x*(x-x0) + normal_y*(y-y0) + normal_z*(z-z0) ) > 0;
        if (sign_of_atoms_to_keep == -1) {
            if (point_is_on_right_side_of_plane) {
                mts0_io->remove_atom(i);
                n_removed_atoms++;
            }
        } else {
            if (!point_is_on_right_side_of_plane) {
                mts0_io->remove_atom(i);
                n_removed_atoms++;
            }
        }
    }

    cout << "Removed " << n_removed_atoms << " atoms (" << setprecision(2) << 100*n_removed_atoms/double(n_atoms) << "%)." << endl;

    save_to_xyz_and_mts0(argv, mts0_io);
}
