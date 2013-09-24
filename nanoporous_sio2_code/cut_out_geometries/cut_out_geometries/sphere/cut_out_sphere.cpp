#include <iostream>
#include <iomanip> // setprecision, setw
#include <algorithm> // std::min

#include <common_functions.h> // calculate_squared_vector_difference_using_minimum_image_convention
#include <mts0_io.h>

int main(int args, char *argv[]) {
    if (args != 13) {
        cout << "Run program with " << endl;
        cout << "./executable.x  nx  ny  nz  mts0_directory_in  mts0_directory_out  xyz_output_filename  print_mts0_files[0|1]  print_xyz_file[0|1]  ";
        cout << "sphere_center_x/system_size_x  sphere_center_y/system_size_y  sphere_center_z/system_size_z  sphere_radius/min(system_size)" << endl;
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
    double sphere_center_x = atof(argv[9]);
    double sphere_center_y = atof(argv[10]);
    double sphere_center_z = atof(argv[11]);
    double sphere_radius = atof(argv[12]);

    Mts0_io* mts0_io = new Mts0_io(nx, ny, nz);
    mts0_io->load_atoms(mts0_directory_in);
    int n_atoms = mts0_io->get_number_of_atoms();

    vector<double> system_size = mts0_io->get_lx_ly_lz();
    vector<double> half_system_size(3);
    half_system_size[0] = system_size[0]/2.0;
    half_system_size[1] = system_size[1]/2.0;
    half_system_size[2] = system_size[2]/2.0;
    vector<double> sphere_center(3);
    sphere_center[0] = system_size[0]*sphere_center_x;
    sphere_center[1] = system_size[1]*sphere_center_y;
    sphere_center[2] = system_size[2]*sphere_center_z;

    double min_system_length = min(
            min(system_size[0], system_size[1]),
            min(system_size[1], system_size[2])
        );
    sphere_radius = sphere_radius*system_size[0];

    cout << "System size   = " << system_size[0] << " " << system_size[1] << " " << system_size[2] << endl;
    cout << "Sphere center = " << sphere_center[0] << " " << sphere_center[1] << " " << sphere_center[2] << endl;
    cout << "Sphere radius = " << sphere_radius << endl;

    int n_removed_atoms = 0;
    for (int i = 0; i < mts0_io->positions.size(); i++) {
        double dr2 = calculate_squared_vector_difference_using_minimum_image_convention(mts0_io->positions[i], sphere_center, system_size, half_system_size);
        if (dr2 < sphere_radius*sphere_radius) {
            mts0_io->remove_atom(i);
            n_removed_atoms++;
        }
    }

    cout << "Removed " << n_removed_atoms << " atoms (" << setprecision(2) << 100*n_removed_atoms/double(n_atoms) << "%)." << endl;

    save_to_xyz_and_mts0(argv, mts0_io);
}
