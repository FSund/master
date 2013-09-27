#include "creator.h"

Creator::Creator(string mts0_directory_in, int nx, int ny, int nz, double Si_O_distance):
    SI_O_DISTANCE(Si_O_distance)
{
    mts0_io = new Mts0_io(nx,ny,nz);
    cout << "Loading mts0 files ..." << endl;
    mts0_io->load_atoms(mts0_directory_in);
    cout << "Done loading mts0 files." << endl;

    system_size = mts0_io->get_lx_ly_lz();
    half_system_size.resize(3);
    half_system_size[0] = system_size[0]/2.0;
    half_system_size[1] = system_size[1]/2.0;
    half_system_size[2] = system_size[2]/2.0;

    cout << "System size = " << system_size[0] << " " << system_size[1] << " " << system_size[2] << endl;
}

void Creator::make_test_system_with_four_types_of_tetrahedra() {

    vector<Atom> tetrahedron_atoms;
    vector<double> tetrahedron_center(3);
    vector<double> normal_vector(3);
    int n_oxygen;
    double a = 5.0;

    //
    tetrahedron_center[0] = system_size[0]/2.0 + a;
    tetrahedron_center[1] = system_size[1]/2.0 + a;
    tetrahedron_center[2] = system_size[2]/2.0;
    cout << "Tetrahedron #1 center = " << tetrahedron_center[0] << " " << tetrahedron_center[1] << " " << tetrahedron_center[2] << endl;

//    normal_vector[0] = 1.1;
//    normal_vector[1] = -0.67;
//    normal_vector[2] = 1.9;
    normal_vector[0] = 0.0;
    normal_vector[1] = 0.0;
    normal_vector[2] = 1.0;
    normal_vector = normalize_vector(normal_vector);

    n_oxygen = 1;
    tetrahedron_atoms = make_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE, n_oxygen);
    add_atoms(tetrahedron_atoms);

    //
    tetrahedron_center[0] = system_size[0]/2.0 + a;
    tetrahedron_center[1] = system_size[1]/2.0 - a;
    tetrahedron_center[2] = system_size[2]/2.0;
    cout << "Tetrahedron #2 center = " << tetrahedron_center[0] << " " << tetrahedron_center[1] << " " << tetrahedron_center[2] << endl;

//    normal_vector[0] = 2.1;
//    normal_vector[1] = 0.67;
//    normal_vector[2] = -1.9;
    normal_vector[0] = 0.0;
    normal_vector[1] = 0.0;
    normal_vector[2] = 1.0;
    normal_vector = normalize_vector(normal_vector);

    n_oxygen = 2;
    tetrahedron_atoms = make_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE, n_oxygen);
    add_atoms(tetrahedron_atoms);

    //
    tetrahedron_center[0] = system_size[0]/2.0 - a;
    tetrahedron_center[1] = system_size[1]/2.0 - a;
    tetrahedron_center[2] = system_size[2]/2.0;
    cout << "Tetrahedron #3 center = " << tetrahedron_center[0] << " " << tetrahedron_center[1] << " " << tetrahedron_center[2] << endl;

//    normal_vector[0] = -0.1;
//    normal_vector[1] = 2.243;
//    normal_vector[2] = -1.54;
    normal_vector[0] = 0.0;
    normal_vector[1] = 0.0;
    normal_vector[2] = 1.0;
    normal_vector = normalize_vector(normal_vector);

    n_oxygen = 3;
    tetrahedron_atoms = make_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE, n_oxygen);
    add_atoms(tetrahedron_atoms);

    //
    tetrahedron_center[0] = system_size[0]/2.0 - a;
    tetrahedron_center[1] = system_size[1]/2.0 + a;
    tetrahedron_center[2] = system_size[2]/2.0;
    cout << "Tetrahedron #4 center = " << tetrahedron_center[0] << " " << tetrahedron_center[1] << " " << tetrahedron_center[2] << endl;

//    normal_vector[0] = 1.1;
//    normal_vector[1] = 1.243;
//    normal_vector[2] = 0.454;
    normal_vector[0] = 0.0;
    normal_vector[1] = 0.0;
    normal_vector[2] = 1.0;
    normal_vector = normalize_vector(normal_vector);

    n_oxygen = 4;
    tetrahedron_atoms = make_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE, n_oxygen);
    add_atoms(tetrahedron_atoms);
}

void Creator::make_test_system_with_close_tetrahedra()
{
    vector<Atom> tetrahedron_atoms;
    vector<double> tetrahedron_center(3);
    vector<double> normal_vector(3);
    int n_oxygen;
    double a = 8.0;

    //
    tetrahedron_center[0] = system_size[0]/a;
    tetrahedron_center[1] = system_size[1]/a;
    tetrahedron_center[2] = system_size[2]/a;

    normal_vector[0] = -1.0;
    normal_vector[1] = -1.0;
    normal_vector[2] = -1.0;
    normal_vector = normalize_vector(normal_vector);

    n_oxygen = 3;
    tetrahedron_atoms = make_inverse_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE, n_oxygen);
    add_atoms(tetrahedron_atoms);

    //
    double b = 1.0;
    tetrahedron_center[0] += b*SI_O_DISTANCE;
    tetrahedron_center[1] += b*SI_O_DISTANCE;
    tetrahedron_center[2] += b*SI_O_DISTANCE;

    normal_vector[0] = 1.0;
    normal_vector[1] = 1.0;
    normal_vector[2] = 1.0;
    normal_vector = normalize_vector(normal_vector);

    n_oxygen = 3;
    tetrahedron_atoms = make_inverse_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE, n_oxygen);
    add_atoms(tetrahedron_atoms);
}

vector<Atom> Creator::make_tetrahedron(const vector<double> &tetrahedron_center, const vector<double> &normal_vector, double Si_O_distance, int n_oxygen) {
    if (n_oxygen > 4 || n_oxygen < 1) {
        cout << "Error: the tetrahedron needs between 1 and 4 oxygen, aborting!" << endl;
        exit(1);
    }

    vector<Atom> atoms(1 + n_oxygen); // one Silicon and four Oxygen

    // Silicon
    atoms[0].type = SI_TYPE;
    atoms[0].position = tetrahedron_center;

    // First oxygen
    atoms[1].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[1].position[i] = tetrahedron_center[i] - normal_vector[i]*Si_O_distance;
    }

    if (n_oxygen == 1) return atoms;

    // Vector from silicon to a point in the plane spanned by the three missing oxygen atoms
    vector<double> silicon_to_point_in_plane(3);
    double distance_from_silicon_to_plane = Si_O_distance*cos(PI - TETRAHEDRAL_ANGLE_RAD);
    silicon_to_point_in_plane[0] = normal_vector[0]*distance_from_silicon_to_plane;
    silicon_to_point_in_plane[1] = normal_vector[1]*distance_from_silicon_to_plane;
    silicon_to_point_in_plane[2] = normal_vector[2]*distance_from_silicon_to_plane;

    // Vector from point in plane to second oxygen (we can choose this positions freely, it only needs to lie in the plane)
    vector<double> point_in_plane_to_second_oxygen = find_vector_in_plane(normal_vector);
    point_in_plane_to_second_oxygen = normalize_vector(point_in_plane_to_second_oxygen);
    double distance_from_point_in_plane_to_second_oxygen = Si_O_distance*sin(PI - TETRAHEDRAL_ANGLE_RAD);
    point_in_plane_to_second_oxygen[0] *= distance_from_point_in_plane_to_second_oxygen;
    point_in_plane_to_second_oxygen[1] *= distance_from_point_in_plane_to_second_oxygen;
    point_in_plane_to_second_oxygen[2] *= distance_from_point_in_plane_to_second_oxygen;

    // Second oxygen
    atoms[2].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[2].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] + point_in_plane_to_second_oxygen[i];
    }

    if (n_oxygen == 2) return atoms;

    // Unit vector from point in plane, normal to vector from point in plane to second oxygen
    vector<double> unit_vector_in_plane = vector_cross_product(point_in_plane_to_second_oxygen, normal_vector);
    unit_vector_in_plane = normalize_vector(unit_vector_in_plane);

    // Unit vector from point in plane, parallel to vector from point in plane to second oxygen
    vector<double> unit_vector_from_point_to_second_oxygen = normalize_vector(point_in_plane_to_second_oxygen);

    double Lx = distance_from_point_in_plane_to_second_oxygen*sin(60.0*PI/180.0);
    double Ly = distance_from_point_in_plane_to_second_oxygen*cos(60.0*PI/180.0);

    // third oxygen
    atoms[3].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[3].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] - unit_vector_from_point_to_second_oxygen[i]*Ly + unit_vector_in_plane[i]*Lx;
    }
    if (n_oxygen == 3) return atoms;

    // fourth oxygen
    atoms[4].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[4].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] - unit_vector_from_point_to_second_oxygen[i]*Ly - unit_vector_in_plane[i]*Lx;
    }

    return atoms;
}

vector<Atom> Creator::make_inverse_tetrahedron(const vector<double> &tetrahedron_center, const vector<double> &normal_vector, double Si_O_distance, int n_oxygen) {
    if (n_oxygen > 4 || n_oxygen < 1) {
        cout << "Error: the tetrahedron needs between 1 and 4 oxygen, aborting!" << endl;
        exit(1);
    }

    vector<Atom> atoms(1 + n_oxygen); // one Silicon and four Oxygen

    // Silicon
    atoms[0].type = SI_TYPE;
    atoms[0].position = tetrahedron_center;

    // Vector from silicon to a point in the plane spanned by the three missing oxygen atoms
    vector<double> silicon_to_point_in_plane(3);
    double distance_from_silicon_to_plane = Si_O_distance*cos(PI - TETRAHEDRAL_ANGLE_RAD);
    silicon_to_point_in_plane[0] = normal_vector[0]*distance_from_silicon_to_plane;
    silicon_to_point_in_plane[1] = normal_vector[1]*distance_from_silicon_to_plane;
    silicon_to_point_in_plane[2] = normal_vector[2]*distance_from_silicon_to_plane;

    // Vector from point in plane to second oxygen (we can choose this positions freely, it only needs to lie in the plane)
    vector<double> point_in_plane_to_first_oxygen = find_vector_in_plane(normal_vector);
    point_in_plane_to_first_oxygen = normalize_vector(point_in_plane_to_first_oxygen);
    double distance_from_point_in_plane_to_first_oxygen = Si_O_distance*sin(PI - TETRAHEDRAL_ANGLE_RAD);
    point_in_plane_to_first_oxygen[0] *= distance_from_point_in_plane_to_first_oxygen;
    point_in_plane_to_first_oxygen[1] *= distance_from_point_in_plane_to_first_oxygen;
    point_in_plane_to_first_oxygen[2] *= distance_from_point_in_plane_to_first_oxygen;

    // first oxygen
    atoms[1].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[1].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] + point_in_plane_to_first_oxygen[i];
    }

    if (n_oxygen == 1) return atoms;

    // Unit vector from point in plane, normal to vector from point in plane to first oxygen
    vector<double> unit_vector_in_plane = vector_cross_product(point_in_plane_to_first_oxygen, normal_vector);
    unit_vector_in_plane = normalize_vector(unit_vector_in_plane);

    // Unit vector from point in plane, parallel to vector from point in plane to first oxygen
    vector<double> unit_vector_from_point_to_second_oxygen = normalize_vector(point_in_plane_to_first_oxygen);

    double Lx = distance_from_point_in_plane_to_first_oxygen*sin(60.0*PI/180.0);
    double Ly = distance_from_point_in_plane_to_first_oxygen*cos(60.0*PI/180.0);

    // second oxygen
    atoms[2].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[2].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] - unit_vector_from_point_to_second_oxygen[i]*Ly + unit_vector_in_plane[i]*Lx;
    }
    if (n_oxygen == 2) return atoms;

    // third oxygen
    atoms[3].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[3].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] - unit_vector_from_point_to_second_oxygen[i]*Ly - unit_vector_in_plane[i]*Lx;
    }
    if (n_oxygen == 3) return atoms;

    // fourth oxygen
    atoms[4].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[4].position[i] = tetrahedron_center[i] - normal_vector[i]*Si_O_distance;
    }

    return atoms;
}

void Creator::remove_all_atoms() {
    for (int atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        mts0_io->remove_atom(atom_index);
    }
    cout << "Removed all atoms in the system." << endl;
}

void Creator::add_atoms(const vector<Atom> &atoms) {
    for (int i = 0; i < atoms.size(); i++) {
        add_atom(atoms[i]);
    }
}

void Creator::add_atom(const Atom &atom) {
    if (atom.type == 0) {
        cout << "! Error: you forgot to set the type of an atom before adding it, skipping the atom!" << endl;
        return;
    }
    mts0_io->add_atom(atom.type, atom.position, atom.velocity);
}

double Creator::vector_dot_product(const vector<double> &u, const vector<double> &v) {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

vector<double> Creator::vector_cross_product(const vector<double> &u, const vector<double> &v) {
    vector<double> cross_product(3);
    cross_product[0] = u[1]*v[2] - u[2]*v[1];
    cross_product[1] = u[2]*v[0] - u[0]*v[2];
    cross_product[2] = u[0]*v[1] - u[1]*v[0];

    return cross_product;
}

vector<double> Creator::normalize_vector(const vector<double> &u) {
    for (int i = 0; i < 3; i++) {
        if (!isfinite(u[i])) {
            cout << "Error: Trying to normalize non-finite vector. Aborting!" << endl;
            cout << u[i] << endl;
            exit(1);
        }
    }
    double inverse_norm = 1.0/sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    vector<double> normalized_vector = u;
    normalized_vector[0] *= inverse_norm;
    normalized_vector[1] *= inverse_norm;
    normalized_vector[2] *= inverse_norm;

    return normalized_vector;
}

vector<double> Creator::calculate_vector_difference_using_minimum_image_convention(const vector<double> &position1, const vector<double> &position2) {
    vector<double> dr_vector(3);
    dr_vector[0] = position1[0] - position2[0];
    dr_vector[1] = position1[1] - position2[1];
    dr_vector[2] = position1[2] - position2[2];

    // minimum image convention
    if      (dr_vector[0] >  half_system_size[0]) dr_vector[0] -= system_size[0];
    else if (dr_vector[0] < -half_system_size[0]) dr_vector[0] += system_size[0];
    if      (dr_vector[1] >  half_system_size[1]) dr_vector[1] -= system_size[1];
    else if (dr_vector[1] < -half_system_size[1]) dr_vector[1] += system_size[1];
    if      (dr_vector[2] >  half_system_size[2]) dr_vector[2] -= system_size[2];
    else if (dr_vector[2] < -half_system_size[2]) dr_vector[2] += system_size[2];

    return dr_vector;
}

vector<double> Creator::find_vector_in_plane(const vector<double> &normal_vector)
{
    // Check if any of the components of the normal vector is close or equal to zero
    vector<bool> small_component_in_normal_vector(3, false);
    int number_of_small_components_in_normal_vector = 0;
    double epsilon = 10e-10;
    for (int i = 0; i < 3; i++) {
        if (abs(normal_vector[i]) < epsilon) {
            small_component_in_normal_vector[i] = true;
            number_of_small_components_in_normal_vector++;
        }
    }

    // Creating vector in the plane defined by normal_vector
    vector<double> u(3);
    switch(number_of_small_components_in_normal_vector) {
    case 0:
        u[0] = 1.0;
        u[1] = 1.0;
        u[2] = -((normal_vector[0] + normal_vector[1])/normal_vector[2]);
        break;
    case 1:
        if (small_component_in_normal_vector[0]) {
            u[0] = 1.0;
            u[0] = 1.0;
            u[0] = 1.0 - ((normal_vector[0] + normal_vector[1])/normal_vector[2]);
        } else {
            u[0] = 1.0 -((normal_vector[1] + normal_vector[2])/normal_vector[0]);
            u[0] = 1.0;
            u[0] = 1.0;
        }
        break;
    case 2:
        u[0] = ( small_component_in_normal_vector[0] ? 1.0 : (normal_vector[1] + normal_vector[2])/normal_vector[0] );
        u[1] = ( small_component_in_normal_vector[1] ? 1.0 : (normal_vector[0] + normal_vector[2])/normal_vector[1] );
        u[2] = ( small_component_in_normal_vector[2] ? 1.0 : (normal_vector[0] + normal_vector[1])/normal_vector[2] );
        break;
    }

    return u;
}

void Creator::ensure_within_periodic_boundaries(vector<double> &position) {
    if      (position[0] < 0.0)            position[0] -= floor(position[0]/system_size[0])*system_size[0];
    else if (position[0] > system_size[0]) position[0] -= floor(position[0]/system_size[0])*system_size[0];
    if      (position[1] < 0.0)            position[1] -= floor(position[1]/system_size[1])*system_size[1];
    else if (position[1] > system_size[1]) position[1] -= floor(position[1]/system_size[1])*system_size[1];
    if      (position[2] < 0.0)            position[2] -= floor(position[2]/system_size[2])*system_size[2];
    else if (position[2] > system_size[2]) position[2] -= floor(position[2]/system_size[2])*system_size[2];
}
