#include "creator.h"

Creator::Creator(string mts0_directory_in, int nx, int ny, int nz)
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

//    // needed for Mathilde's functions
//    //calculate the linear combination a*n + b*plane_vec for two missing oxygen method
//    weight_normal = SI_O_DISTANCE*sin(TETRAHEDRAL_ANGLE_RAD*0.5);
//    weight_vec_sum_existing_oxygen = SI_O_DISTANCE*cos(TETRAHEDRAL_ANGLE_RAD*0.5);

//    //calculate needed lengths for three missing oxygen method
//    length_from_silicon_to_point_in_plane = SI_O_DISTANCE*cos(PI - TETRAHEDRAL_ANGLE_RAD);
//    length_from_point_in_plane_to_new_oxygen = SI_O_DISTANCE*sin(PI - TETRAHEDRAL_ANGLE_RAD);
}

//vector<Atom> Creator::make_tetrahedron_mathildes_function(const vector<double> &tetrahedron_center, const vector<double> &normal_vector, double Si_O_distance) {
//    vector<Atom> atoms(2);

//    // Silicon
//    atoms[0].type = SI_TYPE;
//    atoms[0].position = tetrahedron_center;

//    // First oxygen
//    atoms[1].type = A_TYPE;
//    for (int i = 0; i < 3; i++) {
//        atoms[1].position[i] = tetrahedron_center[i] - normal_vector[i]*Si_O_distance;
//    }

//    // Using Mathilde's functions
//    Atom silicon = atoms[0];
//    vector<Atom> existing_oxygen(1, atoms[1]);
//    vector<Atom> missing_oxygen_atoms = find_three_missing_oxygen_atoms(silicon, existing_oxygen);

//    for (int i = 0; i < missing_oxygen_atoms.size(); i++) {
//        atoms.push_back(missing_oxygen_atoms[i]);
//    }
//    return atoms;
//}

void Creator::add_tetrahedra() {
    vector<double> system_size = mts0_io->get_lx_ly_lz();

    vector<double> tetrahedron_center(3);
    tetrahedron_center[0] = system_size[0]/4.0;
    tetrahedron_center[1] = system_size[1]/4.0;
    tetrahedron_center[2] = system_size[2]/4.0;

    vector<double> normal_vector(3);
//    normal_vector[0] = 1.1;
//    normal_vector[1] = -0.67;
//    normal_vector[2] = 1.9;
    normal_vector[0] = 1.0;
    normal_vector[1] = 1.0;
    normal_vector[2] = 1.0;
    normal_vector = normalize_vector(normal_vector);

    vector<Atom> tetrahedron_atoms = make_tetrahedron(tetrahedron_center, normal_vector, SI_O_DISTANCE);
//    vector<Atom> tetrahedron_atoms = make_tetrahedron_mathildes_function(tetrahedron_center, normal_vector, SI_O_DISTANCE);
    for (int i = 0; i < tetrahedron_atoms.size(); i++) {
        ensure_within_periodic_boundaries(tetrahedron_atoms[i].position);
    }
    add_atoms(tetrahedron_atoms);
}

vector<Atom> Creator::make_tetrahedron(const vector<double> &tetrahedron_center, const vector<double> &normal_vector, double Si_O_distance) {
    vector<Atom> atoms(5); // one Silicon and four Oxygen

    // Silicon
    atoms[0].type = SI_TYPE;
    atoms[0].position = tetrahedron_center;

    // First oxygen
    atoms[1].type = A_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[1].position[i] = tetrahedron_center[i] - normal_vector[i]*Si_O_distance;
    }

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

    // Unit vector from point in plane, normal to vector from point in plane to second oxygen
    vector<double> unit_vector_in_plane = vector_cross_product(point_in_plane_to_second_oxygen, normal_vector);
    unit_vector_in_plane = normalize_vector(unit_vector_in_plane);

    // Unit vector from point in plane, parallel to vector from point in plane to second oxygen
    vector<double> unit_vector_from_point_to_second_oxygen = normalize_vector(point_in_plane_to_second_oxygen);

    double Lx = distance_from_point_in_plane_to_second_oxygen*sin(60.0*PI/180.0);
    double Ly = distance_from_point_in_plane_to_second_oxygen*cos(60.0*PI/180.0);

    // Third and fourth oxygen
    atoms[3].type = O_TYPE;
    atoms[4].type = H_TYPE;
    for (int i = 0; i < 3; i++) {
        atoms[3].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] - unit_vector_from_point_to_second_oxygen[i]*Ly + unit_vector_in_plane[i]*Lx;
        atoms[4].position[i] = tetrahedron_center[i] + silicon_to_point_in_plane[i] - unit_vector_from_point_to_second_oxygen[i]*Ly - unit_vector_in_plane[i]*Lx;
    }

    return atoms;
}

void Creator::remove_all_atoms() {
    for (int atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        mts0_io->remove_atom(atom_index);
    }
    cout << "Removed all atoms in the system." << endl;
}

void Creator::remove_all_atoms_in_half_system() {
    vector<double> system_size = mts0_io->get_lx_ly_lz();
    vector<double> half_system_size(3);
    half_system_size[0] = system_size[0]/2.0;
    half_system_size[1] = system_size[1]/2.0;
    half_system_size[2] = system_size[2]/2.0;

    int n_atoms = mts0_io->get_number_of_atoms();
    int n_removed_atoms = 0;

    for (int atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        if (mts0_io->positions[atom_index][0] > half_system_size[0]) {
            mts0_io->remove_atom(atom_index);
            n_removed_atoms++;
        }
    }

    cout << "remove_all_atoms_in_half_system(): removed " << n_removed_atoms << " atoms (" << setprecision(2) << (n_removed_atoms/double(n_atoms))*100.0 << "%)." << endl;
}

void Creator::cut_system_along_space_diagonal() {
    vector<double> system_size = mts0_io->get_lx_ly_lz();
    vector<double> half_system_size(3);
    half_system_size[0] = system_size[0]/2.0;
    half_system_size[1] = system_size[1]/2.0;
    half_system_size[2] = system_size[2]/2.0;

    int n_atoms = mts0_io->get_number_of_atoms();
    int n_removed_atoms = 0;

    for (int atom_index = 0; atom_index < mts0_io->positions.size(); atom_index++) {
        double position_sum = 0.0;
        for (int i = 0; i < 3; i++) {
            position_sum += mts0_io->positions[atom_index][i];
        }
        if (position_sum > system_size[0]) {
            mts0_io->remove_atom(atom_index);
            n_removed_atoms++;
        }
    }

    cout << "cut_system_along_space_diagonal(): removed " << n_removed_atoms << " atoms (" << setprecision(2) << (n_removed_atoms/double(n_atoms))*100.0 << "%)." << endl;
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
            // do something appropriate
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
        u[0] = 1.0 - (!small_component_in_normal_vector[0])*((normal_vector[1] + normal_vector[2])/normal_vector[0]);
        u[1] = 1.0 - (!small_component_in_normal_vector[1])*((normal_vector[0] + normal_vector[2])/normal_vector[1]);
        u[2] = 1.0 - (!small_component_in_normal_vector[2])*((normal_vector[0] + normal_vector[1])/normal_vector[2]);
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

//vector<double> Creator::calculate_normal_to_plane_from_three_vector_points(const vector<double> &point_1, const vector<double> &point_2, const vector<double> &point_3) {
//    vector<double> vector_in_plane1 = calculate_vector_difference_using_minimum_image_convention(point_1, point_2);
//    vector<double> vector_in_plane2 = calculate_vector_difference_using_minimum_image_convention(point_1, point_3);
//    vector<double> normal_to_plane = vector_cross_product(vector_in_plane1, vector_in_plane2);
//    double inverse_length = 1.0/sqrt(normal_to_plane[0]*normal_to_plane[0] + normal_to_plane[1]*normal_to_plane[1] + normal_to_plane[2]*normal_to_plane[2]);
//    normal_to_plane[0] *= inverse_length;
//    normal_to_plane[1] *= inverse_length;
//    normal_to_plane[2] *= inverse_length;

//    return normal_to_plane;
//}

//Atom Creator::find_one_missing_oxygen_atom(Atom silicon, vector<Atom> existing_oxygen) {
//    vector<double> missing_oxygen_position(3);
//    Atom missing_oxygen;
//    vector<double> normal_to_plane = calculate_normal_to_plane_from_three_vector_points(existing_oxygen[0].position, existing_oxygen[1].position, existing_oxygen[2].position);
//    normal_to_plane = normalize_vector(normal_to_plane);

//    vector<double> vector_from_oxygen_to_silicon = calculate_vector_difference_using_minimum_image_convention(silicon.position, existing_oxygen[0].position);
//    int sign = vector_dot_product(normal_to_plane, vector_from_oxygen_to_silicon) > 0 ? 1 : -1;
//    vector<double> insertion_position(3);
//    //find position to atom to be inserted
//    //insertion_positions[particle][0][0] = types[particle];
//    missing_oxygen_position[0] = sign*normal_to_plane[0]*SI_O_DISTANCE + silicon.position[0];
//    missing_oxygen_position[1] = sign*normal_to_plane[1]*SI_O_DISTANCE + silicon.position[1];
//    missing_oxygen_position[2] = sign*normal_to_plane[2]*SI_O_DISTANCE + silicon.position[2];

//    ensure_within_periodic_boundaries(missing_oxygen_position);
//    missing_oxygen.type = A_TYPE;
//    missing_oxygen.position = missing_oxygen_position;
//    return missing_oxygen;
//}

//vector<Atom> Creator::find_two_missing_oxygen_atoms(Atom silicon, vector<Atom> existing_oxygen) {
//    //find one of the missing
//    vector<double> missing_oxygen_position(3);
//    vector<Atom> missing_oxygen(2);
//    //find the normal vector to the plane of silicon and the two existing oxygen
//    vector<double> normal_to_plane = calculate_normal_to_plane_from_three_vector_points(silicon.position, existing_oxygen[0].position, existing_oxygen[1].position);
//    //find the vector between oxygen vectors, in the plane
//    vector<double> vec_between_silicon_and_oxygen_1 = calculate_vector_difference_using_minimum_image_convention(existing_oxygen[0].position, silicon.position);
//    vector<double> vec_between_silicon_and_oxygen_2 = calculate_vector_difference_using_minimum_image_convention(existing_oxygen[1].position, silicon.position);

//    vector<double> vec_sum_existing_oxygen(3);
//    vec_sum_existing_oxygen[0] = vec_between_silicon_and_oxygen_1[0] + vec_between_silicon_and_oxygen_2[0];
//    vec_sum_existing_oxygen[1] = vec_between_silicon_and_oxygen_1[1] + vec_between_silicon_and_oxygen_2[1];
//    vec_sum_existing_oxygen[2] = vec_between_silicon_and_oxygen_1[2] + vec_between_silicon_and_oxygen_2[2];

//    //normalize the two vectors.
//    normal_to_plane = normalize_vector(normal_to_plane);
//    vec_sum_existing_oxygen = normalize_vector(vec_sum_existing_oxygen);

//    //position of the two new vectors
//    missing_oxygen_position[0]=  silicon.position[0] - weight_vec_sum_existing_oxygen*vec_sum_existing_oxygen[0] + weight_normal*normal_to_plane[0];
//    missing_oxygen_position[1]=  silicon.position[1] - weight_vec_sum_existing_oxygen*vec_sum_existing_oxygen[1] + weight_normal*normal_to_plane[1];
//    missing_oxygen_position[2]=  silicon.position[2] - weight_vec_sum_existing_oxygen*vec_sum_existing_oxygen[2] + weight_normal*normal_to_plane[2];
//    missing_oxygen[0].position = missing_oxygen_position;
//    //missing_oxygen[0].atom_type = A_TYPE;
//    missing_oxygen[0].type = O_TYPE; //debug

//    missing_oxygen_position[0]=  silicon.position[0] - weight_vec_sum_existing_oxygen*vec_sum_existing_oxygen[0] - weight_normal*normal_to_plane[0];
//    missing_oxygen_position[1]=  silicon.position[1] - weight_vec_sum_existing_oxygen*vec_sum_existing_oxygen[1] - weight_normal*normal_to_plane[1];
//    missing_oxygen_position[2]=  silicon.position[2] - weight_vec_sum_existing_oxygen*vec_sum_existing_oxygen[2] - weight_normal*normal_to_plane[2];
//    missing_oxygen[1].position = missing_oxygen_position;
//    //missing_oxygen[1].atom_type = A_TYPE;
//    missing_oxygen[1].type = O_TYPE;//debug

//    ensure_within_periodic_boundaries(missing_oxygen_position);

//    return missing_oxygen;
//}

//vector<Atom> Creator::find_three_missing_oxygen_atoms(Atom silicon, vector<Atom> existing_oxygen) {
//    //find one of the missing
//    vector<double> missing_oxygen_position(3);
//    //vector<Atom>existing_oxygen_vec;
//    //existing_oxygen_vec.push_back(existing_oxygen);
//    vector<Atom> missing_oxygen(3);
//    //find vector betweeen silicon and oxygen and use as normal vector to plane. (n)
//    vector<double> vec_between_silicon_oxygen = calculate_vector_difference_using_minimum_image_convention(existing_oxygen[0].position, silicon.position);
//    //normalize
//    double inverse_length_vec_between_silicon_oxygen = 1/sqrt(vec_between_silicon_oxygen[0]*vec_between_silicon_oxygen[0] + vec_between_silicon_oxygen[1]*vec_between_silicon_oxygen[1] + vec_between_silicon_oxygen[2]*vec_between_silicon_oxygen[2]);
//    vec_between_silicon_oxygen[0] *= inverse_length_vec_between_silicon_oxygen;
//    vec_between_silicon_oxygen[1] *= inverse_length_vec_between_silicon_oxygen;
//    vec_between_silicon_oxygen[2] *= inverse_length_vec_between_silicon_oxygen;


//    //find vector between silicon and point in plane using
//    //normalize and multiply with correct length calculated earlier. (p)

//    //find vector between point in plane and the the position of the oxygen atom to be inserted by using nx*rx + ny*ry + nz*rz = 0
//    //choose rx = ry = 1
//    //rz = -(nx+ny)/nz
//    vector<double> vec_between_point_in_plane_new_pos(3);
//    vec_between_point_in_plane_new_pos[0] = 1;
//    vec_between_point_in_plane_new_pos[1] = 1;
//    vec_between_point_in_plane_new_pos[2] = -(vec_between_silicon_oxygen[0] + vec_between_silicon_oxygen[1])/vec_between_silicon_oxygen[2];
//    //normalize and multiply with correct length calculated earlier.
//    double inverse_length_vec_between_point_in_plane_new_pos = 1/sqrt(vec_between_point_in_plane_new_pos[0]*vec_between_point_in_plane_new_pos[0] + vec_between_point_in_plane_new_pos[1]*vec_between_point_in_plane_new_pos[1]+ vec_between_point_in_plane_new_pos[2]*vec_between_point_in_plane_new_pos[2]);
//    vec_between_point_in_plane_new_pos[0] *= inverse_length_vec_between_point_in_plane_new_pos*length_from_point_in_plane_to_new_oxygen;
//    vec_between_point_in_plane_new_pos[1] *= inverse_length_vec_between_point_in_plane_new_pos*length_from_point_in_plane_to_new_oxygen;
//    vec_between_point_in_plane_new_pos[2] *= inverse_length_vec_between_point_in_plane_new_pos*length_from_point_in_plane_to_new_oxygen;

//    //find the position of the new oxygen atom by si+|p|n+r
//    missing_oxygen_position[0] = silicon.position[0] - length_from_silicon_to_point_in_plane*vec_between_silicon_oxygen[0] + vec_between_point_in_plane_new_pos[0];
//    missing_oxygen_position[1] = silicon.position[1] - length_from_silicon_to_point_in_plane*vec_between_silicon_oxygen[1] + vec_between_point_in_plane_new_pos[1];
//    missing_oxygen_position[2] = silicon.position[2] - length_from_silicon_to_point_in_plane*vec_between_silicon_oxygen[2] + vec_between_point_in_plane_new_pos[2];
//    ensure_within_periodic_boundaries(missing_oxygen_position);
//    //call the function that finds two missing oxygen atoms
//    missing_oxygen[0].type = H_TYPE;
//    missing_oxygen[0].position = missing_oxygen_position;
//    //now we only lack two oxygen atoms

//    existing_oxygen.push_back(missing_oxygen[0]);
//    vector<Atom> two_missing_oxygen = find_two_missing_oxygen_atoms(silicon,existing_oxygen);
//    missing_oxygen[1] = two_missing_oxygen[0];
//    missing_oxygen[2] = two_missing_oxygen[1];
//    return missing_oxygen;
//}
