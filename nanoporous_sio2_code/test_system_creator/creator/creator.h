#include <iostream>
#include <iomanip> // setw, setprecision etc.
#include <fstream> // ofstream
#include <sstream> // ostringstream
#include <string>
#include <vector>
#include "../mts0_io/mts0_io.h"
// #include <time.h>
#include <cstdlib> // srand, rand
#include <cmath> // floor, cos, sin

using namespace std;

#define SI_TYPE 1
#define A_TYPE 2
#define H_TYPE 3
#define O_TYPE 4
#define NA_TYPE 5
#define CL_TYPE 6
#define X_TYPE 7
// // Masses in amu
// #define SI_MASS 28.0855
// #define NA_MASS 22.989769
// #define CL_MASS 35.453
// #define O_MASS 15.9994
// #define H_MASS 1.00794

// char *type[] = {(char*)"Not in use", (char*)"Si",(char*)"A ",(char*)"H ",(char*)"O ",(char*)"Na",(char*)"Cl",(char*)"X "};

const double SI_O_DISTANCE = 1.8;
const double PI = atan(1.0)*4.0;
const double TETRAHEDRAL_ANGLE_RAD = 109.4712*PI/180.0;

struct Atom;

class Creator{
public:
    Creator(string mts0_directory_in, int nx, int ny, int nz);

    void remove_all_atoms();
    void remove_all_atoms_in_half_system();
    void cut_system_along_space_diagonal();
    void add_tetrahedra();

    Mts0_io* mts0_io;

private:
    vector<double> system_size;
    vector<double> half_system_size;

    double weight_normal;
    double weight_vec_between_existing_oxygen;
    double weight_vec_sum_existing_oxygen;
    double length_from_silicon_to_point_in_plane;
    double length_from_point_in_plane_to_new_oxygen;

    vector<Atom> make_tetrahedron(const vector<double> &center, const vector<double> &normal_vector, double Si_O_distance);
    inline void add_atoms(const vector<Atom> &atoms);
    inline void add_atom(const Atom &atom);

    // Vector, periodic boundaries and minimum image convention stuff
    double vector_dot_product(const vector<double> &u, const vector<double> &v);
    vector<double> vector_cross_product(const vector<double> &u, const vector<double> &v);
    inline vector<double> normalize_vector(const vector<double>& u);
    vector<double> find_vector_in_plane(const vector<double>& normal_vector);
    void ensure_within_periodic_boundaries(vector<double> &position);

    // from Mathilde
    vector<Atom> make_tetrahedron_mathildes_function(const vector<double> &tetrahedron_center, const vector<double> &normal_vector, double Si_O_distance);
    Atom find_one_missing_oxygen_atom(Atom silicon, vector<Atom> existing_oxygen);
    vector<Atom> find_two_missing_oxygen_atoms(Atom silicon, vector<Atom> existing_oxygen);
    vector<Atom> find_three_missing_oxygen_atoms(Atom silicon, vector<Atom> existing_oxygen);
    vector<double> calculate_normal_to_plane_from_three_vector_points(const vector<double> &point_1,const vector<double> &point_2,const vector<double> &point_3);
    vector<double> calculate_vector_difference_using_minimum_image_convention(const vector<double> &position1, const vector<double> &position2);
};

struct Atom{
    Atom(){
        type = 0;
        position = vector<double>(3);
        velocity = vector<double>(3,0.0);
    }
    Atom(int type, vector<double> position, vector<double> velocity=vector<double>(3,0.0)):
        type(type),
        position(position),
        velocity(velocity)
    {
    }

    int type;
    vector<double> position;
    vector<double> velocity;
};
