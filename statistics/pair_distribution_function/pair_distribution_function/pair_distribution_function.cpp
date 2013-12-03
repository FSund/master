/*
pair_distribution_function.cpp

Created by Jørgen Trømborg and Anders Hafreager on 24.04.13.
Copyright (c) 2013 Universitetet i Oslo. All rights reserved.

You have to choose time frames that are not correlated.

This program requires that you have timeframes saved as folder_name_base/000xxx.

nx, ny, nz                   - number of nodes in each dimension
out_file                 - output file (full path)
folder_name_base         - state-file directory (usually base_code/dump/)
folder_name_iterator_start - first timeframe
step                         - timeframe step
stop                         - last timeframe
num_bins                     - number of bins, g(r)
r_max                    - maximum r in g(r)
<atom_type_1 atom_type_2>  - tag pairs (implemented atom types are Si, A, H, O, Na, Cl, X)
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <mts0_io.h>
#include <time.h>
#include <stdlib.h> /* srand, rand */
#include <map>
#include <math.h>
#include <sstream>

using namespace std;

#define NUM_TYPES 8 // 0 isn't used, see list below

#define SI_TYPE 1
#define A_TYPE 2
#define H_TYPE 3
#define O_TYPE 4
#define NA_TYPE 5
#define CL_TYPE 6
#define X_TYPE 7

vector<vector<int> > get_list_of_list_of_atoms_of_each_type(Mts0_io *mts0_io) {
    vector<vector<int> > list_of_list_of_atoms_of_each_type(NUM_TYPES);

    for(int n=0; n<mts0_io->get_number_of_atoms(); n++) {
        int atom_type = mts0_io->atom_types[n];
        list_of_list_of_atoms_of_each_type[atom_type].push_back(n);
    }

    return list_of_list_of_atoms_of_each_type;
}

void calculate_pair_distribution_functions_from_timeframe(Mts0_io *mts0_io, const int &num_bins, const double &r_max, const vector<int> &tags, vector<vector<double> > &pair_distribution_functions, vector<int> &pair_distribution_functions_number_of_entries, vector<vector<int> > &list_of_list_of_atoms_of_each_type) {
    vector<vector<int> > neighbor_list = mts0_io->create_neighbor_list(r_max);
    vector<double> system_size = mts0_io->get_lx_ly_lz();
    double Lx = system_size[0];
    double Ly = system_size[1];
    double Lz = system_size[2];

    int number_of_tag_pairs = tags.size()/2;

    // Loop through the different atom pairs
    for(int tag_pair_index=0;tag_pair_index<number_of_tag_pairs;tag_pair_index++) {
        int wanted_atom_type0 = tags[2*tag_pair_index+0];
        int wanted_atom_type1 = tags[2*tag_pair_index+1];

        for(int i=0;i<list_of_list_of_atoms_of_each_type[wanted_atom_type0].size(); i++) {
            int atom_index0 = list_of_list_of_atoms_of_each_type[wanted_atom_type0][i];
            vector<int> this_atoms_neighbor_list = neighbor_list[atom_index0];

            for(int j=0; j<this_atoms_neighbor_list.size(); j++) {
                int atom_index1 = this_atoms_neighbor_list[j];
                int atom_type1 = mts0_io->atom_types[atom_index1];

                if(atom_type1 == wanted_atom_type1) {
                    double dx = mts0_io->positions[atom_index0][0] - mts0_io->positions[atom_index1][0];
                    double dy = mts0_io->positions[atom_index0][1] - mts0_io->positions[atom_index1][1];
                    double dz = mts0_io->positions[atom_index0][2] - mts0_io->positions[atom_index1][2];

                    if(dx > Lx/2) dx -= Lx; if(dx < -Lx/2) dx += Lx;
                    if(dy > Ly/2) dy -= Ly; if(dy < -Ly/2) dy += Ly;
                    if(dz > Lz/2) dz -= Lz; if(dz < -Lz/2) dz += Lz;

                    double dr = sqrt(dx*dx + dy*dy + dz*dz);
                    int bin_index = dr/r_max*num_bins;
                    // cout << bin_index << endl;
                    if(bin_index < num_bins) {
                        pair_distribution_functions[tag_pair_index][bin_index]++;
                        pair_distribution_functions_number_of_entries[tag_pair_index]++;
                    }
                }
            }
        }
    }

    neighbor_list.clear();
}

vector<vector<double> > calculate_pair_distribution_functions(Mts0_io *mts0_io, const int &num_bins, const double &bin_size, const double &r_max, string folder_name_base, const vector<int> &folder_name_numbers, const vector<int> &tags, vector<string> column_headers) {
    int number_of_tag_pairs     = tags.size()/2;
    int n_time_frames_available = folder_name_numbers.size();
    vector<vector<int> > list_of_list_of_atoms_of_each_type;
    vector<vector<double> > pair_distribution_functions(number_of_tag_pairs);
    vector<int>             pair_distribution_functions_number_of_entries(number_of_tag_pairs,0);
    char *file_path = new char[5000];

    for(int i=0;i<number_of_tag_pairs;i++) pair_distribution_functions[i].resize(num_bins,0);

    for (int i=0; i<n_time_frames_available; i++) {
        sprintf(file_path,"%s%06d/mts0/",folder_name_base.c_str(), folder_name_numbers[i]);
        cout << "Working on folder " << file_path << endl;
        mts0_io->load_atoms(file_path);
        list_of_list_of_atoms_of_each_type = get_list_of_list_of_atoms_of_each_type(mts0_io);
        calculate_pair_distribution_functions_from_timeframe(mts0_io, num_bins, r_max, tags, pair_distribution_functions, pair_distribution_functions_number_of_entries, list_of_list_of_atoms_of_each_type);
    }

    // Normalize
    for(int tag_pair_index=0; tag_pair_index<number_of_tag_pairs; tag_pair_index++) {
        cout << "Count " << column_headers[tag_pair_index+1] << ": " << pair_distribution_functions_number_of_entries[tag_pair_index] << endl;
        int atom_type_0 = tags[2*tag_pair_index+0];
        int atom_type_1 = tags[2*tag_pair_index+1];

        for(int bin_index=0; bin_index<num_bins; bin_index++) {
            double inner_radius = bin_size*bin_index;
            double outer_radius = bin_size*(bin_index+1);
            double shell_volume = 4.0/3*M_PI*(pow(outer_radius,3) - pow(inner_radius,3));

            // pair_distribution_functions[tag_pair_index][bin_index] /= pair_distribution_functions_number_of_entries[tag_pair_index]*shell_volume;
            int number_of_atom_type_0 = list_of_list_of_atoms_of_each_type[atom_type_0].size();
            double number_density_of_atom_type_1 = list_of_list_of_atoms_of_each_type[atom_type_1].size() / mts0_io->get_volume();

            pair_distribution_functions[tag_pair_index][bin_index] /= n_time_frames_available*number_of_atom_type_0*shell_volume*number_density_of_atom_type_1;
        }
    }

    pair_distribution_functions_number_of_entries.clear();
    list_of_list_of_atoms_of_each_type.clear();
    delete file_path;

    return pair_distribution_functions;
}

void write_pair_distribution_functions_to_file(vector<double> &bin_centers, vector<vector<double> > &pair_distribution_functions, string out_file, vector<string> &column_headers) {
    FILE *file = fopen(out_file.c_str(),"w");

    if(file==NULL) {
        cout << "Error in write_pair_distribution_functions_to_file(): Could not open file " << out_file << endl;
        exit(1);
    }
    int number_of_tag_pairs = pair_distribution_functions.size();

    for(int column=0; column<column_headers.size(); column++) {
        fprintf(file,"%15s",column_headers[column].c_str());
    }
    fprintf(file,"\n");

    for(int row=0; row<bin_centers.size(); row++) {
        fprintf(file, "%12.4e",bin_centers[row]);
        for(int tag_pair_index=0; tag_pair_index<number_of_tag_pairs; tag_pair_index++) {
            fprintf(file, "%12.4e",pair_distribution_functions[tag_pair_index][row]);
        }
        fprintf(file,"\n");
    }

    fclose(file);
}

int main(int args, char *argv[]) {
    int index_of_first_tag = 11;
    map<string,int> tag_to_int_map;

    tag_to_int_map["Si"] = SI_TYPE;
    tag_to_int_map["A"]  = A_TYPE;
    tag_to_int_map["H"]  = H_TYPE;
    tag_to_int_map["O"]  = O_TYPE;
    tag_to_int_map["Na"] = NA_TYPE;
    tag_to_int_map["Cl"] = CL_TYPE;
    tag_to_int_map["X"]  = X_TYPE;

    if (args<index_of_first_tag+1){
        cout << "Run program with " << endl;
        cout << "./pair_distribution_function.x nx ny nz out_file folder_name_base folder_name_iterator_start step stop num_bins r_max <atom_type_1 atom_type_2> <atom_type_1 atom_type_2> ..." << endl;
        exit(1);
    }
    int nx = atoi(argv[1]); int ny = atoi(argv[2]); int nz = atoi(argv[3]);

    string out_file                = argv[4];
    string folder_name_base        = argv[5];
    int folder_name_iterator_start = atoi(argv[6]);
    int folder_name_iterator_step  = atoi(argv[7]);
    int folder_name_iterator_stop  = atoi(argv[8]);
    int num_bins                   = atoi(argv[9]);
    double r_max                   = atof(argv[10]);
    
    vector<int> tags;
    for (int i = index_of_first_tag; i<args; i++) {
        if (tag_to_int_map.count(argv[i]) <= 0) {
            cout << "Error, unknown atom type: " << argv[i] << endl;
            return 1;
        }
        tags.push_back(tag_to_int_map[argv[i]]); 
    }

    if(tags.size() % 2 != 0) { cout << "Input error, tags must come in pairs." << endl; exit(1); }
    int number_of_tag_pairs = tags.size()/2;
    
    // In general, folder names can increment with steps != 1
    int n_time_frames_available = (folder_name_iterator_stop-folder_name_iterator_start)/folder_name_iterator_step+1;
    vector<int> folder_name_numbers(n_time_frames_available);
    for (int i=0; i<n_time_frames_available; i++) {
        folder_name_numbers[i] = folder_name_iterator_start+i*folder_name_iterator_step;
    }
    if (folder_name_numbers[n_time_frames_available-1] != folder_name_iterator_stop){
        std::cout << "Wrong implementation or input of folder_name_iterator_start-step-stop. Expected start+n*step=stop. Aborting." << std::endl;
        exit(1);
    }

    vector<string> column_headers;
    column_headers.push_back("bin_ctrs[A]");
    for(int tag_pair_index=0; tag_pair_index<number_of_tag_pairs; tag_pair_index++) {
        stringstream tmp;
        string atom0 = argv[index_of_first_tag + 2*tag_pair_index];
        string atom1 = argv[index_of_first_tag + 2*tag_pair_index+1];
        tmp << "g(" << atom0 << "-" << atom1 << ",r)";
        column_headers.push_back(tmp.str());
    }

    double bin_size = r_max / num_bins;
    vector<double> bin_centers(num_bins);
    for(int bin_index=0; bin_index<num_bins; bin_index++) {
        bin_centers[bin_index] = bin_index*bin_size + bin_size/2;
    }

    Mts0_io *mts0_io = new Mts0_io(nx,ny,nz);
    vector<vector<double> > pair_distribution_functions = calculate_pair_distribution_functions(mts0_io, num_bins, bin_size, r_max, folder_name_base, folder_name_numbers, tags, column_headers);


    write_pair_distribution_functions_to_file(bin_centers, pair_distribution_functions, out_file, column_headers);
    
    return 0;
}
