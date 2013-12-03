/*
pair_distribution_function.cpp

Based on program created by Jørgen Trømborg and Anders Hafreager on 24.04.13. Copyright (c) 2013 Universitetet i Oslo. All rights reserved.
Modified by Filip Sund.

You have to choose time frames that are not correlated.

This program requires that you have timeframes saved as folder_name_base/000xxx.

nx, ny, nz      			 - number of nodes in each dimension
out_file   				 - output file (full path)
folder_name_base  		 - state-file directory (usually base_code/dump/)
folder_name_iterator_start - first timeframe
step 						 - timeframe step
stop 						 - last timeframe
num_bins 					 - number of bins, g(r)
r_max 					 - maximum r in g(r)
<atom_type_1 atom_type_2>  - tag pairs (implemented atom types are Si, A, H, O, Na, Cl, X)
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <mts0_io.h>
// #include <time.h>
#include <cstdlib> /* srand, rand */
#include <map>
#include <cmath>
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

typedef unsigned int uint;

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

	if (args < index_of_first_tag+1){
		cout << "Run program with " << endl;
		cout << "./pair_distribution_function.x  nx  ny  nz  out_file  folder_name_base  folder_name_iterator_start  step  stop  num_bins  r_max  <atom_type_1 atom_type_2>  <atom_type_1 atom_type_2>  ..." << endl;
		exit(1);
	}
	uint nx = atoi(argv[1]);
    uint ny = atoi(argv[2]);
    uint nz = atoi(argv[3]);

	string out_file          	    = argv[4];
	string folder_name_base 	    = argv[5];
    uint folder_name_iterator_start = atoi(argv[6]);
    uint folder_name_iterator_step  = atoi(argv[7]);
    uint folder_name_iterator_stop  = atoi(argv[8]);
    uint num_bins                   = atoi(argv[9]);
    double r_max                    = atof(argv[10]);
	
	vector<uint> tags;
    for (uint i = index_of_first_tag; i<args; i++) {
        if (tag_to_int_map.count(argv[i]) <= 0) {
            cout << "Error, unknown atom type: " << argv[i] << ", aborting!" << endl;
            return 1;
        }
        tags.push_back(tag_to_int_map[argv[i]]); 
    }

	if (tags.size() % 2 != 0) {
        cout << "Input error, tags must come in pairs." << endl; exit(1);
    }
	uint number_of_tag_pairs = tags.size()/2;
    
    // In general, folder names can increment with steps != 1
    uint n_time_frames_available = (folder_name_iterator_stop - folder_name_iterator_start)/folder_name_iterator_step + 1;
    vector<uint> folder_name_numbers(n_time_frames_available);
    for (uint i = 0; i < n_time_frames_available; i++) {
        folder_name_numbers[i] = folder_name_iterator_start + i*folder_name_iterator_step;
    }
    if (folder_name_numbers[n_time_frames_available-1] != folder_name_iterator_stop) {
        std::cout << "Wrong implementation or input of folder_name_iterator_start-step-stop. Expected start+n*step=stop. Aborting." << std::endl;
        return 1;
    }

    vector<string> column_headers;
	column_headers.push_back("bin_ctrs[A]");
	for (uint tag_pair_index = 0; tag_pair_index < number_of_tag_pairs; tag_pair_index++) {
		stringstream tmp;
		string atom0 = argv[index_of_first_tag + 2*tag_pair_index];
		string atom1 = argv[index_of_first_tag + 2*tag_pair_index+1];
		tmp << "g(" << atom0 << "-" << atom1 << ",r)";
		column_headers.push_back(tmp.str());
	}

	double bin_size = r_max/num_bins;
	vector<double> bin_centers(num_bins);
	for (uint bin_index = 0; bin_index < num_bins; bin_index++) {
		bin_centers[bin_index] = bin_index*bin_size + bin_size/2;
	}

	Mts0_io *mts0_io = new Mts0_io(nx,ny,nz);
    uint n_slices = 1;
	vector<vector<vector<double> > > pair_distribution_functions = calculate_pair_distribution_functions(
        mts0_io, n_slices, num_bins, bin_size, r_max, folder_name_base, folder_name_numbers, tags, column_headers);

    // uint slice = 0;
	// write_pair_distribution_functions_to_file(bin_centers, pair_distribution_functions[slice], out_file, column_headers);
	
    return 0;
}
