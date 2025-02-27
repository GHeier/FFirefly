/**
 * @file save_data.cpp
 *
 * @brief File contains different functions to save calculated data. save() is for Gap Function
 *
 * @author Griffin Heier
 */

#include <iomanip>
#include <sstream>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <math.h>
#include <unordered_map>
#include <memory>

#include "calculations.h"
#include "frequency_inclusion.hpp"
#include "potential.h"
#include "susceptibility.h"
#include "vec.h"
#include "cfg.h"
#include "matrix.hpp"
#include "eigenvec.hpp"

/*
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

#include "cfg.h"
#include "fermi_surface.h"
#include "frequency_inclusion.hpp"
#include "calculations.h"
#include "potential.h"
#include "save_data.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "utilities.h"

*/

using std::endl;
using std::cout;
using std::vector;
using std::string;

namespace fs = std::filesystem;

void save(string file_name, float T, vector<Vec> FS, Eigenvector* solutions) {
    std::ofstream file(file_name);
    //file.open(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) 
        file << solutions[i].eigenvalue << " ";
    file << endl;
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec q = FS[i]; 
        if(not q.cartesian) 
            q.to_cartesian();
        file << q; 
        for (unsigned int j = 0; j < num_eigenvalues_to_save; j++) {
            file << solutions[j].eigenvector[i] << " ";
        }
        file << endl;
    }
}

void save_with_freq(string file_name, float T, vector<vector<Vec>> &freq_FS, Eigenvector* solutions) {
    std::ofstream file(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

    int size = matrix_size_from_freq_FS(freq_FS);
    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) {
        file << solutions[i].eigenvalue << " ";
    }
    file << endl;
    int ind = 0;
    printf("Size: %d\n", size);
    printf("Eigenvector size: %d\n", solutions[0].size);
    for (unsigned int i = 0; i < freq_FS.size(); i++) {
        for (unsigned int j = 0; j < freq_FS[i].size(); j++) {
            Vec q = freq_FS[i][j]; 
            if(not q.cartesian) 
                q.to_cartesian();
            file << q; 
            for (unsigned int k = 0; k < num_eigenvalues_to_save; k++) {
                if (ind >= size) {
                    printf("Breaking\n");
                }
                file << solutions[k][ind] << " ";
            }
            ind++;
            file << endl;
        }
    }
}





void save_FS(vector<Vec> FS, int iteration, float mu) {
    // completely unnecessary check for existence of FS directory
    if (!fs::exists("FS")) {
        fs::create_directory("FS");
    }

    // cleaning out fs directory on the first iteration
    if (iteration == 0) {
        for (const auto& entry : fs::directory_iterator("FS")) {
            try {
                fs::remove(entry);
            } catch (const fs::filesystem_error& e) {
                cerr << "Error deleting file: " << entry.path() << " - " << e.what() << endl;
            }
        }
    }

    // naming file. include mu for graphing func later

    int mu_prec = 5; // float precision for mu n filename

    ostringstream filename;
    filename << "FS/FS_prec_" << mu_prec << "_mu_" << fixed << setprecision(mu_prec) << mu << "_iter_" << iteration << ".dat";


    std::ofstream file(filename.str());
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
        return;
    }

    // write data
    for (Vec k : FS) {
        file << k << endl;
    }
    
    file.close();
}


void save_FS(vector<Vec> FS) {
    // completely unnecessary check for existence of FS directory
    if (!fs::exists("FS")) {
        fs::create_directory("FS");
    }

    // cleaning out fs directory on the first iteration
    for (const auto& entry : fs::directory_iterator("FS")) {
            try {
                fs::remove(entry);
            } catch (const fs::filesystem_error& e) {
                cerr << "Error deleting file: " << entry.path() << " - " << e.what() << endl;
            }
        }

    std::ofstream file;
    file.open("FS/FS.dat"); // for now manually delete the FS files on succesive iterations
    
    // writing file
    for (Vec k : FS)
        file << k << endl;
}

void save_potential_vs_q(vector<Vec> &FS, Matrix &P, string filename) {
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            float V = -P(i,j) * pow(vp(k2)/k2.area,0.5) * pow(vp(k1)/k1.area,0.5);
            Vec q = k1 - k2; if (q.cartesian == false) q.to_cartesian();
            file << q << V << endl;
        }
    }
}

void save_chi_vs_q(const vector<vector<vector<float>>> &cube, vector<Vec> &FS, string filename) {
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            Vec q = k1 - k2;
            float chi = calculate_chi_from_cube(cube, q);
            file << q << chi << endl;
        }
    }
}


void fermi_velocity_average(float mu)
{
	std::ifstream vp_file("vp_file.dat");
	if(!vp_file)
	{
		std::cerr << "Error: Could not open fermi velocity file (vp_file.dat) for averaging.\n";
		return;
	}
	double sum = 0.0, average = 0.0;
	int count = 0;
	double value;
	
	while(vp_file >> value)
	{
		sum += value;
		count ++;
	}


	string input_var = "y";
	int error_var = 0;

	if(count==0)
	{
		std::cerr << "Error: fermi velocity file (vp_file.dat) is empty. Set average to 0 (y/n)?\n";
		cin >> input_var;
		error_var += 1;

		if(input_var == "n")
		{
			return;
		}
	}

	if(input_var == "y")
	{
		if(error_var == 1)
		{
			average = 0.0;
		}
		else
		{
			average = sum/count;
		}
	}
	
	
	vp_file.close();
	
	std::cout << "Average: " << average << std::endl;

	std::ofstream vp_avgs("average_fermi_velocities.dat", std::ios::app);
	vp_avgs << average << "\t" << mu << endl;
	vp_avgs.close();
}



void interaction_potential_average(float mu)
{
	std::ifstream V_file("V_file.dat");
	if(!V_file)
	{
		std::cerr << "Error: Could not open interaction potential file (V_file.dat) for averaging.\n";
		return;
	}

	double value;
	int count = 0;
	double sum = 0.0, average = 0.0;

	while(V_file >> value)
	{
		sum += value;
		count++;
	}

	string input_var = "y";
	int error_var = 0;

	if(count==0) 
	{
		std::cerr << "Error: interaction potential file (V_file.dat) is empty. Set average to 0 (y/n)?\n";
		cin >> input_var;
		error_var += 1;

		if(input_var == "n")
		{
			return;
		}
	}

	if(input_var == "y")
	{
		if(error_var == 1)
		{
			average = 0.0;
		}
		else
		{
			average = sum/count;
		}
	}

	std::cout << "Average: " << average << std::endl;
	
	V_file.close();

	std::ofstream V_file_avgs("average_interaction_potentials.dat", std::ios::app);
	V_file_avgs << average << "\t" << mu << endl;
	V_file_avgs.close();
}


/*
===============

DATA MANAGEMENT

===============
*/

namespace fs = std::filesystem;

void archive_and_clear_data(const std::string& data_dir)
{
	std::string old_data_dir = data_dir + "/old_" + data_dir;
	if (!fs::exists(old_data_dir))
	{
		fs::create_directories(old_data_dir);
	}

	// timestamp for archive directory. this is done so all archives have unique names
	auto now = std::chrono::system_clock::now();
	auto in_time_t = std::chrono::system_clock::to_time_t(now);
	std::ostringstream oss;
	oss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
	std::string timestamp = oss.str();

	std::string archive_dir = old_data_dir + "/" + data_dir + "_" + timestamp; // archive directory is old_<fermi_velocity or interaction_potential>_data/<fermi_velocity or interaction_potential>_data_<time>
	fs::create_directories(archive_dir);


	// moving files
	for (const auto& entry : fs::directory_iterator(data_dir))
	{
		if (fs::is_regular_file(entry))
		{
			std::string file_path = entry.path().string();
			std::string file_name = entry.path().filename().string();

			fs::copy(file_path, archive_dir + "/" + file_name, fs::copy_options::overwrite_existing);
			if (entry.path().extension() == ".dat")
			{
				fs::remove(file_path);
			}
		}
	}
	std::cout << "Archived data to " << archive_dir << " and cleared .dat files in " << data_dir << "\n";
}

/*
void save_vp_dir(float mu)
{
    std::string directory = "./fermi_velocity_data";

    std::ostringstream vp_data;

    vp_data << directory << "_mu=" << std::fixed << std::setprecision(3) << mu << ".dat";

    std::ofstream outfile(vp_data.str());

        return;
    }

    for(int i = 0; i < layer.size(); i++)
    {
    	vec k = layer[i];
    	outfilr << vp(k) << endl;
    }


    outfile.close();

    std::cout << "Data saved to " << filename.str() << std::endl;
}

void save_V_dir(float mu)
{
    std::string directory = "./interaction_potential_data";

    std::ostringstream vp_data;

    vp_data << directory << "_mu=" << std::fixed << std::setprecision(3) << mu << ".dat";

    std::ofstream outfile(vp_data.str());

        return;
    }

    for(int i = 0; i < layer.size(); i++)
    {
	for(int j = 0; j < layer.size(); j++)
	{
	    Vec k1 = layer[i];
	    Vec k2 = layer[j];
	    float V_val = V(k1, k2, 0, 0, nothing);
	    outfile << V_val << endl;
	}
    }


    outfile.close();

    std::cout << "Data saved to " << filename.str() << std::endl;
}

*/
