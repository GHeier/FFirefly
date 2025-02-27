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


using std::cout;
using std::endl;
using std::sort;
using std::vector;
using std::unordered_map;

void save_vp_dir(float mu, const std::vector<Vec>& layer)
{
    std::string directory = "./fermi_velocity_data/"; // directory data is saved to

    std::string filename_start = "fermi_velocity_data"; // name of the data file without the mu label

    unordered_map<float, vector<vector<vector<float>>>> nothing;

    std::ostringstream vp_data;
    
    vp_data << directory << filename_start << "_mu=" << std::fixed << std::setprecision(10) << mu << ".dat"; // saves data to the directory and adds mu to the end of the filename for plotting

    std::ofstream outfile(vp_data.str());

    for(int i = 0; i < layer.size(); i++)
    {
    	Vec k = layer[i]; // was once vec i believe this was an error
    	outfile << vp(k) << endl;
    }
    
    
    outfile.close();

    std::cout << "Data saved to " << vp_data.str() << std::endl;
}

void save_V_dir(float mu, const std::vector<Vec>& layer)
{
    std::string directory = "./interaction_potential_data/"; // name of the directory this is saved to

    std::string filename_start = "interaction_potential_data"; // name of the start of the file name without mu (this is appended later)

    unordered_map<float, vector<vector<vector<float>>>> nothing;

    std::ostringstream V_data;
    
    V_data << directory << filename_start << "_mu=" << std::fixed << std::setprecision(3) << mu << ".dat";

    std::ofstream outfile(V_data.str());

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

    std::cout << "Data saved to " << V_data.str() << std::endl;
}


