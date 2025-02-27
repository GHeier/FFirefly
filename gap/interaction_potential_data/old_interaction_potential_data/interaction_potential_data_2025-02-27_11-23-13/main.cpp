#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <random>
#include <algorithm>

namespace fs = std::filesystem;

double compute_average(const std::string& filename, int num_samples) //computes average of a random number of samples
{
	printf("Starting average\n");
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 0.0;
    }

    std::vector<double> values;
    double value;

    while (file >> value) {
        values.push_back(value);
    }
    printf("File read\n");

    // If there are fewer than num_samples data points, use all of them
    int count = std::min(static_cast<int>(values.size()), num_samples);

    // shuffle and seletc the first 'count' values randomly
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(values.begin(), values.end(), gen);

    double sum = 0.0;
    count = values.size();

    std::cout << filename << " contains " << count << " points." << std::endl;

    for (int i = 0; i < count; ++i) {
        sum += values[i];
    }

    return (count > 0) ? (sum / count) : 0.0;
}

std::string extract_mu(const std::string& filename) {
    std::size_t start = filename.find("mu=") + 3;
    std::size_t end = filename.find(".dat", start);
    return filename.substr(start, end - start);
}

int main() 
{
    int num_samples = 1000; // samples for average

    std::ofstream output_file("average_interaction_potential.dat");

    for (const auto& entry : fs::directory_iterator(".")) {
        std::string filename = entry.path().filename().string();

        if (filename.find("interaction_potential_data_mu=") == 0 && filename.ends_with(".dat")) {
            std::string mu = extract_mu(filename);
            double average = compute_average(filename, num_samples);
	
	    std::cout << "file: " << filename << std::endl;
	    std::cout << "Average taken for interaction_potential_data_mu=" << mu << "." << std::endl;
	    std::cout << "Average equals " << average << std::endl;


            output_file << mu << "\t" << average << std::endl;
        }
    }

    std::cout << "Averages written to average_interaction_potential.dat" << std::endl;
    return 0;
}

