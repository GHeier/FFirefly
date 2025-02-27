#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;

double compute_average(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 0.0;
    }

    double sum = 0.0;
    int count = 0;
    double value;

    while (file >> value) {
        sum += value;
        count++;
    }

    return (count > 0) ? (sum / count) : 0.0;
}

std::string extract_mu(const std::string& filename) {
    std::size_t start = filename.find("mu=") + 3;
    std::size_t end = filename.find(".dat", start);
    return filename.substr(start, end - start);
}

int main() {
    std::ofstream output_file("average_fermi_velocity.dat");

    for (const auto& entry : fs::directory_iterator(".")) {
        std::string filename = entry.path().filename().string();

        if (filename.find("fermi_velocity_data_mu=") == 0 && filename.ends_with(".dat")) {
            std::string mu = extract_mu(filename);
            double average = compute_average(filename);

            output_file << mu << "\t" << average << std::endl;
        }
    }

    std::cout << "Averages written to average_fermi_velocity.dat" << std::endl;
    return 0;
}

