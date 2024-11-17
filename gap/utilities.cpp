/**
 * @file utilities.cpp
 *
 * @brief Simple utilities, generally for display purposes unrelated to the codebase
 *
 * @author Griffin Heier
 */

#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>

#include "cfg.h"

using namespace std;

void progress_bar(float progress, string message) {
    int barWidth = 70;

    std::cout << message << " [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(2) << progress * 100.0 << " %  \r";
    std::cout.flush();
}

void progress_percent(float progress) {
    std::cout << std::setprecision(8) << progress << "    \r";
    std::cout.flush();
}

float round_val(float number, int decimal_places) {
    return round(number * pow(10, decimal_places)) / pow(10, decimal_places);
}

string get_SC_filename() {
    string bcs = "";
    if (FS_only) bcs = "_FS_only";
    std::ostringstream out;
    out.precision(1);
    out << std::fixed << "data/" + potential_name << dim << "D" 
        << "_mu=" << mu << "_U=" << U << "_wc=" << wc 
        << "_n=" << n << bcs << ".dat";
    string file_name = std::move(out).str();
    return file_name;
}

