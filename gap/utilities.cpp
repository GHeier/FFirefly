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

float round(float number, int decimal_places) {
    return round(number * pow(10, decimal_places)) / pow(10, decimal_places);
}
