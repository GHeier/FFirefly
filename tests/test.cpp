#include <unordered_map>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    unordered_map<double, int> map;
    for (int i = 0; i < 3; i++) {
        double point1 = 1.0 * i / 3;
        for (int j = 0; j < 3; j++) {
            double point2 = 1.0 * j / 3;
            map[point1 + point2] = i * 3 + j;
        }
   }
}
