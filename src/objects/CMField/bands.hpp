#include "fields.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class Bands {
    public:
        vector<Field_R> fields;
        bool file_found;

        Bands();
        ~Bands();
        float operator()(int n, Vec k);
};
