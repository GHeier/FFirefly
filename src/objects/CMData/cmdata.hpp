#pragma once

#include <vector>
#include "../vec.hpp"

using namespace std;

class CMData {
    public:
        vector<Vec> points;
        vector<float> w_points;
        vector<complex<Vec>> values;
        int dimension;
        bool with_w;
        bool with_n;
        vector<int> n_inds;
        bool is_complex;
        bool is_vector;

        bool filled;

        CMData();
        ~CMData();
        CMData(vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool with_n, bool is_complex, bool is_vector);

        CMData(string filename);

        string get_header(int dimension, bool with_w, bool with_n, bool is_complex, bool is_vector);

        void save_hdf5(const std::string& filename) const;
        void load_hdf5(const std::string& filename);
};

CMData load(string filename);
void save_to_file(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool with_n, bool is_complex, bool is_vector);
void save(CMData &data, string filename);

