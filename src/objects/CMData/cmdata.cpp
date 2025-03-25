#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include "cmdata.hpp"

CMData::CMData() {
    points = vector<Vec>();
    w_points = vector<float>();
    values = vector<complex<Vec>>();
    dimension = 0;
    is_complex = false;
    is_vector = false;
    with_w = false;
    with_n = false;

    filled = false;
}

CMData::CMData(vector<Vec> points, vector<complex<Vec>> values, int dimension, bool with_w, bool with_n, bool is_complex, bool is_vector) {
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    this->is_vector = is_vector;
    this->with_w = with_w;
    this->with_n = with_n;

    filled = true;
}

CMData::CMData(string filename) {
    *this = load(filename);
}

CMData load(string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Failed to open the file.");
    }
    string line;
    int ind = 0;
    int dimension = 0;
    bool w_col = false;
    bool n_col = false;
    bool is_complex = false;
    bool is_vector = false;
    vector<Vec> points;
    vector<complex<Vec>> values;
    while (getline(file, line)) {
        istringstream iss(line);
        if (ind == 0) {
            string word;
            while (iss >> word) {
                if (word == "x") dimension++;
                if (word == "y") dimension++;
                if (word == "z") dimension++;
                if (word == "w") w_col = true;
                if (word == "n") n_col = true;
                if (word == "Re(f)") is_complex = true;
                if (word == "fx") is_vector = true;
                if (word == "Re(fx)") {is_vector = true; is_complex = true;}
            }
            dimension += w_col;
            ind++;
            continue;
        }
        Vec point;
        float value;
        float ivalue;
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;
            point(i) = temp;
        }
        iss >> value;
        point.dimension = dimension;
        points.push_back(point);
        complex<Vec> vval = complex<Vec>(value, 0);
        if (is_complex) {
            iss >> ivalue;
            vval = complex<Vec>(value, ivalue);
        }
        values.push_back(vval);
    }
    return CMData(points, values, dimension, w_col, n_col, is_complex, is_vector);
}

string get_header(int dimension, bool with_w, bool is_complex, bool is_vector) {
   // Determine fheader
    string fheader = "    f    ";
    if (is_complex) {
        fheader = "   Re(f)         Im(f)";
    }
    if (is_vector && !is_complex) {
        fheader = "   fx    ";
    }
    if (is_vector && is_complex) {
        fheader = "    Re(fx)      Im(fx) ";
    }

    // Build header
    string header = "    x         ";
    if (dimension > 1) {
        header += "    y         ";
        if (is_vector && !is_complex) {
            fheader += "   fy        ";
        } else if (is_vector && is_complex) {
            fheader += "    Re(fy)    Im(fy) ";
        }
    }
    if (dimension > 2) {
        header += "    z        ";
        if (is_vector && !is_complex) {
            fheader += "   fz        ";
        }
        if (is_vector && is_complex) {
            fheader += "     Re(fz)    Im(fz) ";
        }
    }
    if (with_w) {
        header += "   w        ";
    }

    header += fheader;
    return header;
}


void save_to_file(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector) {
    ofstream file(filename);
    dimension -= with_w;
    string header = get_header(dimension, with_w, is_complex, is_vector);
    file << fixed << setprecision(6);
    file << header << endl;
    for (int i = 0; i < points.size(); i++) {
        Vec point = points[i];
        complex<Vec> value = values[i];
        file << point(0) << "      ";
        if (dimension > 1) file << point(1) << "      ";
        if (dimension > 2) file << point(2) << "      ";
        if (with_w) file << point(dimension) << "      ";
        if (is_complex) {
            file << value.real()(0) << "      " << value.imag()(0) << "      ";
            if (is_vector) {
                file << value.real()(1) << "      " << value.imag()(1) << "      ";
                if (dimension > 2) {
                    file << value.real()(2) << "      " << value.imag()(2) << "      ";
                }
            }
        } else {
            file << value.real()(0) << "      ";
            if (is_vector) {
                file << value.real()(1) << "      ";
                if (dimension > 2) {
                    file << value.real()(2) << "      ";
                }
            }
        }
        file << endl;
    }
}

void combine_points_and_w(vector<Vec> &points, vector<float> &w_points) {
    for (int i = 0; i < points.size(); i++) {
        points[i].dimension++;
        points[i](points[i].dimension-1) = w_points[i % w_points.size()];
    }
}

void save(CMData &data, string filename) {
    if (data.with_w) combine_points_and_w(data.points, data.w_points);
    save_to_file(filename, data.points, data.values, data.dimension + data.with_w, data.with_w, data.is_complex, data.is_vector);
}
