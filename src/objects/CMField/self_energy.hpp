#include "fields.hpp"

#include <complex>

using namespace std;

class Self_Energy {
  public:
    CMField field;
    bool file_found;

    Self_Energy();
    //  ~Vertex();
    // float operator()(Vec k, float w = 0, string label1 = "", string label2 =
    // "");
    complex<float> operator()(Vec k, complex<float> w, string label1 = "",
                              string label2 = "");
    complex<float> operator()(Vec k, float w, string label1 = "",
                              string label2 = "");
};

