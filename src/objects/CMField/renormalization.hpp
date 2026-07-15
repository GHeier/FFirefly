#include "fields.hpp"

using namespace std;

class Renormalization {
  public:
    Field_R field;
    bool file_found;

    Renormalization();
    float operator()(Vec k, string label1 = "",
                              string label2 = "");
};
