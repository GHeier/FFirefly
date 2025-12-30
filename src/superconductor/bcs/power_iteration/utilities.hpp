#pragma once

#include <string>
using namespace std;

void progress_bar(float progress, string message = "");
void progress_percent(float progress);
float round_val(float number, int decimal_places);
string get_SC_filename();

