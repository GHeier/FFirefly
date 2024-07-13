#include <string>
#pragma once
#ifndef UTILITIES_H_
#define UTILITIES_H_

using namespace std;

void progress_bar(double progress, string message = "");
void progress_percent(double progress);
double round(double number, int decimal_places);

#endif
