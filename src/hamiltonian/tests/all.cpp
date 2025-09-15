#include "all.hpp"
#include "../../config/load/c_config.h"
#include "../../config/load/jl_interface.h"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../../superconductor/cfg.hpp"
#include "../../superconductor/save_data.hpp"
#include <cmath>

using namespace std;

extern "C" void load_cpp_config_wrapper() { load_cpp_config(); }

extern "C" bool hamiltonian_tests() {
  printf("\nRunning Hamiltonian tests\n");
  load_cpp_cfg();
  int num_tests = 2;
  bool all_tests[num_tests] = {DOS_test_2d(), DOS_test_3d()};

  return print_test_results(all_tests, num_tests, "Hamiltonian tests");
}

bool DOS_test_3d() {
  // Set global variables for testing
  set_global(nbnd, 1);
  set_global(dimension, 3);
  band.push_back("fermi_gas");
  eff_mass.push_back(0.5);
  set_global(fermi_energy, 1.0);
  set_global(k_mesh, {18, 18, 18});

  float answer = 0.02533;
  printv("Starting 3D Density of States Test\n");

  vector<Vec> FS = get_FS(fermi_energy);
  float DOS = get_DOS(FS);
  printv("Fermi Energy: %.3f\n", fermi_energy);

  if (fabs(answer - DOS) < 0.001) {
    printv("3D Density of States Test Passed\n");
    return true;
  } else {
    printf(
        "3D Density of States Test Failed\nCalculated: %.3f\nExpected: %.3f\n",
        DOS, answer);
    return false;
  }
}

bool DOS_test_2d() {
  // Set global variables for testing
  set_global(nbnd, 1);
  set_global(dimension, 2);
  band.push_back("fermi_gas");
  eff_mass.push_back(0.5);
  set_global(fermi_energy, 1.0);
  set_global(k_mesh, {16, 16, 16});

  float answer = 0.079577;
  printv("Starting 2D Density of States Test\n");

  vector<Vec> FS = get_FS(fermi_energy);
  float DOS = get_DOS(FS);
  printv("Fermi Energy: %.3f\n", fermi_energy);

  if (fabs(answer - DOS) < 0.001) {
    printv("2D Density of States Test Passed\n");
    return true;
  } else {
    printf(
        "2D Density of States Test Failed\nCalculated: %.3f\nExpected: %.3f\n",
        DOS, answer);
    return false;
  }
}
