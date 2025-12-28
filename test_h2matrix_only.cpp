/**
 * Run only H2Matrix tests to verify they work
 */
#include <iostream>
#include <cstdio>

// Simple implementation of print_test_results for standalone testing
extern "C" bool print_test_results(bool all_tests[], int num_tests, const char* test_name) {
    int passed = 0;
    for (int i = 0; i < num_tests; i++) {
        if (all_tests[i]) passed++;
    }

    if (passed == num_tests) {
        printf("\033[1;32mAll %d %s passed!\n\033[0m", num_tests, test_name);
        return true;
    } else {
        printf("\033[1;31m - %d/%d %s passed\n\033[0m", passed, num_tests, test_name);
        for (int i = 0; i < num_tests; i++) {
            if (!all_tests[i]) {
                printf("\033[1;31m   - Test %d failed\n\033[0m", i+1);
            }
        }
        return false;
    }
}

#include "src/objects/tests/h2matrix_tests.hpp"

int main() {
    std::cout << "Running H2Matrix tests only..." << std::endl;
    bool result = h2matrix_tests();

    if (result) {
        std::cout << "\n✓ All H2Matrix tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ Some H2Matrix tests FAILED!" << std::endl;
        return 1;
    }
}
