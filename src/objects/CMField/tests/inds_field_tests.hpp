#ifndef INDS_FIELD_TESTS_HPP
#define INDS_FIELD_TESTS_HPP

/**
 * Test declarations for comprehensive inds field testing
 */

bool inds_tests();
void run_all_inds_tests();

// Individual test functions
bool test_scalar_1d_w();
bool test_scalar_2d_w();
bool test_scalar_save_load();
bool test_vector_field();
bool test_single_element_vector();
bool test_uniform_matrix();
bool test_nonuniform_matrix();
bool test_matrix_save_load();
bool test_nonuniform_matrix_save_load();
bool test_uniform_tensor3();
bool test_nonuniform_tensor3();
bool test_singleband_vertex();
bool test_twoband_vertex();
bool test_nonuniform_tensor4();
bool test_vertex_save_load();
bool test_inds_hdf5_storage();
bool test_large_tensor();
bool test_total_index_size();

#endif // INDS_FIELD_TESTS_HPP
