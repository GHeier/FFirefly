module Config_module
    use iso_c_binding
    implicit none

    integer :: n, m, l, w_pts, dim, num_eigenvalues_to_save
    real :: max_freq, mu, k_max, t, tn, tnn, U, wc
    character(len=50) :: file
