module confighub
    use iso_c_binding
    implicit none
    ! Global variables

    ![CONTROL]
    character(kind=c_char), bind(C, name="c_category") :: c_category(50)
    character(len=50) :: category
    character(kind=c_char), bind(C, name="c_calculation") :: c_calculation(50)
    character(len=50) :: calculation
    character(kind=c_char), bind(C, name="c_outdir") :: c_outdir(50)
    character(len=50) :: outdir
    character(kind=c_char), bind(C, name="c_prefix") :: c_prefix(50)
    character(len=50) :: prefix
    character(kind=c_char), bind(C, name="c_verbosity") :: c_verbosity(50)
    character(len=50) :: verbosity
    character(kind=c_char), bind(C, name="c_datfile_in") :: c_datfile_in(50)
    character(len=50) :: datfile_in
    character(kind=c_char), bind(C, name="c_datfile_out") :: c_datfile_out(50)
    character(len=50) :: datfile_out

    ![SYSTEM]
    character(kind=c_char), bind(C, name="c_interaction") :: c_interaction(50)
    character(len=50) :: interaction
    integer(c_int), bind(C, name="c_dimension") :: c_dimension
    integer :: dimension
    integer(c_int), bind(C, name="c_ibrav") :: c_ibrav
    integer :: ibrav
    real(c_float), bind(C, name="c_fermi_energy") :: c_fermi_energy
    real :: fermi_energy
    real(c_float), bind(C, name="c_onsite_U") :: c_onsite_U
    real :: onsite_U

    ![MESH]
    integer(c_int), bind(C, name="c_k_mesh") :: c_k_mesh(3)
    integer :: k_mesh(3)
    integer(c_int), bind(C, name="c_q_mesh") :: c_q_mesh(3)
    integer :: q_mesh(3)
    integer(c_int), bind(C, name="c_w_pts") :: c_w_pts
    integer :: w_pts

    ![CELL]
    real(c_float), bind(C, name="c_cell") :: c_cell(3,3)
    real :: cell(3,3)
    real(c_float), bind(C, name="c_brillouin_zone") :: c_brillouin_zone(3,3)
    real :: brillouin_zone(3,3)

    ![BANDS]
    character(kind=c_char), bind(C, name="c_bands") :: c_bands(50)
    character(len=50) :: bands
    real(c_float), bind(C, name="c_eff_mass") :: c_eff_mass
    real :: eff_mass

    ![SUPERCONDUCTOR]
    logical(c_bool), bind(C, name="c_FS_only") :: c_FS_only
    logical :: FS_only
    real(c_float), bind(C, name="c_bcs_cutoff_frequency") :: c_bcs_cutoff_frequency
    real :: bcs_cutoff_frequency
    integer(c_int), bind(C, name="c_num_eigenvalues_to_save") :: c_num_eigenvalues_to_save
    integer :: num_eigenvalues_to_save
    integer(c_int), bind(C, name="c_frequency_pts") :: c_frequency_pts
    integer :: frequency_pts

    ![RESPONSE]
    logical(c_bool), bind(C, name="c_dynamic") :: c_dynamic
    logical :: dynamic
    ! End of global variables

    interface
    ! Global functions

    ![CONTROL]
        function get_category() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_category
            end function get_category

        function get_calculation() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_calculation
            end function get_calculation

        function get_outdir() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_outdir
            end function get_outdir

        function get_prefix() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_prefix
            end function get_prefix

        function get_verbosity() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_verbosity
            end function get_verbosity

        function get_datfile_in() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_datfile_in
            end function get_datfile_in

        function get_datfile_out() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_datfile_out
            end function get_datfile_out


    ![SYSTEM]
        function get_interaction() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_interaction
            end function get_interaction


    ![MESH]

    ![CELL]

    ![BANDS]
        function get_bands() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_bands
            end function get_bands


    ![SUPERCONDUCTOR]

    ![RESPONSE]
