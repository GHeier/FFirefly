module confighub
    use iso_c_binding
    implicit none
    ! Global variables

![CONTROL]
    character(len=50) :: category
    character(len=50) :: calculation
    character(len=50) :: outdir
    character(len=50) :: prefix
    character(len=50) :: verbosity
    character(len=50) :: datfile_in
    character(len=50) :: datfile_out

![SYSTEM]
    character(len=50) :: interaction
    integer(c_int), bind(C, name="c_dimension") :: c_dimension
    integer :: dimension
    integer(c_int), bind(C, name="c_ibrav") :: c_ibrav
    integer :: ibrav
    integer(c_int), bind(C, name="c_nbnd") :: c_nbnd
    integer :: nbnd
    real(c_float), bind(C, name="c_fermi_energy") :: c_fermi_energy
    real :: fermi_energy
    real(c_float), bind(C, name="c_Temperature") :: c_Temperature
    real :: Temperature
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

![BRILLOUIN_ZONE]
    real(c_float), bind(C, name="c_brillouin_zone") :: c_brillouin_zone(3,3)
    real :: brillouin_zone(3,3)

![BANDS]
    character(len=50) :: band(50,50)
    real(c_float), bind(C, name="c_eff_mass") :: c_eff_mass(50)
    real :: eff_mass(50)
    real(c_float), bind(C, name="c_t0") :: c_t0(50)
    real :: t0(50)
    real(c_float), bind(C, name="c_t1") :: c_t1(50)
    real :: t1(50)
    real(c_float), bind(C, name="c_t2") :: c_t2(50)
    real :: t2(50)
    real(c_float), bind(C, name="c_t3") :: c_t3(50)
    real :: t3(50)
    real(c_float), bind(C, name="c_t4") :: c_t4(50)
    real :: t4(50)
    real(c_float), bind(C, name="c_t5") :: c_t5(50)
    real :: t5(50)
    real(c_float), bind(C, name="c_t6") :: c_t6(50)
    real :: t6(50)
    real(c_float), bind(C, name="c_t7") :: c_t7(50)
    real :: t7(50)
    real(c_float), bind(C, name="c_t8") :: c_t8(50)
    real :: t8(50)
    real(c_float), bind(C, name="c_t9") :: c_t9(50)
    real :: t9(50)
    real(c_float), bind(C, name="c_t10") :: c_t10(50)
    real :: t10(50)

![SUPERCONDUCTOR]
    character(len=50) :: method
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


![BRILLOUIN_ZONE]


![BANDS]
        function get_band() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_band
    end function get_band













![SUPERCONDUCTOR]
        function get_method() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_method
    end function get_method





![RESPONSE]

    ! End of global functions

    end interface

contains

    function get_string(c_string) result(fortran_string)
        type(c_ptr), intent(in) :: c_string
        character(len=:), allocatable :: fortran_string
        character(kind=c_char, len=1), pointer :: c_string_ptr(:)
        integer :: i, length
        call c_f_pointer(c_string, c_string_ptr, [1000])  ! 1000 is a safe buffer size; adjust as needed
        length = 0
        do i = 1, size(c_string_ptr)
            if (c_string_ptr(i) == c_null_char) exit
            length = length + 1
        end do
        allocate(character(len=length) :: fortran_string)
        do i = 1, length
            fortran_string(i:i) = c_string_ptr(i)
        end do
    end function get_string

    subroutine load_f90_config()
        ! Load variables

![CONTROL]
        category = get_string(get_category())
        calculation = get_string(get_calculation())
        outdir = get_string(get_outdir())
        prefix = get_string(get_prefix())
        verbosity = get_string(get_verbosity())
        datfile_in = get_string(get_datfile_in())
        datfile_out = get_string(get_datfile_out())

![SYSTEM]
        interaction = get_string(get_interaction())
        dimension = c_dimension
        ibrav = c_ibrav
        nbnd = c_nbnd
        fermi_energy = c_fermi_energy
        Temperature = c_Temperature
        onsite_U = c_onsite_U

![MESH]
        k_mesh = c_k_mesh
        q_mesh = c_q_mesh
        w_pts = c_w_pts

![CELL]
        cell = c_cell

![BRILLOUIN_ZONE]
        brillouin_zone = c_brillouin_zone

![BANDS]
        band = get_string(get_band())
        eff_mass = c_eff_mass
        t0 = c_t0
        t1 = c_t1
        t2 = c_t2
        t3 = c_t3
        t4 = c_t4
        t5 = c_t5
        t6 = c_t6
        t7 = c_t7
        t8 = c_t8
        t9 = c_t9
        t10 = c_t10

![SUPERCONDUCTOR]
        method = get_string(get_method())
        FS_only = c_FS_only
        bcs_cutoff_frequency = c_bcs_cutoff_frequency
        num_eigenvalues_to_save = c_num_eigenvalues_to_save
        frequency_pts = c_frequency_pts

![RESPONSE]
        dynamic = c_dynamic
        ! End of loading variables
    end subroutine load_f90_config

end module confighub
