module ffirefly
    use iso_c_binding
    implicit none
    interface
        function epsilon(n, k) bind(C, name="epsilon_c")
            import :: C_INT, C_DOUBLE
            implicit none
            integer(C_INT), value :: n
            real(C_DOUBLE), dimension(3) :: k
            real(C_DOUBLE) :: epsilon
        end function epsilon

        function Vs(k1_c, k2_c, spin1_c, spin2_c) bind(C, name="Vs_c")
            import :: C_FLOAT, C_DOUBLE, C_CHAR
            real(C_DOUBLE), dimension(3), intent(in) :: k1_c
            real(C_DOUBLE), dimension(3), intent(in) :: k2_c
            character(C_CHAR), intent(in) :: spin1_c
            character(C_CHAR), intent(in) :: spin2_c
            real(C_FLOAT) :: V
        end function Vs

        function V(k1_c, k2_c) bind(C, name="V_c")
            import :: C_FLOAT, C_DOUBLE
            real(C_DOUBLE), dimension(3), intent(in) :: k1_c
            real(C_DOUBLE), dimension(3), intent(in) :: k2_c
            real(C_FLOAT) :: V
        end function V
    end interface
    ! Global variables

![CONTROL]
    character(len=50) :: category
    character(len=50) :: calculation
    character(len=50) :: method
    character(len=50) :: outdir
    character(len=50) :: indir
    character(len=50) :: prefix
    character(len=50) :: verbosity
    logical(c_bool), bind(C, name="c_automatic_file_read") :: c_automatic_file_read
    logical :: automatic_file_read
    logical(c_bool), bind(C, name="c_write_result") :: c_write_result
    logical :: write_result
    character(len=50) :: filetype

![SYSTEM]
    character(len=50) :: interaction
    integer(c_int), bind(C, name="c_dimension") :: c_dimension
    integer :: dimension
    character(len=50) :: celltype
    integer(c_int), bind(C, name="c_nbnd") :: c_nbnd
    integer :: nbnd
    integer(c_int), bind(C, name="c_natoms") :: c_natoms
    integer :: natoms
    real(c_float), bind(C, name="c_fermi_energy") :: c_fermi_energy
    real :: fermi_energy
    real(c_float), bind(C, name="c_Temperature") :: c_Temperature
    real :: Temperature
    real(c_float), bind(C, name="c_onsite_U") :: c_onsite_U
    real :: onsite_U
    real(c_float), bind(C, name="c_cutoff_energy") :: c_cutoff_energy
    real :: cutoff_energy

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

![ATOMIC_POSITIONS]
    character(len=50) :: atom
    real(c_float), bind(C, name="c_position") :: c_position(3)
    real :: position(3)

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
    logical(c_bool), bind(C, name="c_FS_only") :: c_FS_only
    logical :: FS_only
    integer(c_int), bind(C, name="c_num_eigenvalues_to_save") :: c_num_eigenvalues_to_save
    integer :: num_eigenvalues_to_save
    integer(c_int), bind(C, name="c_frequency_pts") :: c_frequency_pts
    integer :: frequency_pts
    character(len=50) :: projections

![RESPONSE]
    logical(c_bool), bind(C, name="c_dynamic") :: c_dynamic
    logical :: dynamic

![MANY_BODY]
    logical(c_bool), bind(C, name="c_self_consistent") :: c_self_consistent
    logical :: self_consistent
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
        function get_method() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_method
    end function get_method
        function get_outdir() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_outdir
    end function get_outdir
        function get_indir() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_indir
    end function get_indir
        function get_prefix() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_prefix
    end function get_prefix
        function get_verbosity() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_verbosity
    end function get_verbosity


        function get_filetype() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_filetype
    end function get_filetype

![SYSTEM]
        function get_interaction() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_interaction
    end function get_interaction

        function get_celltype() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_celltype
    end function get_celltype







![MESH]




![CELL]


![BRILLOUIN_ZONE]


![ATOMIC_POSITIONS]
        function get_atom() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_atom
    end function get_atom


![BANDS]
        function get_band() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_band
    end function get_band













![SUPERCONDUCTOR]



        function get_projections() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_projections
    end function get_projections

![RESPONSE]


![MANY_BODY]

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
        method = get_string(get_method())
        outdir = get_string(get_outdir())
        indir = get_string(get_indir())
        prefix = get_string(get_prefix())
        verbosity = get_string(get_verbosity())
        automatic_file_read = c_automatic_file_read
        write_result = c_write_result
        filetype = get_string(get_filetype())

![SYSTEM]
        interaction = get_string(get_interaction())
        dimension = c_dimension
        celltype = get_string(get_celltype())
        nbnd = c_nbnd
        natoms = c_natoms
        fermi_energy = c_fermi_energy
        Temperature = c_Temperature
        onsite_U = c_onsite_U
        cutoff_energy = c_cutoff_energy

![MESH]
        k_mesh = c_k_mesh
        q_mesh = c_q_mesh
        w_pts = c_w_pts

![CELL]
        cell = c_cell

![BRILLOUIN_ZONE]
        brillouin_zone = c_brillouin_zone

![ATOMIC_POSITIONS]
        atom = get_string(get_atom())
        position = c_position

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
        FS_only = c_FS_only
        num_eigenvalues_to_save = c_num_eigenvalues_to_save
        frequency_pts = c_frequency_pts
        projections = get_string(get_projections())

![RESPONSE]
        dynamic = c_dynamic

![MANY_BODY]
        self_consistent = c_self_consistent
        ! End of loading variables
    end subroutine load_f90_config

end module ffirefly
