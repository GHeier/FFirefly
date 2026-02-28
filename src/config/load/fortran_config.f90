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
    logical(c_bool), bind(C, name="c_debug") :: c_debug
    logical :: debug
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
    real(c_float), bind(C, name="c_fermi_energy") :: c_fermi_energy
    real :: fermi_energy
    real(c_float), bind(C, name="c_num_electrons") :: c_num_electrons
    real :: num_electrons
    logical(c_bool), bind(C, name="c_mu_from_n") :: c_mu_from_n
    logical :: mu_from_n
    real(c_float), bind(C, name="c_Temperature") :: c_Temperature
    real :: Temperature
    real(c_float), bind(C, name="c_cutoff_energy") :: c_cutoff_energy
    real :: cutoff_energy
    real(c_float), bind(C, name="c_smearing") :: c_smearing
    real :: smearing
    real(c_float), bind(C, name="c_mixing") :: c_mixing
    real :: mixing
    integer(c_int), bind(C, name="c_max_iters") :: c_max_iters
    integer :: max_iters
    real(c_float), bind(C, name="c_qp_weight") :: c_qp_weight
    real :: qp_weight

![HAMILTONIAN]
    character(len=50) :: hamiltonian

![HUBBARD]
    real(c_float), bind(C, name="c_U0") :: c_U0
    real :: U0
    real(c_float), bind(C, name="c_U1") :: c_U1
    real :: U1
    real(c_float), bind(C, name="c_J0") :: c_J0
    real :: J0
    real(c_float), bind(C, name="c_J1") :: c_J1
    real :: J1

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

![BASIS]
    character(len=50) :: states(50)
    real(c_float), bind(C, name="c_positions") :: c_positions(50,3)
    real :: positions(50,3)

![BANDS]
    character(len=50) :: band
    real(c_float), bind(C, name="c_eff_mass") :: c_eff_mass
    real :: eff_mass
    real(c_float), bind(C, name="c_t0") :: c_t0
    real :: t0
    real(c_float), bind(C, name="c_t1") :: c_t1
    real :: t1
    real(c_float), bind(C, name="c_t2") :: c_t2
    real :: t2
    real(c_float), bind(C, name="c_t3") :: c_t3
    real :: t3
    real(c_float), bind(C, name="c_t4") :: c_t4
    real :: t4
    real(c_float), bind(C, name="c_t5") :: c_t5
    real :: t5
    real(c_float), bind(C, name="c_t6") :: c_t6
    real :: t6
    real(c_float), bind(C, name="c_t7") :: c_t7
    real :: t7
    real(c_float), bind(C, name="c_t8") :: c_t8
    real :: t8
    real(c_float), bind(C, name="c_t9") :: c_t9
    real :: t9
    real(c_float), bind(C, name="c_t10") :: c_t10
    real :: t10

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











![HAMILTONIAN]
        function get_hamiltonian() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_hamiltonian
    end function get_hamiltonian

![HUBBARD]





![MESH]




![CELL]


![BRILLOUIN_ZONE]


![BASIS]
        function get_states() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_states
    end function get_states


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

    function get_string_from_array(c_string_array, index) result(fortran_string)
        type(c_ptr), intent(in) :: c_string_array
        integer, intent(in) :: index
        character(len=:), allocatable :: fortran_string
        type(c_ptr) :: c_string_ptr
        type(c_ptr), pointer :: array_of_ptrs(:)
        character(kind=c_char, len=1), pointer :: char_array(:)
        integer :: i, length

        ! Convert the c_ptr to an array of c_ptr
        call c_f_pointer(c_string_array, array_of_ptrs, [50])

        ! Get the pointer to the specific string
        c_string_ptr = array_of_ptrs(index + 1)

        ! Convert that string
        call c_f_pointer(c_string_ptr, char_array, [1000])
        length = 0
        do i = 1, size(char_array)
            if (char_array(i) == c_null_char) exit
            length = length + 1
        end do
        allocate(character(len=length) :: fortran_string)
        do i = 1, length
            fortran_string(i:i) = char_array(i)
        end do
    end function get_string_from_array

    subroutine load_f90_config()
        integer :: i, j
        ! Load variables

![CONTROL]
        category = get_string(get_category())
        calculation = get_string(get_calculation())
        method = get_string(get_method())
        outdir = get_string(get_outdir())
        debug = c_debug
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
        fermi_energy = c_fermi_energy
        num_electrons = c_num_electrons
        mu_from_n = c_mu_from_n
        Temperature = c_Temperature
        cutoff_energy = c_cutoff_energy
        smearing = c_smearing
        mixing = c_mixing
        max_iters = c_max_iters
        qp_weight = c_qp_weight

![HAMILTONIAN]
        hamiltonian = get_string(get_hamiltonian())

![HUBBARD]
        U0 = c_U0
        U1 = c_U1
        J0 = c_J0
        J1 = c_J1

![MESH]
        k_mesh = c_k_mesh
        q_mesh = c_q_mesh
        w_pts = c_w_pts

![CELL]
        cell = c_cell

![BRILLOUIN_ZONE]
        brillouin_zone = c_brillouin_zone

![BASIS]
        do i = 1, nbnd
            states(i) = get_string_from_array(get_states(), i-1)
        end do
        do i = 1, nbnd
            do j = 1, 3
                positions(i,j) = c_positions(i,j)
            end do
        end do

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
