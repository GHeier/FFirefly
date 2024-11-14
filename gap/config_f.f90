module confighub
    use iso_c_binding
    implicit none
    integer(c_int), bind(C, name="c_ibrav") :: c_ibrav 
    integer(c_int), bind(C, name="c_k_mesh") :: c_k_mesh(3)
    integer(c_int), bind(C, name="c_q_mesh") :: c_q_mesh(3)
    integer(c_int), bind(C, name="c_w_pts") :: c_w_pts
    integer(c_int), bind(C, name="c_dim") :: c_dim

    real(c_float), bind(C, name="c_mu") :: c_fermi_energy
    real(c_float), bind(C, name="c_U") :: c_onsite_U
    real(c_float), bind(C, name="c_wc") :: c_bcs_cutoff_frequency
    real(c_float), bind(C, name="c_cell") :: c_cell(3,3)
    real(c_float), bind(C, name="c_brillouin_zone") :: c_brillouin_zone(3,3)
    real(c_float), bind(C, name="c_max_freq") :: c_max_freq

    character(kind=c_char), bind(C, name="c_calculation_type") :: c_calculation_type(50)
    character(kind=c_char), bind(C, name="c_outdir") :: c_outdir(50)
    character(kind=c_char), bind(C, name="c_prefix") :: c_prefix(50)
    character(kind=c_char), bind(C, name="c_verbosity") :: c_verbosity(50)
    character(kind=c_char), bind(C, name="c_interaction") :: c_interaction(50)

    integer :: ibrav, k_mesh(3), q_mesh(3), w_pts, dim
    real :: fermi_energy, onsite_U, bcs_cutoff_frequency, cell(3,3), brillouin_zone(3,3), max_freq
    character(len=50) :: calculation_type, outdir, prefix, verbosity, interaction

    interface
        function get_calculation_type() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_calculation_type
        end function get_calculation_type

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

        function get_interaction() bind(C)
            use iso_c_binding
            type(c_ptr) :: get_interaction
        end function get_interaction

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
        ibrav = c_ibrav
        w_pts = c_w_pts
        dim = c_dim
        k_mesh = c_k_mesh
        q_mesh = c_q_mesh

        fermi_energy = c_fermi_energy
        onsite_U = c_onsite_U
        bcs_cutoff_frequency = c_bcs_cutoff_frequency
        cell = c_cell
        brillouin_zone = c_brillouin_zone
        max_freq = c_max_freq

        calculation_type = get_string(get_calculation_type())
        outdir = get_string(get_outdir())
        prefix = get_string(get_prefix())
        verbosity = get_string(get_verbosity())
        interaction = get_string(get_interaction())
    end subroutine load_f90_config

end module confighub
