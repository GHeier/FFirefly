module response
    use iso_c_binding
    use globals
    use mesh
    implicit none
    contains
     subroutine polarization_wrapper() bind(C)
        real(8), ALLOCATABLE :: chi_mesh_static(:,:,:)
        complex(8), ALLOCATABLE :: chi_mesh_dynamic(:,:,:,:)
        real(8) :: k(3), e
        call load_f90_config()
        call get_vals()
        if (dynamic) then
            allocate(chi_mesh_dynamic(ne, nge(1), nge(2), nge(3)))
            chi_mesh_dynamic = get_dynamic_polarization_mesh()
            call save_dynamic_mesh(chi_mesh_dynamic)
            deallocate(chi_mesh_dynamic)
        else
            allocate(chi_mesh_static(nge(1), nge(2), nge(3)))
            chi_mesh_static = get_static_polarization_mesh()
            call save_static_mesh(chi_mesh_static)
            deallocate(chi_mesh_static)
        end if
     end subroutine polarization_wrapper

     subroutine dos_wrapper() bind(C)
        real(8), ALLOCATABLE :: dos_mesh(:)
        real(8) :: emax, emin
        call load_f90_config()
        call get_vals()
        allocate(dos_mesh(w_pts))
        call get_full_DOS_spectrum(emax, emin, dos_mesh)
        call save_DOS_spectrum(dos_mesh, emax, emin)
        deallocate(dos_mesh)
     end subroutine DOS_wrapper
 end module response
