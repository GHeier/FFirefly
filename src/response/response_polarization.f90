module globals
    use, intrinsic :: iso_c_binding  ! Ensure C interoperability
    use ffirefly
    implicit none
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841971
    integer :: nb, ne, nge(3), ngw(3), nke, nkw
    real(8) :: ef, qmesh(3), bvec(3,3), VBZ, T

    contains

        subroutine get_vals()
            bvec = brillouin_zone
            nge = k_mesh
            ngw = q_mesh
            ef = fermi_energy
            ne = w_pts
            nb = 1
            VBZ = bvec(1,1) * bvec(2,2) * bvec(3,3) + bvec(1,2) * bvec(2,3) * bvec(3,1) &
            &   + bvec(1,3) * bvec(2,1) * bvec(3,2) - bvec(1,3) * bvec(2,2) * bvec(3,1) &
            &   + bvec(1,2) * bvec(2,1) * bvec(3,3) - bvec(1,1) * bvec(2,3) * bvec(3,2)
            if (dimension == 2) then
                VBZ = bvec(1,1) * bvec(2,2) - bvec(1,2) * bvec(2,1)
                nge(3) = 2
                ngw(3) = 1
            end if
            nke = product(nge)
            nkw = product(ngw)
            T = Temperature
        end subroutine get_vals
end module globals
module mesh
    use globals
    use libtetrabz
    implicit none
    contains

    function get_kvec(i1, i2, i3, npts) result(kvec)
        integer, intent(in) :: i1, i2, i3, npts(3)
        integer :: ivec(3)
        real(8), dimension(3) :: kvec
        ivec(1:3) = (/i1, i2, i3/)
        ivec(1:3) = modulo(ivec(1:3) + npts(1:3) / 2, npts(1:3)) - npts(1:3) / 2
        kvec(1:3) = dble(ivec(1:3)) / dble(npts(1:3)) - 0.5
    end function get_kvec

    function get_qvec(i1, i2, i3, npts) result(qvec)
        integer, intent(in) :: i1, i2, i3, npts(3)
        real(8), dimension(3) :: qvec, temp_pts
        temp_pts = npts
        if (dimension == 2) then
            temp_pts(3) = 2
        end if
        qvec(1:3) = dble([i1, i2, i3] + 1e-4) / dble(temp_pts(1:3) - 1) - 0.5
    end function get_qvec

    function fill_energy_mesh(qvec) result(eig)
        real(8), dimension(3), intent(in) :: qvec
        integer :: i1, i2, i3, ik
        real(8), dimension(3) :: kvec
        real(8), dimension(nb,nke) :: eig
        ik = 0
        do i3 = 0, nge(3) - 1
            do i2 = 0, nge(2) - 1
                do i1 = 0, nge(1) - 1
                    ik = ik + 1
                    kvec = get_kvec(i1, i2, i3, nge)
                    kvec(1:3) = matmul(bvec(1:3,1:3), kvec(1:3))
                    kvec(1:3) = kvec(1:3) + qvec(1:3)
                    eig(1,ik) = epsilon(1, kvec) - ef
                end do
            end do
        end do
     end function fill_energy_mesh

    function epsilon_mesh() result(eig)
        integer :: i1, i2, i3, ik
        real(8), dimension(3) :: kvec
        real(8), dimension(nb,nke) :: eig
        ik = 0
        do i3 = 0, nge(3) - 1
            do i2 = 0, nge(2) - 1
                do i1 = 0, nge(1) - 1
                    ik = ik + 1
                    kvec = get_kvec(i1, i2, i3, nge)
                    kvec(1:3) = matmul(bvec(1:3,1:3), kvec(1:3))
                    eig(1,ik) = epsilon(1, kvec) 
                end do
            end do
        end do
     end function epsilon_mesh

     subroutine get_full_DOS_spectrum(emax, emin, dos_list)
        real(8) :: dos_list(nkw)
        integer :: i, ltetra
        real(8) :: emax, emin, e0(ne)
        real(8), allocatable :: eig1(:,:), wght_dos(:,:,:)
        allocate(eig1(nb,nke), wght_dos(ne, nb, nkw))

        eig1 = epsilon_mesh()
        emax = maxval(eig1)
        emin = minval(eig1)
        ltetra = 2
        do i = 1, w_pts
            e0(i) = emin + dble(i-1) * (emax - emin) / dble(w_pts - 1)
        end do
        call libtetrabz_dos(ltetra, bvec, nb, nge, eig1, ngw, wght_dos, ne, e0)
        do i = 1, w_pts
            dos_list(i) = sum(wght_dos(i,1:nb,1:nkw)) / 2
        end do
    end subroutine get_full_DOS_spectrum

    function get_static_polarization(eig1, eig2, wght, wght_dos, qvec) result(chi)
        integer :: ltetra = 1
        real(8) :: eig1(:,:), qvec(3)
        real(8) :: eig2(:,:), wght(:,:,:), wght_dos(:,:,:)
        real(8) :: chi, qmag, qtmag, e0(nb)

        ltetra = 1
        qmag = sqrt(sum(qvec(1:3)**2))
        qtmag = sqrt(sum(qvec(1:2)**2))
        if (qmag < 1d-6 .OR. dimension == 2 .AND. qtmag < 1d-6) then
            e0(1) = 0d0
            call libtetrabz_dos(ltetra, bvec, nb, nge, eig1, ngw, wght_dos, ne, e0)
            chi = sum(wght_dos(1:ne,1:nb,1:nkw)) / 2
        else
            ! Calculate the polarization function
            call libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
            chi = sum(wght(1:nb,1:nb,1:nkw))
        end if
    end function get_static_polarization

    function get_static_polarization_mesh() result(chi_mesh)
        real(8) :: chi, chi_mesh(ngw(1), ngw(2), ngw(3))
        integer :: i1, i2, i3
        real(8), dimension(3) :: qvec
        real(8), allocatable :: eig1(:,:), eig2(:,:), wght(:,:,:), wght_dos(:,:,:)
        allocate(eig1(nb,nke), eig2(nb,nke), wght(nb,nb,nkw), wght_dos(ne,nb,nkw))

        print *, 'Calculating static polarization for ', ngw(1), 'x', ngw(2), 'x', ngw(3), ' k-points'
        eig1 = fill_energy_mesh([0d0, 0d0, 0d0])
        do i3 = 0, ngw(3) - 1
            print *, 'Calculating chi iteration ', i3+1, ' of ', ngw(3)
            !$OMP PARALLEL DO PRIVATE(i1, i2, qvec, eig2, chi, wght, wght_dos)
            do i2 = 0, ngw(2) - 1
                do i1 = 0, ngw(1) - 1
                    qvec = get_qvec(i1, i2, i3, ngw)
                    qvec(1:3) = matmul(bvec(1:3,1:3), qvec(1:3))
                    eig2 = fill_energy_mesh(qvec)
                    chi = get_static_polarization(eig1, eig2, wght, wght_dos, qvec)
                    chi_mesh(i1+1, i2+1, i3+1) = 2d0 * chi
                end do
            end do
        end do 
        deallocate(eig1, eig2, wght, wght_dos)
    end function get_static_polarization_mesh

    function get_dynamic_polarization(eig1, eig2, wght) result(chi)
        integer :: ltetra = 1
        real(8) :: eig1(:,:), qvec(3)
        real(8) :: eig2(:,:)
        complex(8) :: wght(:,:,:,:)
        complex(8) :: chi(ne)
        complex(8) :: e0(ne)
        integer :: i
        real(8) :: pi = 3.1415926535897932384626433832795028841971

        
        do i = 1, ne
            e0(i) = cmplx(0, pi * T * 2 * (i-1))
        end do
        ltetra = 1
        ! Calculate the polarization function
        call libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,0)
        chi = 2d0 * SUM(SUM(SUM(wght, DIM=4), DIM=3), DIM=2) * VBZ
    end function get_dynamic_polarization

    function get_dynamic_polarization_mesh() result(chi_mesh)
        complex(8) :: chi(ne), chi_mesh(ne, ngw(1), ngw(2), ngw(3)), max_w
        integer :: i1, i2, i3, i4
        real(8), dimension(3) :: qvec
        real(8), allocatable :: eig1(:,:), eig2(:,:)
        complex(8), allocatable :: wght(:,:,:,:)
        allocate(eig1(nb,nke), eig2(nb,nke), wght(ne,nb,nb,nkw))

        max_w = cmplx(0, pi * T * (2 * (ne-1) + 1))
        print *, 'Max w = ', max_w
        print *, 'Calculating dynamic polarization for ', ngw(1), 'x', ngw(2), 'x', ngw(3), ' k-points and ', ne, ' energy points'
        eig1 = fill_energy_mesh([0d0, 0d0, 0d0])
        do i1 = 0, ngw(1) - 1
            print *, 'Calculating chi iteration ', i1+1, ' of ', ngw(1)
            !$OMP PARALLEL DO PRIVATE(i1, i2, qvec, eig2, chi, wght)
            do i2 = 0, ngw(2) - 1
                do i3 = 0, ngw(3) - 1
                    qvec = get_qvec(i1, i2, i3, ngw)
                    qvec(1:3) = matmul(bvec(1:3,1:3), qvec(1:3))
                    eig2 = fill_energy_mesh(qvec)
                    chi = get_dynamic_polarization(eig1, eig2, wght)
                    chi_mesh(:, i1+1, i2+1, i3+1) = chi
                end do
            end do
        end do 
        deallocate(eig1, eig2, wght)
    end function get_dynamic_polarization_mesh

     subroutine save_static_mesh(chi_mesh)
        real(8) :: chi_mesh(ngw(1), ngw(2), ngw(3))
        integer :: i1, i2, i3
        real(8), dimension(3) :: qvec
        character(len=100) :: header, filename
        filename = trim(adjustl(trim(outdir) // trim(prefix))) // '_chi.dat'
        if (dimension == 3) then
            header = "         x             y             z             f"
        else if (dimension == 2) then
            header = "         x             y             f"
        else if (dimension == 1) then
            header = "         x             f"
        end if
        open(10, file=filename, status='unknown')
        write(10, '(A)') header
        do i1 = 0, ngw(1) - 1
            do i2 = 0, ngw(2) - 1
                do i3 = 0, ngw(3) - 1
                    qvec = get_qvec(i1, i2, i3, ngw)
                    qvec(1:3) = matmul(bvec(1:3,1:3), qvec(1:3))
                    write(10, '(1x, f13.6, 1x, f13.6, 1x, f13.6)') qvec(1), qvec(2), chi_mesh(i1+1, i2+1, i3+1)
                end do
            end do
        end do
        close(10)
        print *, 'Chi mesh saved to ', filename
     end subroutine save_static_mesh

     subroutine save_dynamic_mesh(chi_mesh)
        complex(8) :: chi_mesh(ne, ngw(1), ngw(2), ngw(3))
        integer :: i1, i2, i3, i4
        real(8), dimension(3) :: qvec
        real(8) :: pi = 3.1415926535897932384626433832795028841971
        complex(8) :: e0
        character(len=100) :: header, filename
        filename = trim(adjustl(trim(outdir) // trim(prefix))) // '_chi.dat'
        if (dimension == 3) then
            header = "         x             y             z             w            Re(f)         Im(f)"
        !else if (dimension == 2) then
        !    header = "         x             y             w            Re(f)         Im(f)"
        !else if (dimension == 1) then
        !    header = "         x             w            Re(f)         Im(f)"
        end if

        open(10, file=filename, status='unknown')
        write(10, '(A)') header
        do i4 = 0, ne - 1
            do i3 = 0, ngw(3) - 1
                do i2 = 0, ngw(2) - 1
                    do i1 = 0, ngw(1) - 1
                        qvec = get_qvec(i1, i2, i3, ngw)
                        qvec(1:3) = matmul(bvec(1:3,1:3), qvec(1:3))
                        e0 = cmplx(0, pi * T * (2 * i4))
                        write(10, '(1x, f13.6, 1x, f13.6, 1x, f13.6, 1x, f13.6, 1x, f13.6, 1x, f13.6)') qvec(1), qvec(2), qvec(3), aimag(e0), real(chi_mesh(i4+1, i1+1, i2+1, i3+1)), aimag(chi_mesh(i4+1, i1+1, i2+1, i3+1))
                    end do
                end do
            end do
        end do
        close(10)
        print *, 'Chi mesh saved to ', filename
     end subroutine save_dynamic_mesh

     subroutine save_DOS_spectrum(dos_list, emax, emin)
        real(8) :: dos_list(ne), emax, emin
        integer :: i
        character(len=100) :: header, filename
        filename = trim(adjustl(trim(outdir) // trim(prefix))) // '_DOS.dat'
        open(10, file=filename, status='unknown')
        header = "         w             f"
        write(10, '(A)') header
        do i = 0, w_pts - 1
            write(10, '(1x, f13.6, 1x, f13.6)') emin + dble(i) * (emax - emin) / dble(w_pts - 1), dos_list(i+1)
        end do
        close(10)
        print *, 'DOS spectrum saved to ', filename
     end subroutine save_DOS_spectrum

 end module mesh

