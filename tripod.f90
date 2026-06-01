module tripod
    use unitsPL
    implicit none
    integer, parameter :: Nm_l = 5
    integer, parameter :: Nm_s = 2
    double precision, parameter :: f_fudge = 0.4
    double precision, parameter :: v_frag = 1.0d2 ! cm/s, fragmentation velocity
    double precision, parameter :: q_turb1 = -3.5 ! power law index for the turbulent relative velocity distribution
    double precision, parameter :: q_turb2 = -3.75 ! power law index
    double precision, parameter :: q_drfr = -3.75 ! power law index for the radial drift relative velocity distribution
    double precision, parameter :: q_sweep_tri = -3. ! power law index for the sweep-up relative velocity distribution




    double precision dimension(nrad_max,2) :: Sig_tri
    double precision dimension(nrad_max) :: S
    double precision dimension(nrad_max) :: a_max_tri
    double precision dimension(nrad_max) :: a_min_tri
    double precision dimension(nrad_max) :: q_rec
    double precision dimension(nrad_max,5) :: a_tri
    double precision dimension(nrad_max,5) :: m_tri 
    double precision dimension(nrad_max) :: rho_tri
    double precision dimension(nrad_max) :: rhos_tri
    double precision dimension(nrad_max) :: H_tri
    double precision dimension(nrad_max) :: fill_tri
    double precision dimension(nrad_max) :: St_tri
    double precision dimension(nrad_max,5) :: D_tri
    double precision dimension(nrad_max,5) :: v_rad_tri
    double precision dimension(nrad_max,5,5) :: v_rel_tot_tri
    double precision dimension(nrad_max,5,5) :: v_rel_azi_tri
    double precision dimension(nrad_max,5,5) :: v_rel_brown_tri
    double precision dimension(nrad_max,5,5) :: v_rel_rad_tri
    double precision dimension(nrad_max,5,5) :: v_rel_turb_tri
    double precision dimension(nrad_max,5,5) :: v_rel_vert_tri
    double precision dimension(nrad_max) :: alpha_vert_tri
    double precision dimension(nrad_max) :: alpha_rad_tri
    double precision dimension(nrad_max) :: alpha_turb_tri

    double precision dimension(nrad_max,5) :: D_tri

    ! q and p stuff for determining the slope
    double precision dimension(nrad_max) :: p_frag_tri
    double precision dimension(nrad_max) :: p_fragtrans
    double precision dimension(nrad_max) :: p_drfr_tri
    double precision dimension(nrad_max) :: q_frag_tri
    double precision dimension(nrad_max) :: q_eff_tri


    !timestep stuff 
    double precision :: ts_tri


    !Jacobian stuff




contains

!!!! 
! Initialize the arrays for the tripod module. This subroutine sets all the values to zero at the beginning of the simulation.
!!!

subroutine init_tripod
    implicit none

    integer :: i

    do i = 1, nrad_max
        a_max_tri(i) = 0.0d0
        a_min_tri(i) = 0.0d0
        q_rec(i) = 0.0d0
        q_eq(i) = 0.0d0
        rhos_tri(i) = 0.0d0
        St_tri(i) = 0.0d0
        p_frag_tri(i) = 0.0d0
        p_fragtrans(i) = 0.0d0
        p_drfr_tri(i) = 0.0d0
        q_frag_tri(i) = 0.0d0
        q_eff_tri(i) = 0.0d0
    end do

    ts_tri = 0.0d0

end subroutine init_tripod

!!!


subroutine update_tripod()

    implicit none
    call update_dust(mfp,Sigma,cs,H_gas,nrad_max)
    call integrate_dust()
    
end subroutine update_tripod

!!!!
! this subroutine mimicks the dust.update from tripodpy 
!!!

subroutine update_dust(OmegaK,mfp,Sigma,cs,H_gas,nrad_max)
    implicit none

    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: eta(nrad_max)
    double precision, intent(in) :: T(nrad_max)
    double precision, intent(in) :: mump(nrad_max)
    double precision, intent(in) :: OmegaK(nrad_max)
    double precision, intent(in) :: mfp(nrad_max)
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: cs(nrad_max)
    double precision, intent(in) :: H_gas(nrad_max)
    double precision, intent(in) :: nrad_max
    ! add all the updater funtions here, for example
    ! tripodpy default updater
    !['delta', 'rhos', 'fill', 'backreaction', 'f', 'qrec', 'a', 'm', 'St', 'H', 'rho', 'D', 'eps', 'v', 'p', 'q', 'SigmaFloor', 'S'].
    call calc_q_rec(Sig_tri,a_min_tri,a_max_tri,q_rec,nrad_max)
    call calculate_a(a_min_tri,a_max_tri,q_rec,f_fudge,a_tri,nrad_max,Nm_l)
    call calculate_m(a_tri,rhos_tri,fill_tri,m_tri,nrad_max,Nm_l)
    call st_epstein_stokes1(a_tri,mfp,rhos_tri,Sigma,St_tri,nrad_max,Nm_l)
    call h_dubrulle1995(H_gas,St_tri,alpha_vert_tri,H_tri,nrad_max,Nm_l)
    !calculate the midplane density of the dust, which is needed for the collision rates
    rho_tri = Sigma/(sqrt(2.0d0*pi)*H_tri)
    call d(alpha_rad_tri*cs**2, OmegaK, St_tri, D_tri, nrad_max, Nm_l)

    call vrad(St_tri, eta*R*OmegaK, v_rad_tri, nrad_max, Nm_l)
    ! Relative velocities
    call vrel_azimuthal_drift(eta*R*OmegaK, St_tri, v_rel_azi_tri, nrad_max, Nm_l)
    call vrel_brownian_motion(cs, m_tri, T, v_rel_brown_tri, nrad_max, Nm_l)
    call vrel_radial_drift(v_rad_tri, v_rel_rad_tri, nrad_max, Nm_l)
    call vrel_ormel_cuzzi_2007(alpha_turb_tri, cs, mump, OmegaK, Sigma, St_tri, v_rel_turb_tri, nrad_max, Nm_l)
    call vrel_vertical_settling(H_tri, OmegaK, St_tri, v_rel_vert_tri, nrad_max, Nm_l)
    v_rel_tot_tri = sqrt(v_rel_azi_tri**2 + v_rel_brown_tri**2 + v_rel_rad_tri**2 + v_rel_turb_tri**2 + v_rel_vert_tri**2)

    ! collision outcomes p and q 
    call pfrag(v_rel_tot_tri, v_frag, p_frag_tri, nrad_max, Nm_l)
    call pfrag_trans( St_tri(:,Nm_l), alpha_turb_tri, Sigma, mump, p_fragtrans, Nr)
    call pdriftfrag(v_rel_rad_tri(:,4,5),v_rel_frag_azi_tri(:,4,5),St_tri(:,Nm_l),alpha_rad_tri,Sigma,mump,cs,p_fragtrans,p_drfr_tri,nrad_max)
    call qfrag(p_dr,v_ret_tot_tri(:,4,5),v_frag,St_tri(:,Nm_l),q_turb1,q_turb2,q_drfr,alpha_turb_tri,Sigma,mump,q_frag_tri,nrad_max)
    q_eff_tri = q_frag_tri*p_frag_tri + q_sweep_tri*(1.0d0 - p_frag_tri)


subroutine construct_Jacobian()

    implicit none
    ! construct the Jacobian matrix for the system of equations governing the dust evolution
    ! this will be used for the implicit time-stepping scheme

end subroutine construct_Jacobian

! This subroutine constructs a dummy sparse matrix J (tridiagonal) of size nrad_max x nrad_max
! and solves the linear system J * x = b using SuperLU.
subroutine solve_dummy_superlu()
    implicit none

    ! SuperLU dimensions and variables
    integer :: n, nrhs, ldb, iopt
    ! Note: 'factors' must be an 8-byte integer to hold a C pointer (fptr)
    integer*8 :: factors 
    
    ! Depending on whether SuperLU was compiled with 64-bit integer indices 
    ! (XSDK_INDEX_SIZE=64), these might need to be integer*8.
    ! We use standard integers here as default.
    integer :: nnz_max, nnz, info
    integer, allocatable :: rowind(:), colptr(:)
    real*8, allocatable  :: values(:), b(:)
    
    integer :: j

    n = nrad_max
    nnz_max = 3 * n ! Tridiagonal matrix: main diag, upper diag, lower diag

    allocate(values(nnz_max))
    allocate(rowind(nnz_max))
    allocate(colptr(n + 1))
    allocate(b(n))

    ! Construct a dummy tridiagonal matrix J = diag(-2) + super/subdiag(1)
    ! in Compressed Column Storage (CSC) format.
    nnz = 0
    do j = 1, n
        colptr(j) = nnz + 1
        
        ! sub-diagonal
        if (j > 1) then
            nnz = nnz + 1
            values(nnz) = 1.0d0
            rowind(nnz) = j - 1
        end if
        
        ! diagonal
        nnz = nnz + 1
        values(nnz) = -2.0d0
        rowind(nnz) = j
        
        ! super-diagonal
        if (j < n) then
            nnz = nnz + 1
            values(nnz) = 1.0d0
            rowind(nnz) = j + 1
        end if
    end do
    colptr(n + 1) = nnz + 1

    ! Set up a dummy right-hand side `b`
    b = 1.0d0
    nrhs = 1
    ldb = n

    ! 1. Factorize the matrix J
    iopt = 1
    call c_fortran_dgssv(iopt, n, nnz, nrhs, values, rowind, colptr, &
                         b, ldb, factors, info)

    if (info /= 0) then
        write(*,*) 'SuperLU Factorization failed, info = ', info
        deallocate(values, rowind, colptr, b)
        return
    else
        write(*,*) 'SuperLU Factorization succeeded.'
    end if

    ! 2. Solve the linear system
    iopt = 2
    call c_fortran_dgssv(iopt, n, nnz, nrhs, values, rowind, colptr, &
                         b, ldb, factors, info)

    if (info == 0) then
        write(*,*) 'SuperLU Solve succeeded. First 5 solution elements:'
        write(*,*) b(1:min(5, n))
    else
        write(*,*) 'SuperLU Solve failed, info = ', info
    end if

    ! 3. Free the internal storage allocated by SuperLU
    iopt = 3
    call c_fortran_dgssv(iopt, n, nnz, nrhs, values, rowind, colptr, &
                         b, ldb, factors, info)

    deallocate(values, rowind, colptr, b)

end subroutine solve_dummy_superlu

end module tripod
