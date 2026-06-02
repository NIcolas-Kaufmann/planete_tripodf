module tripod

    use unitsPL
    include 'parameters.h'

    integer, parameter :: Nm_l = 5
    integer, parameter :: Nm_s = 2
    double precision, parameter :: f_fudge = 0.4
    double precision,dimension(nrad_max), parameter :: v_frag = 1.0d2 ! cm/s, fragmentation velocity
    double precision, parameter :: q_turb1 = -3.5 ! power law index for the turbulent relative velocity distribution
    double precision, parameter :: q_turb2 = -3.75 ! power law index
    double precision, parameter :: q_drfr = -3.75 ! power law index for the radial drift relative velocity distribution
    double precision, parameter :: q_sweep_tri = -3. ! power law index for the sweep-up relative velocity distribution
    double precision, parameter :: f_crit = 0.425d0 ! critical mass ratio for fragmentation, this can be adjusted as needed
    double precision, parameter,dimension(nrad_max,Nm_s) :: Sig_floor_tri = 1e-10 ! g/cm^2, floor for the surface density of the dust in each bin to avoid numerical issues, this can be adjusted as needed



    double precision, dimension(nrad_max,2) :: Sig_tri
    double precision, dimension(nrad_max) :: S
    double precision, dimension(nrad_max) :: a_max_tri
    double precision, dimension(nrad_max) :: a_min_tri
    double precision, dimension(nrad_max) :: q_rec
    double precision, dimension(nrad_max,5) :: a_tri
    double precision, dimension(nrad_max,5) :: m_tri 
    double precision, dimension(nrad_max,2) :: rho_tri
    double precision, dimension(nrad_max) :: rhos_tri
    double precision, dimension(nrad_max,5) :: H_tri
    double precision, dimension(nrad_max) :: fill_tri
    double precision, dimension(nrad_max,5) :: St_tri
    double precision, dimension(nrad_max,5) :: D_tri
    double precision, dimension(nrad_max,5) :: v_rad_tri
    double precision, dimension(nrad_max,5,5) :: v_rel_tot_tri
    double precision, dimension(nrad_max,5,5) :: v_rel_azi_tri
    double precision, dimension(nrad_max,5,5) :: v_rel_brown_tri
    double precision, dimension(nrad_max,5,5) :: v_rel_rad_tri
    double precision, dimension(nrad_max,5,5) :: v_rel_turb_tri
    double precision, dimension(nrad_max,5,5) :: v_rel_vert_tri
    double precision, dimension(nrad_max) :: alpha_vert_tri
    double precision, dimension(nrad_max) :: alpha_rad_tri
    double precision, dimension(nrad_max) :: alpha_turb_tri


    ! q and p stuff for determining the slope
    double precision, dimension(nrad_max) :: p_frag_tri
    double precision, dimension(nrad_max) :: p_fragtrans
    double precision, dimension(nrad_max) :: p_drfr_tri
    double precision, dimension(nrad_max) :: q_frag_tri
    double precision, dimension(nrad_max) :: q_eff_tri


    !timestep stuff 
    double precision :: ts_tri

    !Integrated quantities 
    double precision,dimension(nrad_max*3) :: rhs
    double precision,dimension(nrad_max) :: deriv_s_max
    double precision,dimension(nrad_max,Nm_s) :: S_rhs
    double precision,dimension(nrad_max,Nm_s) :: S_coag


    !boundary conditions
    double precision,parameter, dimension(Nm_s) :: inner_bc = [1e-5, 1e-5] ! small non-zero values to avoid numerical issues, these can be adjusted as needed
    double precision, parameter, dimension(Nm_s) :: outer_bc = [1e-5, 1e-5]
    double precision,parameter :: inner_s_bc = 1e-4
    double precision,parameter :: outer_s_bc = 1e-4


contains

!!!! 
! Initialize the arrays for the tripod module. This subroutine sets all the values to zero at the beginning of the simulation.
!!!

subroutine init_tripod()

end subroutine init_tripod

!!!


subroutine update_tripod(R,eta,T,mump,OmegaK,mfp,Sigma,cs,H_gas,dt,area,Ri)

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
    double precision, intent(in) :: dt
    double precision, intent(in) :: area(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    call update_dust(R,eta,T,mump,OmegaK,mfp,Sigma,cs,H_gas)
    call integrate_dust(area,R,Ri,Sigma,dt)
    
    !IO stuff

end subroutine update_tripod

!!!!
! this subroutine mimicks the dust.update from tripodpy 
!!!

subroutine update_dust(R,eta,T,mump,OmegaK,mfp,Sigma,cs,H_gas)
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
    ! add all the updater funtions here, for example
    ! tripodpy default updater
    !['delta', 'rhos', 'fill', 'backreaction', 'f', 'qrec', 'a', 'm', 'St', 'H', 'rho', 'D', 'eps', 'v', 'p', 'q', 'SigmaFloor', 'S'].
    call calc_q_rec(Sig_tri,a_min_tri,a_max_tri,q_rec,nrad_max)
    call calculate_a(a_min_tri,a_max_tri,q_rec,f_fudge,a_tri,nrad_max,Nm_l)
    call calculate_m(a_tri,rhos_tri,fill_tri,m_tri,nrad_max,Nm_l)
    call st_epstein_stokes1(a_tri,mfp,rhos_tri,Sigma,St_tri,nrad_max,Nm_l)
    call h_dubrulle1995(H_gas,St_tri,alpha_vert_tri,H_tri,nrad_max,Nm_l)
    !calculate the midplane density of the dust, which is needed for the collision rates
    rho_tri = Sig_tri/(sqrt(2.0d0*pi)*H_tri(:,[1,3]))
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
    call pfrag_trans( St_tri(:,Nm_l), alpha_turb_tri, Sigma, mump, p_fragtrans, nrad_max)
    call pdriftfrag(v_rel_rad_tri(:,4,5),v_rel_azi_tri(:,4,5),St_tri(:,Nm_l),alpha_rad_tri,Sigma,mump,cs,&
                    p_fragtrans,p_drfr_tri,nrad_max)
    call qfrag(p_drfr_tri,v_rel_tot_tri(:,4,5),v_frag,St_tri(:,Nm_l),q_turb1,q_turb2,q_drfr,alpha_turb_tri,Sigma,mump,q_frag_tri,nrad_max)
    q_eff_tri = q_frag_tri*p_frag_tri + q_sweep_tri*(1.0d0 - p_frag_tri)

    ! perarator sterp in tripodpy -> set the state vector rhs 
    rhs(1:nrad_max*Nm_s) = reshape(transpose(Sig_tri), [nrad_max*Nm_s])
    rhs((nrad_max*Nm_s)+1:(nrad_max*Nm_s)+nrad_max) = a_max_tri*Sig_tri(:,2)
    call smax_deriv(v_rel_tot_tri(:,4,5),rho_tri(:,2),rhos_tri, a_min_tri, a_max_tri,v_frag,Sig_tri,Sig_floor_tri,deriv_s_max,nrad_max,Nm_s)
    S_rhs = 0.0d0

    call s_coag(pi*(a_tri(:,[1,3])+a_tri(:,[3,2]))**2d0,v_rel_tot_tri(:,[1,3],[3,2]),H_tri(:,[1,3]),m_tri(:,[1,3]),Sig_tri,a_min_tri,a_max_tri,q_rec,Sig_floor_tri,S_coag,nrad_max,Nm_s)
    

end subroutine update_dust

subroutine Jacobian(Sigma,R,Ri,area,dt,dat_tot,row_tot,col_tot)
    implicit none


    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    double precision, intent(in) :: area(nrad_max)
    double precision, intent(in) :: dt
    double precision, intent(out), allocatable :: dat_tot(:)
    integer, intent(out), allocatable :: row_tot(:), col_tot(:)
    ! Local variables for the Jacobian construction
    double precision, dimension(nrad_max,2) :: A,B,C
    double precision, dimension((nrad_max-2)*Nm_s*Nm_s) :: dat_coag 
    integer, dimension((nrad_max-2)*Nm_s*Nm_s) :: row_coag, col_coag
    double precision, dimension(nrad_max,Nm_s) :: cross_section_tri
    integer, parameter :: n_dat_tot = (nrad_max-2)*Nm_s*Nm_s 
    integer, parameter :: N_tot = int(nrad_max*Nm_s)
    double precision :: dat_in(Nm_s*3),dat_out(Nm_s*3)
    integer :: row_in(Nm_s*3), col_in(Nm_s*3), row_out(Nm_s*3), col_out(Nm_s*3)

    double precision, allocatable :: dat_hydro(:)
    integer, allocatable :: row_hydro(:), col_hydro(:)
    integer :: i    




    cross_section_tri = pi*(a_tri(:,[1,3])+a_tri(:,[3,2]))**2d0
    call jacobian_coagulation_generator(cross_section_tri,v_rel_tot_tri(:,[1,3],[3,2]),H_tri(:,[1,3]),m_tri(:,[1,3]),Sig_tri,a_min_tri,a_max_tri,q_rec,&
                                        dat_coag,row_coag,col_coag,nrad_max,Nm_s)
    !Fortan idexing the arrays wer constucted with python sttyle indexing in mind 
    row_coag = row_coag + 1
    col_coag = col_coag + 1

    !construct the jacobian for the coagulation part and unravel the arrays C like in the pyhton version
    call jacobian_hydrodynamic_generator(area,D_tri(:,[1,3]),R,Ri,Sigma,v_rad_tri(:,3),A,B,C,nrad_max,Nm_s)

    !transpose the arrays first to match the C style ordering of indices and then reshape them to 1D arrays
    dat_hydro = [RESHAPE(Transpose(A(Nm_s+1:,:)), [SIZE(A(Nm_s+1:,:))]), RESHAPE(Transpose(B), [SIZE(B)]), RESHAPE(Transpose(C(1:nrad_max-Nm_s,:)), [SIZE(C(1:nrad_max-Nm_s,:))])]
    row_hydro = [(i+Nm_s, i=1,N_tot-Nm_s), (i, i=1,N_tot), (i, i=1,N_tot-Nm_s)]
    col_hydro = [(i, i=1,N_tot-Nm_s), (i, i=1,N_tot), (i+Nm_s, i=1,N_tot-Nm_s)]

    !inner boundary
    row_in = [(i, i=1,Nm_s), (i, i=1,Nm_s), (i, i=1,Nm_s)]
    col_in = [(i, i=1,Nm_s), (i+Nm_s, i=1,Nm_s), (i+2*Nm_s, i=1,Nm_s)]
    dat_in = 0.0d0
    !outer boundary
    row_out = [(N_tot-Nm_s+i, i=1,Nm_s), (N_tot-Nm_s+i, i=1,Nm_s), (N_tot-Nm_s+i, i=1,Nm_s)]
    col_out = [(N_tot-3*Nm_s+i, i=1,Nm_s), (N_tot-2*Nm_s+i, i=1,Nm_s), (N_tot-1*Nm_s+i, i=1,Nm_s)]
    dat_out = 0.0d0
    !todo Implement boundaries
    !val
    rhs(1:nm_s) = inner_bc
    rhs(N_tot-Nm_s+1:N_tot) = outer_bc

    dat_tot = [dat_hydro, dat_coag,dat_in, dat_out]
    row_tot = [row_hydro, row_coag, row_in, row_out]
    col_tot = [col_hydro, col_coag, col_in, col_out]



    !se

end subroutine Jacobian

subroutine Y_jacobian(area,R,Ri,Sigma,dt,values_J,rowind_J,colptr_J)    
    implicit none

    double precision, intent(in) :: area(nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: dt
    double precision, intent(out), allocatable :: values_J(:)
    integer, intent(out), allocatable :: rowind_J(:), colptr_J(:)
    ! Local variables for the Jacobian construction
    integer, parameter :: N_tot = int(nrad_max*Nm_s)
    double precision, dimension(nrad_max) :: A,B,C
    double precision, dimension(3) :: dat_in,dat_out
    integer, dimension(3) :: row_in,col_in,row_out,col_out
    integer :: col_diag(nrad_max*nm_s+nrad_max),row_diag(nrad_max*nm_s+nrad_max)
    double precision :: dat_diag(nrad_max*nm_s+nrad_max)
    integer :: i,k,nnz_diag
    logical :: found

    double precision, allocatable :: dat_J(:),dat_total(:),dat_hydro(:)
    integer, allocatable :: row_J(:), col_J(:),row_total(:), col_total(:),row_hydro(:), col_hydro(:)

    !get base JAcobian
    call Jacobian(Sigma,R,Ri,area,dt,dat_J,row_J,col_J)

    call jacobian_hydrodynamic_generator(area,D_tri(:,3),R,Ri,Sigma,v_rad_tri(:,3),A,B,C,nrad_max,Nm_s)

    dat_hydro = [A(2:nrad_max), B, C(1:nrad_max-1)]
    row_hydro = [(i+1, i=1,nrad_max-1), (i, i=1,nrad_max), (i, i=1,nrad_max-1)] + N_tot
    col_hydro = [(i, i=1,nrad_max-1), (i, i=1,nrad_max), (i+1, i=1,nrad_max-1)] + N_tot

    !boudary arrays
    row_in = 1 + N_tot
    col_in = [(i, i=1,3)] + N_tot
    dat_in = 0.0d0
    row_out = nrad_max-1 + N_tot
    col_out = [(nrad_max-3+i, i=1,3)] + N_tot
    dat_out = 0.0d0

    dat_total = [dat_J,dat_hydro, dat_in, dat_out]
    row_total = [row_J,row_hydro, row_in, row_out]
    col_total = [col_J,col_hydro, col_in, col_out]
    deallocate(dat_J, row_J, col_J, dat_hydro, row_hydro, col_hydro)


    ! make the actual integration matrix by subtracting eye *dt 
    nnz_diag = 0
    do i = 1, Nrad_max*Nm_s+nrad_max
        found = .false.
        do k = 1, size(dat_total)
            if (row_total(k) == i .and. col_total(k) == i) then
                dat_total(k) = dat_total(k) - dt
                found = .true.
                exit
            end if
        end do

        if (.not. found) then
            nnz_diag = nnz_diag + 1
            row_diag(nnz_diag) = i
            col_diag(nnz_diag) = i
            dat_diag(nnz_diag) = -dt
        end if
    end do

    if (nnz_diag > 0) then
        dat_total = [dat_total, dat_diag(1:nnz_diag)]
        row_total = [row_total, row_diag(1:nnz_diag)]
        col_total = [col_total, col_diag(1:nnz_diag)]
    end if

    allocate(values_J(SIZE(dat_total)), rowind_J(SIZE(dat_total)), colptr_J(N_tot+nrad_max+1))
    call triplet_to_csc(N_tot+nrad_max, SIZE(dat_total), row_total, col_total, dat_total, &
                        colptr_J, rowind_J, values_J)

    !set the boundary conditions in the rhs vector preliminarary -> only value is implmeneted at the moment 
    rhs(N_tot+1)=inner_s_bc
    rhs(N_tot+nrad_max)=outer_s_bc




end subroutine Y_jacobian

 
subroutine integrate_dust(area,R,Ri,Sigma,dt)
    implicit none

    double precision, intent(in) :: area(nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: dt

    double precision, allocatable :: values_J(:)
    integer, allocatable :: rowind_J(:), colptr_J(:)    
    call Y_jacobian(area,R,Ri,Sigma,dt,values_J,rowind_J,colptr_J)

    !implement the S coag source for the rhs term should work without though
    rhs(1:nrad_max*Nm_s) = rhs(1:nrad_max*Nm_s) + dt * reshape(transpose(S_rhs), [nrad_max*Nm_s])
    rhs((nrad_max*Nm_s)+1:(nrad_max*Nm_s)+nrad_max) = rhs((nrad_max*Nm_s)+1:(nrad_max*Nm_s)+nrad_max) + dt * ((deriv_s_max*Sig_tri(:,2))+(a_max_tri*(S_rhs(:,2)+S_coag(:,2))))

    call solve_superlu(SIZE(values_J), 1, values_J, rowind_J, colptr_J, rhs)
    deallocate(values_J, rowind_J, colptr_J)
    call finalize_integration()
    
end subroutine integrate_dust


subroutine solve_superlu(nzz_max,nrhs,values,rowind,colptr,b)
    implicit none

    integer, intent(in) :: nzz_max, nrhs
    real*8, intent(in) :: values(nzz_max)
    integer, intent(in) :: rowind(nzz_max), colptr(nrad_max*Nm_s+nrad_max+1)
    real*8, intent(inout) :: b(nrad_max*Nm_s+nrad_max)
    ! SuperLU dimensions and variables
    integer :: ldb, iopt,n,nnz
    ! Note: 'factors' must be an 8-byte integer to hold a C pointer (fptr)
    integer*8 :: factors 
    
    ! Depending on whether SuperLU was compiled with 64-bit integer indices 
    ! (XSDK_INDEX_SIZE=64), these might need to be integer*8.
    ! We use standard integers here as default.
    integer :: info

    ! 1. Factorize the matrix J
    ldb = nrad_max*Nm_s+nrad_max
    n = ldb
    nnz = nzz_max
    iopt = 1
    call c_fortran_dgssv(iopt, n, nnz, nrhs, values, rowind, colptr, &
                         b, ldb, factors, info)

    if (info /= 0) then
        write(*,*) 'SuperLU Factorization failed, info = ', info
        !deallocate(values, rowind, colptr, b)
        stop
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


end subroutine solve_superlu


subroutine finalize_integration()
    implicit none

    Sig_tri = reshape(rhs(1:nrad_max*Nm_s), [nrad_max,Nm_s], order=[2,1])
    a_max_tri = rhs(nrad_max*Nm_s+1:nrad_max*Nm_s+nrad_max)/Sig_tri(:,2)
    a_max_tri = max(a_max_tri, 1.5d0 * a_min_tri) ! enforce that a_max is not smaller than a_min to avoid numerical issues, this can be adjusted as needed
    call enforce_f() ! enforce that the fragmentation barrier is respected by adjusting the surface density in the largest bin, this is a simple fix to avoid numerical issues and can be adjusted as needed
end subroutine finalize_integration


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


! boiler plate subroutine to convert from triplet format (row, col, val) to compressed sparse column (CSC) format (colptr, rowind, values)
subroutine triplet_to_csc(ncols, nnz, row, col, val, &
                          colptr, rowind, values)

    implicit none

    integer, intent(in) :: ncols
    integer, intent(in) :: nnz

    integer, intent(in) :: row(nnz)
    integer, intent(in) :: col(nnz)
    real(8), intent(in) :: val(nnz)

    integer, intent(out) :: colptr(ncols+1)
    integer, intent(out) :: rowind(nnz)
    real(8), intent(out) :: values(nnz)

    integer :: counts(ncols)
    integer :: next(ncols)
    integer :: k, c, pos

    !---------------------------------------
    ! Count entries in each column
    !---------------------------------------
    counts = 0

    do k = 1, nnz
        counts(col(k)) = counts(col(k)) + 1
    end do

    !---------------------------------------
    ! Build column pointers
    !---------------------------------------
    colptr(1) = 1

    do c = 1, ncols
        colptr(c+1) = colptr(c) + counts(c)
    end do

    !---------------------------------------
    ! Working copy of column starts
    !---------------------------------------
    next = colptr(1:ncols)

    !---------------------------------------
    ! Fill CSC arrays
    !---------------------------------------
    do k = 1, nnz

        c = col(k)

        pos = next(c)

        rowind(pos) = row(k)
        values(pos) = val(k)

        next(c) = next(c) + 1

    end do

end subroutine triplet_to_csc


subroutine enforce_f()
double precision,dimension(nrad_max) :: delta 
integer :: i

delta = f_crit*sum(Sig_tri,dim=2) -Sig_tri(:,2)
do i = 1, nrad_max
    delta(i) = max(0.0d0, delta(i))
end do
Sig_tri(:,2) = Sig_tri(:,2) + delta
Sig_tri(:,1) = Sig_tri(:,1) - delta
end subroutine enforce_f


subroutine write_output()
    implicit none

    ! This subroutine should handle all the output writing, for example writing the dust surface density and maximum grain size to files for post-processing and visualization. The implementation can be adjusted as needed, for example by using different file formats or adding more output variables.

end subroutine write_output

end module tripod

