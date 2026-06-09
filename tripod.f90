module tripod

    use unitsPL
    include 'parameters.h'

    integer, parameter :: Nm_l = 5
    integer, parameter :: Nm_s = 2
    double precision, parameter :: f_fudge = 0.4
    double precision,dimension(nrad_max), parameter :: v_frag = 1.0d3 ! cm/s, fragmentation velocity
    double precision, parameter :: q_turb1 = -3.5 ! power law index for the turbulent relative velocity distribution
    double precision, parameter :: q_turb2 = -3.75 ! power law index
    double precision, parameter :: q_drfr = -3.75 ! power law index for the radial drift relative velocity distribution
    double precision, parameter :: q_sweep_tri = -3. ! power law index for the sweep-up relative velocity distribution
    double precision, parameter :: f_crit = 0.425d0 ! critical mass ratio for fragmentation, this can be adjusted as needed
    double precision, parameter :: a_lim = 1e-4 ! minimal shrikage size in cm
    double precision, parameter :: f_drift = 0.8
    double precision, parameter,dimension(nrad_max,Nm_s) :: Sig_floor_tri = 1e-10 ! g/cm^2, floor for the surface density of the dust in each bin to avoid numerical issues, this can be adjusted as needed
    double precision, parameter :: cfl_tri = 1d-1



    double precision, dimension(nrad_max,2) :: Sig_tri
    double precision, dimension(nrad_max) :: S
    double precision, dimension(nrad_max) :: a_max_tri
    double precision, dimension(nrad_max) :: a_min_tri
    double precision, dimension(nrad_max) :: q_rec
    double precision, dimension(nrad_max,5) :: a_tri
    double precision, dimension(nrad_max,5) :: m_tri 
    double precision, dimension(nrad_max,2) :: rho_tri
    double precision, dimension(nrad_max,5) :: rhos_tri
    double precision, dimension(nrad_max,5) :: H_tri
    double precision, dimension(nrad_max,5) :: fill_tri
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
    double precision, dimension(nrad_max) :: smax_dot_hyd

    !Integrated quantities 
    double precision,dimension(nrad_max*3) :: rhs
    double precision,dimension(nrad_max) :: deriv_s_max
    double precision,dimension(nrad_max,Nm_s) :: S_rhs
    double precision,dimension(nrad_max,Nm_s) :: S_coag_tri
    double precision, dimension(nrad_max,Nm_s) :: S_hyd_tri
    double precision, dimension(nrad_max,Nm_s) :: S_tot_tri


    !boundary conditions
    character(len=*),parameter ::bd_inner_type = "const_grad"
    character(len=*),parameter :: bd_outer_type = "val"
    double precision, dimension(Nm_s) :: inner_bc = [1e-5, 1e-5] ! small non-zero values to avoid numerical issues, these can be adjusted as needed
    double precision, dimension(Nm_s) :: outer_bc = [1e-11, 1e-11]

    character(len=*),parameter :: s_bd_inner_type = "const_grad"
    double precision,parameter :: inner_s_bc = 1e-4
    double precision,parameter :: outer_s_bc = 7.8431328811173312E-005

    !output name
    character(len=*),parameter :: outfile_name = "gap"

    !grid_stuff 
    double precision, dimension(nrad_max+1) :: Ri_tri

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

subroutine initialize_dust(a_min_ini,a_max_ini,alpha_rad,alpha_vert,alpha_turb,fd2g,rhos,Sigma,P,cs,T,H_gas,mump,eta,mfp,OmegaK,R)

    implicit none 

    double precision, intent(in) :: a_min_ini 
    double precision, intent(in) :: a_max_ini
    double precision, intent(in) :: alpha_rad
    double precision, intent(in) :: alpha_vert
    double precision, intent(in) :: alpha_turb
    double precision, intent(in) :: fd2g 
    double precision, intent(in) :: rhos
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: P(nrad_max)
    double precision, intent(in) :: cs(nrad_max)
    double precision, intent(in) :: T(nrad_max)
    double precision, intent(in) :: H_gas(nrad_max)
    double precision, intent(in) :: mump(nrad_max)
    double precision, intent(in) :: eta(nrad_max)
    double precision, intent(in) :: mfp(nrad_max)
    double precision, intent(in) :: OmegaK(nrad_max)
    double precision, intent(in) :: R(nrad_max)

    !set smin not a parameter in principle but should not be changed 
    a_min_tri = a_min_ini
    alpha_rad_tri = alpha_rad
    alpha_turb_tri = alpha_turb
    alpha_vert_tri = alpha_vert
    rhos_tri = rhos
    !print *, "grid."
    !initalize the ri gird as planete does not have one
    call log_grid_interfaces(nrad_max, R)

    !print *, "amax_ini."
    call a_max_initial(a_max_ini,fd2g,Sigma,cs,P,OmegaK,R,Ri_tri)
    
    !print *, "Sig_ini."
    call Sigma_initial(a_min_tri, a_max_tri,Sigma, fd2g)
    print *, "update_ini."
    call update_dust(R,eta,T,mump,OmegaK,mfp,Sigma,cs,H_gas)
    print *, "setup boundaries"
    !set boundaries 
    inner_bc = Sig_tri(2,:)
    outer_bc = Sig_tri(nrad_max-1,:)
    print *, "update_done."
end subroutine initialize_dust


subroutine a_max_initial(a_max_ini,fd2g,Sigma,cs,P,OmegaK,R,Ri)
    use interpolation
    implicit none 
    
    double precision, intent(in) :: a_max_ini
    double precision, intent(in) :: fd2g
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: cs(nrad_max)
    double precision, intent(in) :: P(nrad_max)
    double precision, intent(in) :: OmegaK(nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    !local variables 
    double precision :: P_i(nrad_max+1)
    double precision :: gam(nrad_max)
    double precision :: ad(nrad_max)
    integer :: i



    call interp1d(ri, r, P, P_i, nrad_max)

    gam = abs((P_i(2:)-P_i(:nrad_max))/(Ri(2:)-Ri(:nrad_max))) * R/P

    ad = 5d-3 * 2d0 /pi * fd2g * Sigma * rhos_tri(:,1) * (OmegaK*R)**2 /cs**2 /gam
    

    do i = 1, nrad_max
        a_max_tri(i) = max(1.5d0*a_min_tri(i),min(ad(i),a_max_ini))
    enddo 

end subroutine


subroutine Sigma_initial(s_min, s_max,gas_Sigma, d2gRatio)
  ! Calculates the initial condition of the dust surface densities.
  !
  ! Parameters
  ! ----------
  ! nrad_max         : number of radial grid cells
  ! Nm_s         : number of size bins (expected to be 2)
  ! q_eff      : effective size distribution slope      (nrad_max)
  ! s_min      : minimum grain size per cell            (nrad_max)
  ! s_max      : maximum grain size per cell            (nrad_max)
  ! SigmaFloor : floor surface density                  (nrad_max, Nm_s)
  ! gas_Sigma  : gas surface density                    (nrad_max)
  ! d2gRatio   : dust-to-gas ratio                      (scalar)
  !
  ! Output
  ! ------
  ! Sig_tri      : initial dust surface density           (nrad_max, Nm_s)

  implicit none

  real(8), intent(in)  :: s_min(nrad_max)
  real(8), intent(in)  :: s_max(nrad_max)
  real(8), intent(in)  :: gas_Sigma(nrad_max)
  real(8), intent(in)  :: d2gRatio

  ! Local variables
  real(8) :: sint, qp4, S0, S1
  integer :: i

  real(8), parameter :: Q4 = -4.0d0

  q_rec = -3.5d0
  do i = 1, nrad_max

    sint = sqrt(s_min(i) * s_max(i))   ! geometric mean of s_min and s_max
    qp4  = q_rec(i) + 4.0d0

    ! ------------------------------------------------------------------
    ! Compute fractional weights S0, S1
    ! ------------------------------------------------------------------
    if (q_rec(i) == Q4) then

      ! q == -4 branch: use logarithmic weights
      if (s_max(i) == s_min(i)) then
        S0 = Sig_floor_tri(i, 1)
        S1 = Sig_floor_tri(i, 2)
      else
        S0 = log(sint      / s_min(i)) / log(s_max(i) / s_min(i))
        S1 = 1.0d0 - S0
      end if

    else

      ! q != -4 branch: use power-law weights
      if (s_max(i) <= 1.5d0 * s_min(i)) then
        S0 = Sig_floor_tri(i, 1)
        S1 = Sig_floor_tri(i, 2)
      else
        S0 = (sint**qp4 - s_min(i)**qp4) / (s_max(i)**qp4 - s_min(i)**qp4)
        S1 = 1.0d0 - S0
      end if

    end if

    ! ------------------------------------------------------------------
    ! Scale by gas surface density and dust-to-gas ratio
    ! ------------------------------------------------------------------
    Sig_tri(i, 1) = d2gRatio * gas_Sigma(i) * S0
    Sig_tri(i, 2) = d2gRatio * gas_Sigma(i) * S1

    ! ------------------------------------------------------------------
    ! Apply floor: where Sig_tri <= SigmaFloor, set to 0.1 * SigmaFloor
    ! ------------------------------------------------------------------
    if (Sig_tri(i, 1) <= Sig_floor_tri(i, 1)) Sig_tri(i, 1) = 0.1d0 * Sig_floor_tri(i, 1)
    if (Sig_tri(i, 2) <= Sig_floor_tri(i, 2)) Sig_tri(i, 2) = 0.1d0 * Sig_floor_tri(i, 2)

  end do

end subroutine Sigma_initial




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
    !local variables 
    double precision :: fi_tot(nrad_max+1,Nm_s), fi_diffusive(nrad_max+1,Nm_s),fi_advective(nrad_max+1,Nm_s)
    double precision :: v_gas(nrad_max),v_rad_vrel(nrad_max,Nm_s)
    ! add all the updater funtions here, for example
    ! tripodpy default updater
    !['delta', 'rhos', 'fill', 'backreaction', 'f', 'qrec', 'a', 'm', 'St', 'H', 'rho', 'D', 'eps', 'v', 'p', 'q', 'SigmaFloor', 'S'].
    !print *, "q_rec."
    call calc_q_rec(Sig_tri,a_min_tri,a_max_tri,q_rec,nrad_max)
    !print *, "a."
    call calculate_a(a_min_tri,a_max_tri,q_rec,f_fudge,a_tri,nrad_max,Nm_l)
    !print *, "m."
    call calculate_m(a_tri,rhos_tri,fill_tri,m_tri,nrad_max,Nm_l)
    !print *, "St"
    call st_epstein_stokes1(a_tri,mfp,rhos_tri,Sigma,St_tri,nrad_max,Nm_l)
    !print *, "H_d."
    call h_dubrulle1995(H_gas,St_tri,alpha_vert_tri,H_tri,nrad_max,Nm_l)
    !calculate the midplane density of the dust, which is needed for the collision rates
    !print *, "rho."
    rho_tri = Sig_tri/(sqrt(2.0d0*pi)*H_tri(:,[1,3]))
    !print *, "D."
    call d(alpha_rad_tri*cs**2, OmegaK, St_tri*f_drift, D_tri, nrad_max, Nm_l)
    D_tri(1:2,:) = 0d0 
    D_tri(nrad_max-2+1:,:) = 0d0
    !print *, "vrad."
    !secoind argument is nu 
    call v_visc(Sigma,alpha_rad_tri*cs*H_gas,R,Ri_tri,v_gas,nrad_max)
    v_gas = 0d0
    !second argument is vdieftmax
    call vrad(St_tri*f_drift, -eta*R*OmegaK,v_gas, v_rad_tri, nrad_max, Nm_l)
    call vrad(St_tri, -eta*R*OmegaK,v_gas, v_rad_vrel, nrad_max, Nm_l)
    ! Relative velocities
    !print *, "vrel azi "
    call vrel_azimuthal_drift(eta*R*OmegaK, St_tri, v_rel_azi_tri, nrad_max, Nm_l)
    !print *, "vrel b "
    call vrel_brownian_motion(cs, m_tri, T, v_rel_brown_tri, nrad_max, Nm_l)
    !print *, "vrel dr "
    call vrel_radial_drift(v_rad_vrel, v_rel_rad_tri, nrad_max, Nm_l)
    !print *, "vrel turb "
    call vrel_ormel_cuzzi_2007(alpha_turb_tri, cs, mump, OmegaK, Sigma, St_tri, v_rel_turb_tri, nrad_max, Nm_l)
    !print *, "vrel set "
    call vrel_vertical_settling(H_tri, OmegaK, St_tri, v_rel_vert_tri, nrad_max, Nm_l)
    !print *, "vrel tot "
    v_rel_tot_tri = sqrt(v_rel_azi_tri**2 + v_rel_brown_tri**2 + v_rel_rad_tri**2 + v_rel_turb_tri**2 + v_rel_vert_tri**2)

    ! collision outcomes p and q 
    !print *, "pfrag "
    call pfrag(v_rel_tot_tri(:,4,5), v_frag, p_frag_tri, nrad_max, Nm_l)
    !print *, "pfragtrans ",shape(p_frag_tri),(v_rel_rad_tri(2,4,5))
    call pfrag_trans(St_tri(:,Nm_l), alpha_turb_tri, Sigma, mump, p_fragtrans, nrad_max)
    !print *, "pdriftfrag "
    call pdriftfrag(v_rel_rad_tri(:,4,5),v_rel_azi_tri(:,4,5),St_tri(:,Nm_l),alpha_rad_tri,Sigma,mump,cs,&
                    p_fragtrans,p_drfr_tri,nrad_max)
    !print *, "q_frag "
    call qfrag(p_drfr_tri,v_rel_tot_tri(:,4,5),v_frag,St_tri(:,Nm_l),q_turb1,q_turb2,q_drfr,alpha_turb_tri,Sigma,mump,q_frag_tri,nrad_max)
    q_eff_tri = q_frag_tri*p_frag_tri + q_sweep_tri*(1.0d0 - p_frag_tri)

    ! perarator sterp in tripodpy -> set the state vector rhs 
    !print *, "rh "
    rhs(1:nrad_max*Nm_s) = reshape(transpose(Sig_tri), [nrad_max*Nm_s])
    rhs((nrad_max*Nm_s)+1:(nrad_max*Nm_s)+nrad_max) = a_max_tri*Sig_tri(:,2)
    !print *, "rhs ", rhs((nrad_max*Nm_s)+nrad_max)
    !print *,"vrel", v_rel_tot_tri(:5,4,5),v_rel_azi_tri(:5,4,5),v_rel_rad_tri(:5,4,5),v_rel_turb_tri(:5,4,5),v_rel_vert_tri(:5,4,5),v_rel_brown_tri(:5,4,5)
    call smax_deriv(v_rel_tot_tri(:,4,5),rho_tri(:,2),rhos_tri(:,3), a_min_tri, a_max_tri,v_frag,Sig_tri,Sig_floor_tri,deriv_s_max,nrad_max,Nm_s)
    S_rhs = 0.0d0
    !print *, " scoag"
    call s_coag(pi*(a_tri(:,[1,3])+a_tri(:,[3,2]))**2d0,v_rel_tot_tri(:,[1,3],[3,2]),H_tri(:,[1,3]),m_tri(:,[1,3]),Sig_tri,a_min_tri,a_max_tri,q_eff_tri,Sig_floor_tri,S_coag_tri,nrad_max,Nm_s)
    !print *, "fi "
    call fi_adv(Sig_tri,v_rad_tri(:,[1,3]),R,Ri_tri,fi_advective,nrad_max,Nm_s)
    !print *, "fi diff"
    call fi_diff(D_tri(:,[1,3]),Sig_tri,Sigma,St_tri(:,[1,3])*f_drift,sqrt(alpha_rad_tri*cs**2),R,Ri_tri,fi_diffusive,nrad_max,Nm_s)
    fi_tot = fi_diffusive + fi_advective
    !print *, "s_hyd "
    call s_hyd(Fi_tot,Ri_tri,S_hyd_tri,nrad_max,Nm_s)
    S_tot_tri = S_coag_tri + S_hyd_tri
    call def_smax_hyd(smax_dot_hyd,Sigma,cs,R,Ri_tri)
    !print * ,"update complete"
    call enforce_f()
    !write(*,*) "test_quantity",Ri_tri(:10)
    !print *, "gas Sigma", Sigma(nrad_max-3:)
    !print *, "gas eta ", eta(nrad_max-3:)
    !print *, "gas T ",T(nrad_max-3:)
    !print *, "gas mu ",mump(nrad_max-3:)
    !print *, "gas mfp ",mfp(nrad_max-3:)
    !print *, "gas cs ",cs(nrad_max-3:)
    !print *, "gas H ",H_gas(nrad_max-3:)
    !print *, "D eta ",D_tri(nrad_max-3:,3)
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
    double precision :: Di,K1,K2
    integer :: row_in(Nm_s*3), col_in(Nm_s*3), row_out(Nm_s*3), col_out(Nm_s*3)

    double precision, allocatable :: dat_hydro(:)
    integer, allocatable :: row_hydro(:), col_hydro(:)
    integer :: i    




    cross_section_tri = pi*(a_tri(:,[1,3])+a_tri(:,[3,2]))**2d0
    call jacobian_coagulation_generator(cross_section_tri,v_rel_tot_tri(:,[1,3],[3,2]),H_tri(:,[1,3]),m_tri(:,[1,3]),Sig_tri,a_min_tri,a_max_tri,q_eff_tri,&
                                        dat_coag,row_coag,col_coag,nrad_max,Nm_s)
    !Fortan idexing the arrays wer constucted with python sttyle indexing in mind 
    row_coag = row_coag + 1
    col_coag = col_coag + 1

    !construct the jacobian for the coagulation part and unravel the arrays C like in the pyhton version
    call jacobian_hydrodynamic_generator(area,D_tri(:,[1,3]),R,Ri,Sigma,v_rad_tri(:,[1,3]),A,B,C,nrad_max,Nm_s)

    !transpose the arrays first to match the C style ordering of indices and then reshape them to 1D arrays
    allocate(dat_hydro(3*N_tot -2*Nm_s),row_hydro(3*N_tot -2*Nm_s),col_hydro(3*N_tot -2*Nm_s))
    dat_hydro = [RESHAPE(Transpose(A(2:,:)), [N_tot - Nm_s]), RESHAPE(Transpose(B), [N_tot]), RESHAPE(Transpose(C(1:nrad_max-1,:)), [N_tot - Nm_s])]
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
    if(bd_inner_type .eq. "val")then
      rhs(1:nm_s) = inner_bc
    elseif(bd_inner_type .eq. "const_grad")then 
      Di = ri(2) / ri(3) * (r(2)- r(1)) / (r(3) - r(1))
      K1 = - r(2) / r(1) * (1. + Di)
      K2 = r(3) / r(1) * Di
      dat_in(Nm_s+1:2*Nm_s) = -K1/dt
      dat_in(2*Nm_s+1:) = -K2/dt
      rhs(:Nm_s) = 0
    endif 
    rhs(N_tot-Nm_s+1:N_tot) = outer_bc


    dat_tot = [dat_hydro, dat_coag,dat_in, dat_out]
    row_tot = [row_hydro, row_coag, row_in, row_out]
    col_tot = [col_hydro, col_coag, col_in, col_out]


end subroutine Jacobian

subroutine Y_jacobian(area,R,Ri,Sigma,dt,values_J_out,rowind_J_out,colptr_J_out)    
  use csc_builder, only: coo_to_csc
  implicit none

    double precision, intent(in) :: area(nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: dt
    double precision, intent(out), allocatable :: values_J_out(:)
    integer, intent(out), allocatable :: rowind_J_out(:) 
    integer, intent(out) ::  colptr_J_out(nrad_max *(Nm_s+1)+1)
    ! Local variables for the Jacobian construction
    integer, parameter :: N_tot = int(nrad_max*Nm_s)
    double precision, dimension(nrad_max) :: A,B,C
    double precision, dimension(3) :: dat_in,dat_out
    integer, dimension(3) :: row_in,col_in,row_out,col_out
    integer :: col_diag(nrad_max*nm_s+nrad_max),row_diag(nrad_max*nm_s+nrad_max)
    double precision :: dat_diag(nrad_max*nm_s+nrad_max)
    double precision :: Di,K1,K2
    integer :: i,k,nnz_diag,nzz_new
    logical :: found

    double precision, allocatable :: dat_J(:),dat_total(:),dat_hydro(:)
    integer, allocatable :: row_J(:), col_J(:),row_total(:), col_total(:),row_hydro(:), col_hydro(:)
    double precision, allocatable :: values_J(:)
    integer, allocatable :: rowind_J(:)
    integer :: colptr_J(3*Nrad_max+1)
  
  

    !get base JAcobian
    print *, "calling jacobian"
    !print *, "convert re jac",allocated(values_J_out),allocated(values_J)
    call Jacobian(Sigma,R,Ri,area,dt,dat_J,row_J,col_J)
    !print *, "convert jac",allocated(values_J_out),allocated(values_J)
    !print *, "calling jac_hyd_smax", shape(area),shape(D_tri(:,3)),v_rad_tri(2,3),shape(A),shape(Ri)
    call jacobian_hydrodynamic_generator(area,D_tri(:,3),R,Ri,Sigma,v_rad_tri(:,3),A,B,C,nrad_max,1)
    !print *, "convert hdy jac",allocated(values_J_out),allocated(values_J)
    dat_hydro = [A(2:nrad_max), B(:), C(1:nrad_max-1)]
    row_hydro = [(i+1, i=1,nrad_max-1), (i, i=1,nrad_max), (i, i=1,nrad_max-1)] + N_tot
    col_hydro = [(i, i=1,nrad_max-1), (i, i=1,nrad_max), (i+1, i=1,nrad_max-1)] + N_tot

    !print *, size(dat_hydro),size(dat_J)
    !boudary arrays
    row_in = 1 + N_tot
    col_in = [(i, i=1,3)] + N_tot
    dat_in = 0.0d0
    row_out = nrad_max + N_tot
    col_out = [(nrad_max-3+i, i=1,3)] + N_tot
    dat_out = 0.0d0

    !set smax voundaries here 
    !print *, "look here",s_bd_inner_type 
    if (s_bd_inner_type .eq. "val") then
      rhs(N_tot+1) = inner_s_bc*inner_bc(2)
    elseif(s_bd_inner_type .eq. "const_grad")then
      !print *, "ping --------------------"
      Di = ri(2) / ri(3) * (r(2)- r(1)) / (r(3) - r(1))
      K1 = - r(2) / r(1) * (1. + Di)
      K2 = r(3) / r(1) * Di
      dat_in(2) = -K1/dt
      dat_in(3) = -K2/dt
      rhs(N_tot+1) = 0
    endif 
    rhs(N_tot+nrad_max)=outer_s_bc*outer_bc(2)
    
    dat_total = [dat_J,dat_hydro, dat_in, dat_out]
    row_total = [row_J,row_hydro, row_in, row_out]
    col_total = [col_J,col_hydro, col_in, col_out]
    !print *, size(dat_total)
    deallocate(dat_J, row_J, col_J, dat_hydro, row_hydro, col_hydro)


    dat_total = dat_total*(-dt)
    ! make the actual integration matrix by subtracting eye *dt 
    nnz_diag = 0
    do i = 1, Nrad_max*Nm_s+nrad_max
        found = .false.
        do k = 1, size(dat_total)
            if (row_total(k) == i .and. col_total(k) == i) then
                dat_total(k) = 1d0  + dat_total(k) 
                found = .true.
                exit
            end if
        end do

        if (.not. found) then
            nnz_diag = nnz_diag + 1
            row_diag(nnz_diag) = i
            col_diag(nnz_diag) = i
            dat_diag(nnz_diag) = 1
        end if
    end do

    if (nnz_diag > 0) then
        !print *, "t---------------------------"
        dat_total = [dat_total, dat_diag(1:nnz_diag)]
        row_total = [row_total, row_diag(1:nnz_diag)]
        col_total = [col_total, col_diag(1:nnz_diag)]
    end if
    !print *, size(dat_total)

    !allocate(values_J(SIZE(dat_total)), rowind_J(SIZE(dat_total)), colptr_J(N_tot+nrad_max+1))


    call coo_to_csc(N_tot+nrad_max,N_tot+nrad_max, SIZE(dat_total), dat_total,row_total, col_total,&
                        colptr_J_out, rowind_J_out, values_J_out,nzz_new)

    deallocate(dat_total,row_total,col_total)
    !print *, "end_of Y"



end subroutine Y_jacobian

 
subroutine integrate_dust(area,R,Ri,Sigma,dt)
    implicit none

    double precision, intent(in) :: area(nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: dt

    double precision, allocatable :: values_J(:)
    integer, allocatable :: rowind_J(:)
    integer ::  colptr_J(nrad_max*(Nm_s+1)+1)    


    call Y_jacobian(area,R,Ri,Sigma,dt,values_J,rowind_J,colptr_J)

    !print *, "calling 1dsa", rhs(1:nrad_max*Nm_s) + dt * reshape(transpose(S_rhs), [nrad_max*Nm_s])
    !implement the S coag source for the rhs term should work without though
    rhs(1:nrad_max*Nm_s) = rhs(1:nrad_max*Nm_s) + dt * reshape(transpose(S_rhs), [nrad_max*Nm_s])
    rhs((nrad_max*Nm_s)+2:(nrad_max*Nm_s)+nrad_max-1) = rhs((nrad_max*Nm_s)+2:(nrad_max*Nm_s)+nrad_max-1) + dt * ((deriv_s_max(2:nrad_max-1)*Sig_tri(2:nrad_max-1,2))+(a_max_tri(2:nrad_max-1)*(S_rhs(2:nrad_max-1,2)+S_coag_tri(2:nrad_max-1,2))))
    !print *, "calling superlu", rhs(nrad_max*Nm_s -6 : nrad_max*Nm_s+2)
    !call print_csc_subblock(nrad_max*Nm_s-2*Nm_s,nrad_max*Nm_s,nrad_max*Nm_s-2*Nm_s,nrad_max*Nm_s,size(values_J),colptr_J,rowind_J,values_J)
    call solve_superlu(SIZE(values_J), 1, values_J, rowind_J, colptr_J, rhs)
    !print *, "dealovate arrays"
    deallocate(values_J, rowind_J)
    !print *, "finalizer"
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
        !write(*,*) 'SuperLU Factorization succeeded.'
    end if

    ! 2. Solve the linear system
    iopt = 2
    call c_fortran_dgssv(iopt, n, nnz, nrhs, values, rowind, colptr, &
                         b, ldb, factors, info)

    if (info == 0) then
        !write(*,*) 'SuperLU Solve succeeded. First 5 solution elements:'
        !write(*,*) b(1:min(5, n))
    else
        !write(*,*) 'SuperLU Solve failed, info = ', info
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
    where(Sig_tri .le. Sig_floor_tri)
      Sig_tri = Sig_floor_tri *0.1d0
    end where 
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
        write(*,*) b(1:min(5, n)),values(:10),colptr(:3)
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
double precision, dimension(nrad_max) :: dadsig1
integer :: i

delta = f_crit*sum(Sig_tri,dim=2) -Sig_tri(:,2)
do i = 1, nrad_max
    delta(i) = max(0.0d0, delta(i))
end do
call dadsig(a_lim,q_rec,f_crit,a_max_tri,a_min_tri,Sig_tri,dadsig1,nrad_max,Nm_s)
Sig_tri(:,2) = Sig_tri(:,2) + delta
Sig_tri(:,1) = Sig_tri(:,1) - delta

do i = 1, nrad_max
  a_max_tri(i) = max(a_lim,a_max_tri(i) + delta(i)*dadsig1(i))
enddo 

call calc_q_rec(Sig_tri,a_min_tri,a_max_tri,q_rec,nrad_max)

end subroutine enforce_f


subroutine write_output(t,i_output)
    implicit none
    double precision, intent(in) :: t
    integer, intent(inout) :: i_output
    
    integer :: i

    call OPEN_OUTPUT_FILE(5200+i_output, 1, .True., .false., outfile_name, 4, i_output)
    ! This subroutine should handle all the output writing, for example writing the dust surface density and maximum grain size to files for post-processing and visualization. The implementation can be adjusted as needed, for example by using different file formats or adding more output variables.
    do i = 1, nrad_max
        write(5200+i_output,*) &
            i, &
            sqrt(Ri_tri(i+1)*Ri_tri(i)),&
            Sig_tri(i,1),&
            Sig_tri(i,2),&
            a_max_tri(i),&
            q_rec(i),&
            t
    enddo 
    close(5200+i_output)
    i_output = i_output +1 
end subroutine write_output

!!!!
!Timestep stuff
!!

subroutine calc_ts_tri(dt)
    implicit none 

    double precision, intent(out) :: dt 
    !local variables 
    double precision :: dt_sig,dt_smax


    call calc_dt_Sigma(dt_sig)
    call calc_dt_smax(dt_smax)

    dt = min(dt_sig,dt_smax) *cfl_tri
end subroutine calc_ts_tri


subroutine calc_dt_Sigma(dt_out)
  ! Calculates the time step due to changes in Sig_tri.
  !
  ! Parameters
  ! ----------
  ! nrad_max         : number of radial grid cells
  ! Nm_s         : number of size bins
  ! Sig_tri      : dust surface density          (nrad_max, Nm_s)
  ! S_tot      : total source term             (nrad_max, Nm_s)
  ! Sig_floor_tri : floor value for Sigma         (nrad_max, Nm_s)
  ! a_max_tri      : maximum grain size per cell   (nrad_max)
  ! deriv_s_max: coagulation size growth rate  (nrad_max)
  ! smax_dot_hyd : hydro source at a_max_tri         (nrad_max)
  ! f_crit     : critical fragmentation ratio  (scalar)
  ! dsig_da    : d(Sig_tri)/d(a) per radial cell (nrad_max)   [from dsigda()]
  !
  ! Output
  ! ------
  ! dt_out     : minimum time step (scalar)

  implicit none

  real(8),  intent(out) :: dt_out

  ! Local variables
  double precision ::dsig_da(nrad_max)
  real(8)  :: dt(nrad_max, Nm_s)
  real(8)  :: dt_pred
  real(8)  :: Sigma_sum, f, numerator, denominator
  logical  :: mask(nrad_max, Nm_s)
  logical  :: mask2(nrad_max)
  logical  :: any_neg
  integer  :: i, j

  real(8), parameter :: LARGE = 1.0d100



  call dsigda(a_lim,q_rec,f_crit,a_max_tri,a_min_tri,Sig_tri,dsig_da,nrad_max,Nm_s)

  ! ----------------------------------------------------------------
  ! Check if any interior S_tot < 0  (Python: sim.dust.S.tot[1:-1, ...])
  ! ----------------------------------------------------------------
  any_neg = .false.
  do j = 1, Nm_s
    do i = 2, nrad_max - 1        ! interior rows only (1-based; skip 1 and nrad_max)
      if (S_tot_tri(i, j) < 0.0d0) then
        any_neg = .true.
        exit
      end if
    end do
    if (any_neg) exit
  end do

  if (.not. any_neg) then
    dt_out = LARGE
    return
  end if

  ! ----------------------------------------------------------------
  ! Build mask: Sig_tri > Sig_triFloor  .AND.  S_tot < 0
  ! then zero out the first and last rows
  ! ----------------------------------------------------------------
  do j = 1, Nm_s
    do i = 1, nrad_max
      mask(i, j) = (Sig_tri(i, j) > Sig_floor_tri(i, j)) .and. (S_tot_tri(i, j) < 0.0d0)
    end do
  end do
  mask(1,  :) = .false.
  mask(nrad_max, :) = .false.

  ! ----------------------------------------------------------------
  ! Build mask2 (radial vector, length nrad_max)
  ! Python: S_tot[:,1]*Sigma[:,0] - S_tot[:,0]*Sigma[:,1] < 0
  !         AND  f = Sigma[:,1]/Sigma.sum(-1)  < 0.43
  ! Note: Python index 0 -> Fortran index 1, Python 1 -> Fortran 2
  ! ----------------------------------------------------------------
  do i = 1, nrad_max
    Sigma_sum = 0.0d0
    do j = 1, Nm_s
      Sigma_sum = Sigma_sum + Sig_tri(i, j)
    end do

    if (Sigma_sum > 0.0d0) then
      f = Sig_tri(i, 2) / Sigma_sum
    else
      f = 1.0d0
    end if

    mask2(i) = ( S_tot_tri(i, 2) * Sig_tri(i, 1) - S_tot_tri(i, 1) * Sig_tri(i, 2) < 0.0d0 ) &
               .and. ( f < 0.43d0 )
  end do
  ! Boundaries cannot be active
  mask2(1)  = .false.
  mask2(nrad_max) = .false.

  ! ----------------------------------------------------------------
  ! Initialise dt to LARGE everywhere
  ! ----------------------------------------------------------------
  dt(:, :) = LARGE

  ! ----------------------------------------------------------------
  ! dt[mask] = |Sig_tri[mask] / S_tot[mask]|
  ! ----------------------------------------------------------------
  do j = 1, Nm_s
    do i = 1, nrad_max
      if (mask(i, j)) then
        if (S_tot_tri(i, j) /= 0.0d0) then
          dt(i, j) = abs(Sig_tri(i, j) / S_tot_tri(i, j))
        end if
      end if
    end do
  end do

  ! ----------------------------------------------------------------
  ! Predictive time step for mask2 cells (second size bin, j=2)
  !
  ! Python numerator:
  !   10 * ( -0.1*a_max_tri * dsig_da + Sig_tri[:,1] - f_crit*Sig_tri.sum(-1) )
  ! Python denominator:
  !   S_tot[:,1]*(f_crit-1) + f_crit*S_tot[:,0]
  !   + dsig_da*deriv_s_max + smax_dot_hyd*dsig_da
  ! ----------------------------------------------------------------
  do i = 2, nrad_max - 1
    if (mask2(i)) then
      Sigma_sum = 0.0d0
      do j = 1, Nm_s
        Sigma_sum = Sigma_sum + Sig_tri(i, j)
      end do

      numerator   = 10.0d0 * ( -0.1d0 * a_max_tri(i) * dsig_da(i) &
                               + Sig_tri(i, 2)                    &
                               - f_crit * Sigma_sum )

      denominator = S_tot_tri(i, 2) * (f_crit - 1.0d0)             &
                  + f_crit * S_tot_tri(i, 1)                        &
                  + dsig_da(i) * deriv_s_max(i)                 &
                  + smax_dot_hyd(i) * dsig_da(i)

      ! Handle NaN / Inf (zero denominator) the same way Python does
      if (denominator == 0.0d0) then
        dt_pred = LARGE
      else
        dt_pred = abs(numerator / denominator)
        if (dt_pred /= dt_pred) dt_pred = LARGE   ! NaN guard
        if (dt_pred > LARGE)    dt_pred = LARGE   ! Inf guard
      end if

      dt(i, 2) = min(dt(i, 2), dt_pred)
    end if
  end do

  ! ----------------------------------------------------------------
  ! Return global minimum
  ! ----------------------------------------------------------------
  dt_out = LARGE
  do j = 1, Nm_s
    do i = 1, nrad_max
      if (dt(i, j) < dt_out) dt_out = dt(i, j)
    end do
  end do

end subroutine calc_dt_Sigma

subroutine calc_dt_smax(dt_out)
  ! Calculates the time step due to changes in smax.
  ! Change of smax during one integration step bound by smin and maximum growth factor.
  !
  ! Parameters
  ! ----------
  ! nrad_max          : number of radial grid cells
  ! Nm_s          : number of size bins
  ! Sigma       : dust surface density          (nrad_max, Nm_s)
  ! S_tot       : total source term             (nrad_max, Nm_s)
  ! a_max_tri       : maximum grain size per cell   (nrad_max)
  ! deriv_s_max   : coagulation size growth rate  (nrad_max)
  ! smax_dot_hyd: hydro source at a_max_tri         (nrad_max)
  !
  ! Output
  ! ------
  ! dt_out      : minimum time step (scalar)

  implicit none

  real(8), intent(out) :: dt_out

  ! Local variables
  real(8)  :: f, Sigma_sum, smax_dot, dt_i
  logical  :: mask2(nrad_max)
  integer  :: i, j

  real(8), parameter :: LARGE = 1.0d100
  real(8), parameter :: EPS   = 1.0d-100

  ! ----------------------------------------------------------------
  ! Build mask2 (radial vector):
  ! Python: S_tot[:,1]*Sigma[:,0] - S_tot[:,0]*Sigma[:,1] < 0
  !         AND  f = Sigma[:,1] / Sigma.sum(-1)  < 0.43
  ! Note: Python index 0 -> Fortran 1, Python index 1 -> Fortran 2
  ! ----------------------------------------------------------------
  do i = 1, nrad_max
    Sigma_sum = 0.0d0
    do j = 1, Nm_s
      Sigma_sum = Sigma_sum + Sig_tri(i, j)
    end do

    if (Sigma_sum > 0.0d0) then
      f = Sig_tri(i, 2) / Sigma_sum
    else
      f = 0.0d0
    end if

    mask2(i) = ( S_tot_tri(i, 2) * Sig_tri(i, 1) - S_tot_tri(i, 1) * Sig_tri(i, 2) < 0.0d0 ) &
               .and. ( f < 0.43 )
  end do

  ! ----------------------------------------------------------------
  ! Loop over interior cells only: Python [1:-1] -> Fortran i=2,nrad_max-1
  !
  ! smax_dot = min( |deriv_s_max[i]| , |deriv_s_max[i] + smax_dot_hyd[i]| )
  ! dt[i]    = a_max_tri[i] / (smax_dot + 1e-100)
  ! dt[i]    = LARGE  if mask2[i] is true
  ! ----------------------------------------------------------------
  dt_out = LARGE

  do i = 2, nrad_max - 1
    if (mask2(i)) then
      ! This cell is excluded — keep LARGE
      cycle
    end if

    smax_dot = min( abs(deriv_s_max(i)), abs(deriv_s_max(i) + smax_dot_hyd(i)) )
    dt_i     = a_max_tri(i) / (smax_dot + EPS)

    if (dt_i < dt_out) dt_out = dt_i
  end do

end subroutine calc_dt_smax


subroutine def_smax_hyd(s_smax_hyd,Sigma,cs,R,Ri)

    implicit none 

    double precision, intent(out) :: s_smax_hyd(nrad_max)
    double precision, intent(in) :: Sigma(nrad_max)
    double precision, intent(in) :: cs (nrad_max)
    double precision, intent(in) :: R(nrad_max)
    double precision, intent(in) :: Ri(nrad_max+1)
    !local variables 
    double precision :: s_temp_hyd(nrad_max,Nm_s)
    double precision :: Fi_tot(nrad_max+1,Nm_s),Fi_adv_int(nrad_max+1,Nm_s),Fi_diff_int(nrad_max+1,Nm_s)
    double precision :: Sig_temp(nrad_max,Nm_s)

    Sig_temp = Sig_tri * spread(a_max_tri, dim=2, ncopies=Nm_s)

    call fi_diff(D_tri(:,[1,3]),Sig_temp,Sigma,St_tri(:,[1,3])*f_drift,sqrt(alpha_rad_tri * cs ** 2),&
                R,Ri,Fi_diff_int,nrad_max,Nm_s)

    call fi_adv(Sig_temp, v_rad_tri(:,[1,3]),R,Ri,Fi_adv_int,nrad_max,Nm_s)
    Fi_tot = Fi_adv_int + Fi_diff_int


    call s_hyd(Fi_tot,Ri,s_temp_hyd,nrad_max, Nm_s)
    s_smax_hyd = (s_temp_hyd(:,2) - S_hyd_tri(:,2)*a_max_tri)/Sig_tri(:,2)


end subroutine

!calculate the Ri grid from R since planete does not calculate these by default

subroutine log_grid_interfaces(nrad_max, r_mid)
  ! Calculates the cell interfaces of a logarithmic radial grid from
  ! the cell centres. The interface array includes both outer boundaries,
  ! so it has size nrad_max + 1.
  !
  ! Interior interfaces are placed at the geometric mean of adjacent
  ! centres. The inner and outer boundary interfaces are extrapolated
  ! by the same log-spacing as the first and last cell pair respectively.
  !
  ! Parameters
  ! ----------
  ! nrad_max : number of grid cells
  ! r_mid    : cell centre radii, size (nrad_max)
  !
  ! Output
  ! ------
  ! r_int    : cell interface radii, size (nrad_max + 1)

  implicit none

  integer, intent(in)  :: nrad_max
  real(8), intent(in)  :: r_mid(nrad_max)
  integer :: i

  ! Interior interfaces: geometric mean of adjacent cell centres
  do i = 2, nrad_max
    Ri_tri(i) = sqrt(r_mid(i - 1) * r_mid(i))
  end do

  ! Inner boundary: extrapolate inward by the same log-ratio as
  ! the first cell pair  ->  r_int(1) = r_mid(1)^2 / r_int(2)
   Ri_tri(1) = r_mid(1) * r_mid(1) / Ri_tri(2)

  ! Outer boundary: extrapolate outward by the same log-ratio as
  ! the last cell pair  ->  Ri_tri(nrad_max+1) = r_mid(nrad_max)^2 / Ri_tri(nrad_max)
  Ri_tri(nrad_max+1) = r_mid(nrad_max) * r_mid(nrad_max) / Ri_tri(nrad_max)

end subroutine log_grid_interfaces


subroutine read_static_gas_disk(fname, nrows,ncols,R,OmegaK,Sigma,cs,H_gas,T,mump,mfp,eta,P)
  implicit none

  character(len=*), intent(in)  :: fname
  integer,          intent(in)  :: nrows, ncols
  double precision, intent(out) :: R(nrows)
  double precision, intent(out) :: OmegaK(nrows)
  double precision, intent(out) :: Sigma(nrows)
  double precision, intent(out) :: cs(nrows)
  double precision, intent(out) :: H_gas(nrows)
  double precision, intent(out) :: T(nrows)
  double precision, intent(out) :: mump(nrows)
  double precision, intent(out) :: mfp(nrows)
  double precision, intent(out) :: eta(nrows)
  double precision, intent(out) :: P(nrows)


  !local variables
  integer :: ierr
  real(8) :: data_out(nrows, ncols)

  call read_csv(fname,nrows,ncols,data_out,ierr)

  !fill the relevant arrays
  R = data_out(:,1) !grid 
  OmegaK = data_out(:,2) !
  Sigma = data_out(:,3)
  cs = data_out(:,4)
  H_gas = data_out(:,5)
  T = data_out(:,6)
  mump = data_out(:,7)
  mfp = data_out(:,8)
  eta = data_out(:,9)
  P = data_out(:,10)

end subroutine

subroutine print_csc_subblock(row_start, row_end, col_start, col_end, &
                               ncol, col_ptr, row_ind, csc_val)
  use iso_fortran_env, only: dp => real64, ip => int32
  implicit none
  integer(ip), intent(in) :: row_start, row_end     ! e.g. N-M, N
  integer(ip), intent(in) :: col_start, col_end     ! e.g. 1, N
  integer(ip), intent(in) :: ncol
  integer(ip), intent(in) :: col_ptr(ncol+1)
  integer(ip), intent(in) :: row_ind(:)
  real(dp),    intent(in) :: csc_val(:)

  integer(ip) :: nrows, ncols, i_col, k
  real(dp), allocatable :: dense(:,:)

  nrows = row_end   - row_start + 1
  ncols = col_end   - col_start + 1
  allocate(dense(nrows, ncols))
  dense = 0.0_dp

  do i_col = col_start, min(col_end, ncol)
    do k = col_ptr(i_col), col_ptr(i_col+1) - 1
      if (row_ind(k) >= row_start .and. row_ind(k) <= row_end) then
        dense(row_ind(k) - row_start + 1, &
              i_col       - col_start + 1) = csc_val(k)
      end if
    end do
  end do

  write(*,'(A,I0,A,I0)') "  rows ", row_start, " to ", row_end
  do k = 1, nrows
    write(*,'(*(ES10.3,1X))') dense(k, :)
  end do

  deallocate(dense)
end subroutine

end module tripod