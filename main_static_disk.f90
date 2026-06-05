program main
    use tripod
    use unitsPL
    implicit None

    !gas variables 
    integer, parameter :: nrows = 100
    double precision :: R(nrows)
    double precision :: OmegaK(nrows)
    double precision :: Sigma(nrows)
    double precision :: cs(nrows)
    double precision::  H_gas(nrows)
    double precision :: T(nrows)
    double precision :: mump(nrows)
    double precision::  mfp(nrows)
    double precision :: eta(nrows)
    double precision :: P(nrows)
    double precision :: area(nrows)

    !parameters for initialisation
    double precision,parameter :: a_min_ini = 5.22875516e-05
    double precision,parameter :: a_max_ini = 1e-4
    double precision,parameter :: alpha = 1e-3
    double precision, parameter :: fd2g = 1e-2
    double precision, parameter :: rhos = 1.67

    !additional bits and bobs
    double precision :: time
    double precision :: timestep,timestep_lim

    !Timestep and snapshot variables
    integer, parameter :: nsnaps = 100
    double precision,parameter :: t_min = 1d2
    double precision,parameter :: t_max = 2d5
    double precision, parameter :: cfl = 1d-1

    real(8) :: snaps(nsnaps)
    integer :: i_output, i_ts
    !
    integer :: i

    snaps = 10.0d0**( log10(t_min) + &
        [(real(i-1,8), i=1,nsnaps)] * &
        (log10(t_max)-log10(t_min))/(nsnaps-1) )
    snaps = snaps * an
    i_output = 1
    i_ts = 1
    time = 0d0

    print *, "starting programm"

    call read_static_gas_disk("test_sim.csv", nrows,10,R,OmegaK,Sigma,cs,H_gas,T,mump,mfp,eta,P)
    call write_output(0d0,i_output)
    call log_grid_interfaces(nrows,R)
    area =pi*(Ri_tri(2:)**2 - Ri_tri(:nrows)**2)
    print *, "1.", nrows .eq. nrad_max
    call initialize_dust(a_min_ini,a_max_ini,alpha,alpha,alpha,fd2g,rhos,Sigma,P,cs,T,H_gas,mump,eta,mfp,OmegaK,R)
    print *, "read the gas disk"
    call write_output(0d0,i_output)
    

    do while(.true.)
        print *, "update",time/an
        call update_dust(R,eta,T,mump,OmegaK,mfp,Sigma,cs,H_gas)
        print *,"wub",deriv_s_max(:5)
        if(any(abs(snaps-time) .lt. epsilon(time)))then 
            print *, "write at",time/an
            call write_output(time,i_output)
        endif 
        print *, "calc_ts",time/an
        call calc_ts_tri(timestep)
        call limit_ts_to_snaps(time,timestep,timestep_lim,snaps,nsnaps)
        !timestep_lim = timestep_lim * cfl
        print *, "integrate_dust",time/an,timestep_lim/an, timestep/an
        !timestep_lim = 1e-1 *an
        call  integrate_dust(area,R,Ri_tri,Sigma,timestep_lim)
        print *,"wub",deriv_s_max(:5)
        time = time + timestep_lim
        i_ts = i_ts +1 
        if(.True.)then 
            if(time > 1e5*an .or. .false.)then 
                call write_output(time,i_output)
                stop
            endif
        endif 
    enddo


end program

subroutine limit_ts_to_snaps(t,dt,dt_lim,snaps,nsnaps)
      implicit none

  real(8), intent(in)  :: t,dt
  double precision, intent(out) :: dt_lim
  integer, intent(in) :: nsnaps
  real(8), intent(in)  :: snaps(nsnaps)

  integer :: i
  real(8) :: target

  target = t + dt

  do i = 1, nsnaps
    if (snaps(i) - t > epsilon(t)) then
      exit
    end if
  end do

  dt_lim = min(dt,snaps(i) - t)

end subroutine limit_ts_to_snaps


