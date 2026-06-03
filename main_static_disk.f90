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

    !parameters for initialisation
    double precision,parameter :: a_min_ini = 5.22875516e-05
    double precision,parameter :: a_max_ini = 1e-4
    double precision,parameter :: alpha = 1e-3
    double precision, parameter :: fd2g = 1e-2
    double precision, parameter :: rhos = 1.67

    !Timestep and snapshot variables
    integer, parameter :: nsnaps = 100
    real(8) :: snaps(nsnaps)
    !
    integer :: i

    snaps = 10.0d0**( log10(1.0d-3) + &
        [(real(i-1,8), i=1,nsnaps)] * &
        (log10(1.0d3)-log10(1.0d-3))/(nsnaps-1) )

    print *, "starting programm"

    call read_static_gas_disk("test_sim.csv", nrows,10,R,OmegaK,Sigma,cs,H_gas,T,mump,mfp,eta,P)

    call log_grid_interfaces(nrows,R,Ri_tri)
    call initialize_dust(a_min_ini,a_max_ini,alpha,alpha,alpha,fd2g,rhos,Sigma,P,cs,T,H_gas,mump,eta,mfp,OmegaK,R)
    print *, "read the gas disk"
    write(*,*)R(:5)/AU,Ri_tri(:5)/AU
    print *, "Done."
    call write_output(2,0d0)
end program