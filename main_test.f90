program main
    use tripod
    implicit none

    print *, "Calling solve_dummy_superlu..."

    call solve_dummy_superlu()

    print *, "Done."
end program main