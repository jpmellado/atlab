module TLab_Arrays
    use TLab_Constants, only: wp
    implicit none
    save

    real(wp), allocatable :: q(:, :)        ! Eulerian fields, flow vartiables
    real(wp), allocatable :: s(:, :)        ! Eulerian fields, scalar variables
    real(wp), allocatable :: txc(:, :)      ! Temporary space for Eulerian fields
    real(wp), allocatable :: wrk1d(:, :)    ! Work arrays (scratch space)
    real(wp), allocatable :: wrk2d(:, :)    ! Work arrays (scratch space)
    real(wp), allocatable :: wrk3d(:)       ! Work arrays (scratch space)

    target q, s, txc, wrk1d, wrk2d, wrk3d

end module TLab_Arrays
