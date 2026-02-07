! Generalized vectorized Thomas algorithms for linear systems

module Thomas
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: Thomas_FactorLU_InPlace
    public :: Thomas_SolveL, Thomas_SolveU

    public :: Thomas3_FactorLU_InPlace
    public :: Thomas3_SolveL, Thomas3_SolveU  ! Particularized for tridiagonal systems

    public :: Thomas5_FactorLU_InPlace
    public :: Thomas5_SolveL, Thomas5_SolveU  ! Particularized for pentadiagonal systems

    public :: Thomas7_FactorLU_InPlace
    public :: Thomas7_SolveL, Thomas7_SolveU  ! Particularized for heptadiagonal systems

    public :: thomas_dt

    ! -----------------------------------------------------------------------
    type :: thomas_dt
        real(wp), allocatable :: L(:, :)
        real(wp), allocatable :: U(:, :)
        procedure(thomas_ice), pointer, nopass :: ptr_solveL
        procedure(thomas_ice), pointer, nopass :: ptr_solveU
    contains
        procedure :: initialize => thomas_initialize_dt
        procedure :: solveL => thomas_solveL_dt
        procedure :: solveU => thomas_solveU_dt
    end type

    abstract interface
        subroutine thomas_ice(L, f)
            import wp
            real(wp), intent(in) :: L(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine thomas_ice
    end interface

contains
    ! #######################################################################
    ! #######################################################################
    subroutine thomas_initialize_dt(self, lhs)
        class(thomas_dt), intent(out) :: self
        real(wp), intent(in) :: lhs(:, :)

        integer ndl
        real(wp), allocatable :: lu_aux(:, :)

        allocate (lu_aux, source=lhs)
        ndl = size(lu_aux, 2)
        call Thomas_FactorLU_InPlace(lu_aux(:, 1:ndl/2), &
                                     lu_aux(:, ndl/2 + 1:ndl))
        allocate (self%L, source=lu_aux(:, 1:ndl/2))
        allocate (self%U, source=lu_aux(:, ndl/2 + 1:ndl))
        deallocate (lu_aux)

        select case (ndl)
        case (3)
            self%ptr_solveL => Thomas3_SolveL
            self%ptr_solveU => Thomas3_SolveU
        case (5)
            self%ptr_solveL => Thomas5_SolveL
            self%ptr_solveU => Thomas5_SolveU
        case (7)
            self%ptr_solveL => Thomas7_SolveL
            self%ptr_solveU => Thomas7_SolveU
        end select

        return
    end subroutine

    subroutine thomas_solveL_dt(self, f)
        class(thomas_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%ptr_solveL(self%L, f)

        return
    end subroutine

    subroutine thomas_solveU_dt(self, f)
        class(thomas_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%ptr_solveU(self%U, f)

        return
    end subroutine

    ! #######################################################################
    ! #######################################################################
    ! LU factorization stage; L has diagonal 1
    ! Original array in L (lower diagonals) and U (center and upper diagonals)
    subroutine Thomas_FactorLU_InPlace(L, U)
        real(wp), intent(inout) :: L(:, :), U(:, :)

        ! -----------------------------------------------------------------------
        integer n, nmax
        integer nd, id, ii

        ! #######################################################################
        nmax = size(L, 1)
        nd = size(L, 2)

        do n = 2, nd
            do id = nd - n + 2, nd
                do ii = id - 1, nd - n + 2, -1
                    L(n, id) = L(n, id) - L(n, ii)*U(n - nd + ii - 1, id - ii + 1)
                end do
                L(n, id) = L(n, id)/U(n - nd + id - 1, 1)
            end do

            do id = 1, nd
                do ii = 1, n - 1
                    U(n, id) = U(n, id) - L(n, nd - ii + 1)*U(n - ii, id + ii)
                end do
            end do

        end do

        do n = nd + 1, nmax - nd + 1
            do id = 1, nd
                do ii = 1, id - 1
                    L(n, id) = L(n, id) - L(n, ii)*U(n - nd + ii - 1, id - ii + 1)
                end do
                L(n, id) = L(n, id)/U(n - nd + id - 1, 1)
            end do

            do id = 1, nd
                do ii = 1, nd - id + 1
                    U(n, id) = U(n, id) - L(n, nd - ii + 1)*U(n - ii, id + ii)
                end do
            end do

        end do

        do n = nmax - nd + 2, nmax
            do id = 1, nd
                do ii = 1, id - 1
                    L(n, id) = L(n, id) - L(n, ii)*U(n - nd + ii - 1, id - ii + 1)
                end do
                L(n, id) = L(n, id)/U(n - nd + id - 1, 1)
            end do

            do id = 1, nmax - n + 1
                do ii = 1, nd - id + 1
                    U(n, id) = U(n, id) - L(n, nd - ii + 1)*U(n - ii, id + ii)
                end do
            end do

        end do

        ! Final operations
        do id = 1, nd
            L(nd + 2 - id:, id) = -L(nd + 2 - id:, id)
        end do

        U(:, 1) = 1.0_wp/U(:, 1)
        do id = 1, nd
            U(:nmax - id, id + 1) = -U(:nmax - id, id + 1)*U(:nmax - id, 1)
        end do

        return
    end subroutine Thomas_FactorLU_InPlace

    ! #######################################################################
    ! #######################################################################
    ! Solve Ly=f, forward elimination
    subroutine Thomas_SolveL(L, f)
        real(wp), intent(in) :: L(:, :)
        real(wp), intent(inout) :: f(:, :)          ! RHS and solution

        ! -----------------------------------------------------------------------
        integer n
        integer nd, id

        ! #######################################################################
        if (size(f, 1) <= 0) return

        nd = size(L, 2)

        do n = 2, nd
            do id = 1, n - 1
                f(:, n) = f(:, n) + f(:, n - id)*L(n, nd - id + 1)
            end do
        end do

        do n = nd + 1, size(L, 1)
            do id = 1, nd
                f(:, n) = f(:, n) + f(:, n - id)*L(n, nd - id + 1)
            end do
        end do

        return
    end subroutine Thomas_SolveL

    ! #######################################################################
    ! #######################################################################
    ! Solve Ux=y, backward substitution
    subroutine Thomas_SolveU(U, f)
        real(wp), intent(in) :: U(:, :)
        real(wp), intent(inout) :: f(:, :)      ! RHS and solution

        ! -----------------------------------------------------------------------
        integer n, nmax
        integer nd, id

        ! #######################################################################
        if (size(f, 1) <= 0) return

        nd = size(U, 2)
        nmax = size(U, 1)

        do n = 0, nd - 1
            f(:, nmax - n) = U(nmax - n, 1)*f(:, nmax - n)
            do id = 2, n + 1
                f(:, nmax - n) = f(:, nmax - n) + U(nmax - n, id)*f(:, nmax - n + id - 1)
            end do
        end do

        do n = nmax - nd, 1, -1
            f(:, n) = U(n, 1)*f(:, n)
            do id = 2, nd
                f(:, n) = f(:, n) + U(n, id)*f(:, n + id - 1)
            end do
        end do

        return
    end subroutine Thomas_SolveU

    ! #######################################################################
    ! #######################################################################
    ! LU factorization stage; L has diagonal 1
    ! Original array in L (lower diagonals) and U (center and upper diagonals)
    subroutine Thomas3_FactorLU_InPlace(L, U)
        real(wp), intent(inout) :: L(:, :), U(:, :)

        ! -----------------------------------------------------------------------
        integer n, nmax

        ! #######################################################################
        nmax = size(L, 1)

        do n = 2, nmax
            L(n, 1) = L(n, 1)/U(n - 1, 1)
            U(n, 1) = U(n, 1) - L(n, 1)*U(n - 1, 2)
        end do

        ! Final operations
        L(2:, 1) = -L(2:, 1)
        U(:, 1) = 1.0_wp/U(:, 1)
        U(:nmax - 1, 2) = -U(:nmax - 1, 2)*U(:nmax - 1, 1)

        return
    end subroutine Thomas3_FactorLU_InPlace

    ! #######################################################################
    ! #######################################################################
    ! Solve Ly=f, forward elimination
    subroutine Thomas3_SolveL(L, f)
        real(wp), intent(in) :: L(:, :)
        real(wp), intent(inout) :: f(:, :)          ! RHS and solution

        ! -------------------------------------------------------------------
        integer :: n

        ! ###################################################################
        if (size(f, 1) <= 0) return

        do n = 2, size(L, 1)
            f(:, n) = f(:, n) + L(n, 1)*f(:, n - 1)
        end do

        return
    end subroutine Thomas3_SolveL

    ! #######################################################################
    ! #######################################################################
    ! Solve Ux=y, backward substitution
    subroutine Thomas3_SolveU(U, f)
        real(wp), intent(in) :: U(:, :)
        real(wp), intent(inout) :: f(:, :)      ! RHS and solution

        ! -------------------------------------------------------------------
        integer :: n, nmax

        ! ###################################################################
        if (size(f, 1) <= 0) return

        nmax = size(U, 1)

        n = nmax
        f(:, n) = f(:, n)*U(n, 1)

        do n = nmax - 1, 1, -1
            f(:, n) = U(n, 1)*f(:, n) + U(n, 2)*f(:, n + 1)
        end do

        return
    end subroutine Thomas3_SolveU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas5_FactorLU_InPlace(L, U)
        real(wp), intent(inout) :: L(:, :), U(:, :)

        ! -----------------------------------------------------------------------
        integer n, nmax

        ! #######################################################################
        nmax = size(L, 1)

        n = 2
        L(n, 2) = L(n, 2)/U(n - 1, 1)

        U(n, 1) = U(n, 1) - L(n, 2)*U(n - 1, 2)
        U(n, 2) = U(n, 2) - L(n, 2)*U(n - 1, 3)

        do n = 3, nmax - 1
            L(n, 1) = L(n, 1)/U(n - 2, 1)
            L(n, 2) = (L(n, 2) - L(n, 1)*U(n - 2, 2))/U(n - 1, 1)

            U(n, 1) = U(n, 1) - L(n, 2)*U(n - 1, 2) - L(n, 1)*U(n - 2, 3)
            U(n, 2) = U(n, 2) - L(n, 2)*U(n - 1, 3)

        end do

        n = nmax
        L(n, 1) = L(n, 1)/U(n - 2, 1)
        L(n, 2) = (L(n, 2) - L(n, 1)*U(n - 2, 2))/U(n - 1, 1)

        U(n, 1) = U(n, 1) - L(n, 2)*U(n - 1, 2) - L(n, 1)*U(n - 2, 3)

        ! Final operations
        L(3:, 1) = -L(3:, 1)
        L(2:, 2) = -L(2:, 2)

        U(:, 1) = 1.0_wp/U(:, 1)
        U(:nmax - 1, 2) = -U(:nmax - 1, 2)*U(:nmax - 1, 1)
        U(:nmax - 2, 3) = -U(:nmax - 2, 3)*U(:nmax - 2, 1)

        return
    end subroutine Thomas5_FactorLU_InPlace

    ! #######################################################################
    ! #######################################################################
    ! Solve Ly=f, forward elimination
    subroutine Thomas5_SolveL(L, f)
        real(wp), intent(in) :: L(:, :)
        real(wp), intent(inout) :: f(:, :)          ! RHS and solution

        ! -----------------------------------------------------------------------
        integer n

        ! #######################################################################
        if (size(f, 1) <= 0) return

        n = 2
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2)

        do n = 3, size(L, 1)
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1)
        end do

        return
    end subroutine Thomas5_SolveL

    ! #######################################################################
    ! #######################################################################
    ! Solve Ux=y, backward substitution
    subroutine Thomas5_SolveU(U, f)
        real(wp), intent(in) :: U(:, :)
        real(wp), intent(inout) :: f(:, :)      ! RHS and solution

        ! -----------------------------------------------------------------------
        integer n, nmax

        ! #######################################################################
        if (size(f, 1) <= 0) return

        nmax = size(U, 1)

        n = nmax
        f(:, n) = f(:, n)*U(n, 1)

        n = nmax - 1
        f(:, n) = U(n, 1)*f(:, n) + U(n, 2)*f(:, n + 1)

        do n = nmax - 2, 1, -1
            f(:, n) = U(n, 1)*f(:, n) + U(n, 2)*f(:, n + 1) + U(n, 3)*f(:, n + 2)
        end do

        return
    end subroutine Thomas5_SolveU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas7_FactorLU_InPlace(L, U)
        real(wp), intent(inout) :: L(:, :), U(:, :)

        ! -----------------------------------------------------------------------
        integer n, nmax

        ! #######################################################################
        nmax = size(L, 1)

        n = 2
        L(n, 3) = L(n, 3)/U(n - 1, 1)

        U(n, 1) = U(n, 1) - L(n, 3)*U(n - 1, 2)
        U(n, 2) = U(n, 2) - L(n, 3)*U(n - 1, 3)
        U(n, 3) = U(n, 3) - L(n, 3)*U(n - 1, 4)

        n = 3
        L(n, 2) = L(n, 2)/U(n - 2, 1)
        L(n, 3) = (L(n, 3) - L(n, 2)*U(n - 2, 2))/U(n - 1, 1)

        U(n, 1) = U(n, 1) - L(n, 3)*U(n - 1, 2) - L(n, 2)*U(n - 2, 3)
        U(n, 2) = U(n, 2) - L(n, 3)*U(n - 1, 3) - L(n, 2)*U(n - 2, 4)
        U(n, 3) = U(n, 3) - L(n, 3)*U(n - 1, 4)

        do n = 4, nmax - 2
            L(n, 1) = L(n, 1)/U(n - 3, 1)
            L(n, 2) = (L(n, 2) - L(n, 1)*U(n - 3, 2))/U(n - 2, 1)
            L(n, 3) = (L(n, 3) - L(n, 1)*U(n - 3, 3) - L(n, 2)*U(n - 2, 2))/U(n - 1, 1)

            U(n, 1) = U(n, 1) - L(n, 3)*U(n - 1, 2) - L(n, 2)*U(n - 2, 3) - L(n, 1)*U(n - 3, 4)
            U(n, 2) = U(n, 2) - L(n, 3)*U(n - 1, 3) - L(n, 2)*U(n - 2, 4)
            U(n, 3) = U(n, 3) - L(n, 3)*U(n - 1, 4)

        end do

        n = nmax - 1
        L(n, 1) = L(n, 1)/U(n - 3, 1)
        L(n, 2) = (L(n, 2) - L(n, 1)*U(n - 3, 2))/U(n - 2, 1)
        L(n, 3) = (L(n, 3) - L(n, 1)*U(n - 3, 3) - L(n, 2)*U(n - 2, 2))/U(n - 1, 1)

        U(n, 1) = U(n, 1) - L(n, 3)*U(n - 1, 2) - L(n, 2)*U(n - 2, 3) - L(n, 1)*U(n - 3, 4)
        U(n, 2) = U(n, 2) - L(n, 3)*U(n - 1, 3) - L(n, 2)*U(n - 2, 4)

        n = nmax
        L(n, 1) = L(n, 1)/U(n - 3, 1)
        L(n, 2) = (L(n, 2) - L(n, 1)*U(n - 3, 2))/U(n - 2, 1)
        L(n, 3) = (L(n, 3) - L(n, 1)*U(n - 3, 3) - L(n, 2)*U(n - 2, 2))/U(n - 1, 1)

        U(n, 1) = U(n, 1) - L(n, 3)*U(n - 1, 2) - L(n, 2)*U(n - 2, 3) - L(n, 1)*U(n - 3, 4)

        ! Final operations
        L(4:, 1) = -L(4:, 1)
        L(3:, 2) = -L(3:, 2)
        L(2:, 3) = -L(2:, 3)

        U(:, 1) = 1.0_wp/U(:, 1)
        U(:nmax - 1, 2) = -U(:nmax - 1, 2)*U(:nmax - 1, 1)
        U(:nmax - 2, 3) = -U(:nmax - 2, 3)*U(:nmax - 2, 1)
        U(:nmax - 3, 4) = -U(:nmax - 3, 4)*U(:nmax - 3, 1)

        return
    end subroutine Thomas7_FactorLU_InPlace

    ! #######################################################################
    ! #######################################################################
    ! Solve Ly=f, forward elimination
    subroutine Thomas7_SolveL(L, f)
        real(wp), intent(in) :: L(:, :)
        real(wp), intent(inout) :: f(:, :)          ! RHS and solution

        ! -----------------------------------------------------------------------
        integer n

        ! #######################################################################
        if (size(f, 1) <= 0) return

        ! n = 1
        ! f(:, n) = f(:, n)*L(n, 3) ! Normalize first eqn. See Thomas7_FactorLU

        n = 2
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 3)

        n = 3
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 3) + f(:, n - 2)*L(n, 2)

        do n = 4, size(L, 1)
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 3) + f(:, n - 2)*L(n, 2) + f(:, n - 3)*L(n, 1)
        end do

        return
    end subroutine Thomas7_SolveL

    ! #######################################################################
    ! #######################################################################
    ! Solve Ux=y, backward substitution
    subroutine Thomas7_SolveU(U, f)
        real(wp), intent(in) :: U(:, :)
        real(wp), intent(inout) :: f(:, :)      ! RHS and solution

        ! -----------------------------------------------------------------------
        integer n, nmax

        ! #######################################################################
        if (size(f, 1) <= 0) return

        nmax = size(U, 1)

        n = nmax
        f(:, n) = f(:, n)*U(n, 1)

        n = nmax - 1
        f(:, n) = U(n, 1)*f(:, n) + U(n, 2)*f(:, n + 1)

        n = nmax - 2
        f(:, n) = U(n, 1)*f(:, n) + U(n, 2)*f(:, n + 1) + U(n, 3)*f(:, n + 2)

        do n = nmax - 3, 1, -1
            f(:, n) = U(n, 1)*f(:, n) + U(n, 2)*f(:, n + 1) + U(n, 3)*f(:, n + 2) + U(n, 4)*f(:, n + 3)
        end do

        return
    end subroutine Thomas7_SolveU

end module Thomas
