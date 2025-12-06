!########################################################################
!# HISTORY
!#
!# 2011/11/01 - C. Ansorge
!#              Created
!#
!########################################################################
!#
!# Transposition using Cache-Blocking and OpenMP
!# transpose of the matrix a and place the transposed matrix in b
!# routine trans below is faster than TRANSPOSE routine from f90
!#
!########################################################################
subroutine TLab_Transpose(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    real(wp), intent(in) :: a(ma, *) ! Input array
    real(wp), intent(out) :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
#ifdef HLRS_HAWK
    parameter(jb=16, kb=8)
#else
    parameter(jb=64, kb=64)
#endif

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
#ifdef USE_MKL
    call MKL_DOMATCOPY('c', 't', nra, nca, 1.0_wp, a, ma, b, mb)
#else
    !use own implementation

    kk = 1; jj = 1

    do k = 1, nca - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, nca
        do j = 1, nra
            b(k, j) = a(j, k)
        end do
    end do

    ! do k = 1, nca
    do k = 1, last_k - 1
        do j = last_j, nra
            b(k, j) = a(j, k)
        end do
    end do

#endif

    return
end subroutine TLab_Transpose

!########################################################################
!########################################################################
subroutine TLab_AddTranspose(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    real(wp), intent(in) :: a(ma, *)    ! Input array
    real(wp), intent(inout) :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
#ifdef HLRS_HAWK
    parameter(jb=16, kb=8)
#else
    parameter(jb=64, kb=64)
#endif

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
    kk = 1; jj = 1

    do k = 1, nca - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = b(kk, jj) + a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, nca
        do j = 1, nra
            b(k, j) = b(k, j) + a(j, k)
        end do
    end do

    ! do k = 1, nca
    do k = 1, last_k - 1
        do j = last_j, nra
            b(k, j) = b(k, j) + a(j, k)
        end do
    end do

    return
end subroutine TLab_AddTranspose

!########################################################################
!########################################################################
subroutine TLab_SubtractTranspose(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    real(wp), intent(in) :: a(ma, *)    ! Input array
    real(wp), intent(inout) :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
#ifdef HLRS_HAWK
    parameter(jb=16, kb=8)
#else
    parameter(jb=64, kb=64)
#endif

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
    kk = 1; jj = 1

    do k = 1, nca - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = b(kk, jj) - a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, nca
        do j = 1, nra
            b(k, j) = b(k, j) - a(j, k)
        end do
    end do

    ! do k = 1, nca
    do k = 1, last_k - 1
        do j = last_j, nra
            b(k, j) = b(k, j) - a(j, k)
        end do
    end do

    return
end subroutine TLab_SubtractTranspose

!########################################################################
!########################################################################
subroutine TLab_Transpose_INT1(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    integer(1), intent(in) :: a(ma, *) ! Input array
    integer(1), intent(out) :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
    parameter(jb=32, kb=32)

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
    kk = 1; jj = 1

    do k = 1, nca - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, nca
        do j = 1, nra
            b(k, j) = a(j, k)
        end do
    end do

    ! do k = 1, nca
    do k = 1, last_k - 1
        do j = last_j, nra
            b(k, j) = a(j, k)
        end do
    end do

    return
end subroutine TLab_Transpose_INT1

!########################################################################
!########################################################################
subroutine TLab_Transpose_COMPLEX(a, nra, nca, ma, b, mb)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(in) :: nra      ! Number of rows in a
    integer(wi), intent(in) :: nca      ! Number of columns in b
    integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    complex(wp), intent(in) :: a(ma, *) ! Input array
    complex(wp), intent(out) :: b(mb, *) ! Transposed array

! -------------------------------------------------------------------
    integer(wi) jb, kb
#ifdef HLRS_HAWK
    parameter(jb=16, kb=8)
#else
    parameter(jb=64, kb=64)
#endif

    integer(wi) k, j, jj, kk
    integer(wi) last_k, last_j

! -------------------------------------------------------------------
    kk = 1; jj = 1

    do k = 1, nca - kb + 1, kb; 
        do j = 1, nra - jb + 1, jb; 
            do jj = j, j + jb - 1
                do kk = k, k + kb - 1
                    b(kk, jj) = a(jj, kk)
                end do
            end do
        end do
    end do

    last_k = kk
    last_j = jj

    do k = last_k, nca
        do j = 1, nra
            b(k, j) = a(j, k)
        end do
    end do

    ! do k = 1, nca
    do k = 1, last_k - 1
        do j = last_j, nra
            b(k, j) = a(j, k)
        end do
    end do

    return
end subroutine TLab_Transpose_COMPLEX
