!# Transposition using Cache-Blocking and OpenMP
!# transpose of the matrix a and place the transposed matrix in b
module TLab_Transpose
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: TLab_Transpose_Initialize
    public :: TLab_Transpose_Real
    public :: TLab_AddTranspose, TLab_SubtractTranspose
    public :: TLab_Transpose_Complex
    ! public :: TLab_Transpose_INT1

    public :: trans_x, trans_y

    ! -----------------------------------------------------------------------
    ! optimized block sizes; default values
    type :: block_dt
        sequence
#ifdef HLRS_HAWK
        integer :: jb = 16
        integer :: kb = 8
#else
        integer :: jb = 64
        integer :: kb = 64
#endif
    end type block_dt
    type(block_dt) trans_x, trans_y

contains
    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Transpose_Initialize()
        use TLab_Constants, only: lfile, fmt_r
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Pointers, only: tmp1, tmp2
        use TLab_WorkFlow, only: TLab_Write_ASCII

        ! -----------------------------------------------------------------------
        ! Scanning values
#ifdef HLRS_HAWK
        integer :: jb_scan(3) = [8, 16, 32]
        integer :: kb_scan(3) = [4, 8, 16]
#else
        integer :: jb_scan(3) = [32, 64, 128]
        integer :: kb_scan(3) = [32, 64, 128]
#endif
        ! ###################################################################
        if (.not. associated(tmp1)) return
        if (.not. associated(tmp2)) return
        tmp1(:) = 0.0_wp                ! particular value not important

        if (imax > 1) then
            call TLab_Write_ASCII(lfile, 'Optimizing cache-blocking for array transposition along X...')
            call scan(imax, jmax*kmax, trans_x)
        end if

        if (jmax > 1) then
            call TLab_Write_ASCII(lfile, 'Optimizing cache-blocking for array transposition along Y...')
            call scan(imax*jmax, kmax, trans_y)
        end if

        return
    contains
        ! ###################################################################
        subroutine scan(nsize, msize, trans_loc)
            integer(wi) nsize, msize
            type(block_dt) trans_loc

            integer j, k, opt(2)
            integer jb_bak, kb_bak

            real(wp) :: times(size(jb_scan), size(kb_scan))
            real(wp) :: start, finish
            integer :: nsample = 10
            integer n

            character(len=256) line
            character(len=32) str

            ! ###################################################################
            jb_bak = trans_loc%jb             ! backup default values
            kb_bak = trans_loc%kb

            do k = 1, size(kb_scan)
                trans_loc%kb = kb_scan(k)
                line = 'Scanning times'

                do j = 1, size(jb_scan)
                    trans_loc%jb = jb_scan(j)
                    call cpu_time(start)
                    do n = 1, nsample
                        call TLab_Transpose_Real(tmp1, imax, jmax*kmax, imax, tmp2, jmax*kmax)
                    end do
                    call cpu_time(finish)
                    times(j, k) = finish - start
                    write (str, fmt_r) times(j, k)
                    line = trim(adjustl(line))//' '//trim(adjustl(str))
                end do
                call TLab_Write_ASCII(lfile, trim(adjustl(line)))
                ! print *, trim(adjustl(line))
            end do

            opt = minloc(times)
            trans_loc%jb = jb_scan(opt(1))
            trans_loc%kb = kb_scan(opt(2))

            ! ! revert to default values; used in testing
            ! trans_loc%jb = jb_bak
            ! trans_loc%kb = kb_bak

            write (str, *) trans_loc%jb
            line = 'Transposition block '//trim(adjustl(str))
            write (str, *) trans_loc%kb
            line = trim(adjustl(line))//'x'//trim(adjustl(str))//'.'
            call TLab_Write_ASCII(lfile, trim(adjustl(line)))
            ! print *, trim(adjustl(line))

            return
        end subroutine Scan

    end subroutine TLab_Transpose_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Transpose_Real(a, nra, nca, ma, b, mb, locBlock)
        integer(wi), intent(in) :: nra      ! Number of rows in a
        integer(wi), intent(in) :: nca      ! Number of columns in b
        integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
        integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
        real(wp), intent(in) :: a(ma, nca)  ! Input array
        real(wp), intent(out) :: b(mb, nra) ! Transposed array
        type(block_dt), optional :: locBlock

        ! -------------------------------------------------------------------
        integer(wi) jb, kb
        integer(wi) k, j, jj, kk
        integer(wi) last_k, last_j

        ! ###################################################################
#ifdef USE_MKL
        call MKL_DOMATCOPY('c', 't', nra, nca, 1.0_wp, a, ma, b, mb)
#else
        ! ###################################################################
        ! use own implementation
        if (present(locBlock)) then
            jb = locBlock%jb
            kb = locBlock%kb
        else
            jb = trans_x%jb
            kb = trans_x%kb
        end if

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

        do k = 1, last_k - 1
            do j = last_j, nra
                b(k, j) = a(j, k)
            end do
        end do

#endif

        return
    end subroutine TLab_Transpose_Real

    !########################################################################
    !########################################################################
    subroutine TLab_AddTranspose(a, nra, nca, ma, b, mb, locBlock)
        integer(wi), intent(in) :: nra      ! Number of rows in a
        integer(wi), intent(in) :: nca      ! Number of columns in b
        integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
        integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
        real(wp), intent(in) :: a(ma, nca)    ! Input array
        real(wp), intent(inout) :: b(mb, nra) ! Transposed array
        type(block_dt), optional :: locBlock

        ! -------------------------------------------------------------------
        integer(wi) jb, kb
        integer(wi) k, j, jj, kk
        integer(wi) last_k, last_j

        ! ###################################################################
        if (present(locBlock)) then
            jb = locBlock%jb
            kb = locBlock%kb
        else
            jb = trans_x%jb
            kb = trans_x%kb
        end if

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

        do k = 1, last_k - 1
            do j = last_j, nra
                b(k, j) = b(k, j) + a(j, k)
            end do
        end do

        return
    end subroutine TLab_AddTranspose

    !########################################################################
    !########################################################################
    subroutine TLab_SubtractTranspose(a, nra, nca, ma, b, mb, locBlock)
        integer(wi), intent(in) :: nra      ! Number of rows in a
        integer(wi), intent(in) :: nca      ! Number of columns in b
        integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
        integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
        real(wp), intent(in) :: a(ma, *)    ! Input array
        real(wp), intent(inout) :: b(mb, *) ! Transposed array
        type(block_dt), optional :: locBlock

        ! -------------------------------------------------------------------
        integer(wi) jb, kb
        integer(wi) k, j, jj, kk
        integer(wi) last_k, last_j

        ! ###################################################################
        if (present(locBlock)) then
            jb = locBlock%jb
            kb = locBlock%kb
        else
            jb = trans_x%jb
            kb = trans_x%kb
        end if

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

        do k = 1, last_k - 1
            do j = last_j, nra
                b(k, j) = b(k, j) - a(j, k)
            end do
        end do

        return
    end subroutine TLab_SubtractTranspose

    !########################################################################
    !########################################################################
    subroutine TLab_Transpose_Complex(a, nra, nca, ma, b, mb, locBlock)
        integer(wi), intent(in) :: nra      ! Number of rows in a
        integer(wi), intent(in) :: nca      ! Number of columns in b
        integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
        integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
        complex(wp), intent(in) :: a(ma, *) ! Input array
        complex(wp), intent(out) :: b(mb, *) ! Transposed array
        type(block_dt), optional :: locBlock

        ! -------------------------------------------------------------------
        integer(wi) jb, kb
        integer(wi) k, j, jj, kk
        integer(wi) last_k, last_j

        ! ###################################################################
        if (present(locBlock)) then
            jb = locBlock%jb/2
            kb = locBlock%kb
        else
            jb = trans_x%jb/2
            kb = trans_x%kb
        end if

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

        do k = 1, last_k - 1
            do j = last_j, nra
                b(k, j) = a(j, k)
            end do
        end do

        return
    end subroutine TLab_Transpose_Complex

    ! !########################################################################
    ! !########################################################################
    ! subroutine TLab_Transpose_INT1(a, nra, nca, ma, b, mb)
    !     integer(wi), intent(in) :: nra      ! Number of rows in a
    !     integer(wi), intent(in) :: nca      ! Number of columns in b
    !     integer(wi), intent(in) :: ma       ! Leading dimension on the input matrix a
    !     integer(wi), intent(in) :: mb       ! Leading dimension on the output matrix b
    !     integer(1), intent(in) :: a(ma, *) ! Input array
    !     integer(1), intent(out) :: b(mb, *) ! Transposed array

    !     ! -------------------------------------------------------------------
    !     integer(wi) jb, kb

    !     integer(wi) k, j, jj, kk
    !     integer(wi) last_k, last_j

    !     ! ###################################################################
    !     jb = trans_x%jb*4
    !     kb = trans_x%kb*2

    !     kk = 1; jj = 1

    !     do k = 1, nca - kb + 1, kb;
    !         do j = 1, nra - jb + 1, jb;
    !             do jj = j, j + jb - 1
    !                 do kk = k, k + kb - 1
    !                     b(kk, jj) = a(jj, kk)
    !                 end do
    !             end do
    !         end do
    !     end do

    !     last_k = kk
    !     last_j = jj

    !     do k = last_k, nca
    !         do j = 1, nra
    !             b(k, j) = a(j, k)
    !         end do
    !     end do

    !     do k = 1, last_k - 1
    !         do j = last_j, nra
    !             b(k, j) = a(j, k)
    !         end do
    !     end do

    !     return
    ! end subroutine TLab_Transpose_INT1

end module TLab_Transpose
