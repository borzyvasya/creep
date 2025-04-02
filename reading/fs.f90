    program fs
    use task
    implicit none
    integer, parameter :: ifl = 10, N = 10, MULT = 10
    real :: A(N, N), B(N), X(N), Y
    integer :: I, J, NU
    character(len=20) :: fname = "fin.txt"  
    integer :: file_unit

    file_unit = ifl  

    call input(A, B, N, fname, file_unit)

    call output(A, B, N, NU)

    close(file_unit)  


    end program fs

    subroutine output(A, B, N, NU)
        implicit none
        integer, intent(in) :: N, NU
        real, intent(in) :: A(N, N), B(N)
        integer :: I, J
        character(len=10) :: NR
        character(len=20) :: FMT
        
        write(*, '(/,A)') "Matrix A:"
        do I = 1, N
            write(*, '(10(F6.2x))') (A(I,J), J=1,N)
        end do

        write(*, '(A)') "Vector B:"
        write(*, '(10(F6.2,2X))') (B(I), I=1,N)

        
    end subroutine
    
        subroutine input(A, B, N, filename, file_unit)
        implicit none
        integer, intent(in) :: N, file_unit
        real, intent(out) :: A(N, N), B(N)
        character(len=20), intent(in) :: filename
        integer :: I, J

        open(unit=file_unit, file=filename)

        do I = 1, N
            read(file_unit, *) (A(I, J), J = 1, N)
        end do

        read(file_unit, *) (B(I), I = 1, N)

        close(file_unit)
    end subroutine