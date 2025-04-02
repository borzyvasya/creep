module task
    implicit none
    
    interface
        subroutine output(A, B, N, NU)
            integer, intent(in) :: N, NU
            real, intent(in) :: A(N, N), B(N)
        end subroutine output
        
        subroutine input(A, B, N, filename, file_unit)
            integer, intent(in) :: N, file_unit
            real, intent(out) :: A(N, N), B(N)
            character(len=20), intent(in) :: filename
        end subroutine input
    end interface

end module task
