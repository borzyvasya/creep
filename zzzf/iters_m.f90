module iters_m
    implicit none
    
    interface
        subroutine generate_matrix(A, b, N)
            integer, intent(in) :: N
            double precision, intent(out) :: A(N,N), b(N)
        end subroutine

        subroutine print_matrix(A, b, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N)
        end subroutine

        subroutine simple_iteration(A, b, x, epsilon, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), epsilon
            double precision, intent(out) :: x(N)
        end subroutine

        subroutine gauss_seidel(A, b, x, epsilon, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), epsilon
            double precision, intent(out) :: x(N)
        end subroutine

        subroutine sor(A, b, x, epsilon, omega, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), epsilon, omega
            double precision, intent(out) :: x(N)
        end subroutine

        subroutine output_solutions(A, b, x, method_name, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), x(N)
            character(len=*), intent(in) :: method_name
        end subroutine
    end interface

end module iters_m
