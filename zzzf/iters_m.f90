module IterativeMethodsModule
    implicit none

    interface
        subroutine generate_matrix(A, b, N)
            integer, intent(in) :: N
            double precision, intent(out) :: A(N,N), b(N)
        end subroutine generate_matrix

        subroutine print_matrix(A, b, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N)
        end subroutine print_matrix

        subroutine diagonal(A, diag)
            double precision, intent(in) :: A(:,:)
            double precision, intent(out) :: diag(:)
        end subroutine diagonal

        subroutine simple_iteration(A, b, x, epsilon, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N, N), b(N), epsilon
            double precision, intent(out) :: x(N)
        end subroutine simple_iteration

        subroutine gauss_seidel(A, b, x, epsilon, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), epsilon
            double precision, intent(out) :: x(N)
        end subroutine gauss_seidel

        subroutine sor(A, b, x, epsilon, omega, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), epsilon, omega
            double precision, intent(out) :: x(N)
        end subroutine sor

        subroutine output_solutions(A, b, x, method_name, N)
            integer, intent(in) :: N
            double precision, intent(in) :: A(N,N), b(N), x(N)
            character(len=*), intent(in) :: method_name
        end subroutine output_solutions
    end interface

end module IterativeMethodsModule
