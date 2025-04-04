program IterativeMethods
    use iters_m
    implicit none

    integer, parameter :: N = 10
    double precision, parameter :: EPS = 1.0d-6, OMEGA = 1.3d0
    double precision :: A(N,N), b(N), x(N)

    call generate_matrix(A, b, N)
    write(*, *) "Matrix A and b:"
    call print_matrix(A, b, N)

    write(*, *) "--------------------------------------------------"
    write(*, *) "Simple Iterations (Jacobi Method):"
    call simple_iteration(A, b, x, EPS, N)
    call output_solutions(A, b, x, "Simple Iterations", N)

    write(*, *) "--------------------------------------------------"
    write(*, *) "Gauss-Seidel Method:"
    call gauss_seidel(A, b, x, EPS, N)
    call output_solutions(A, b, x, "Gauss-Seidel", N)

    write(*, *) "--------------------------------------------------"
    write(*, *) "Successive Over-Relaxation (SOR) Method:"
    call sor(A, b, x, EPS, OMEGA, N)
    call output_solutions(A, b, x, "SOR", N)
    write(*, *) "--------------------------------------------------"

end program IterativeMethods
    
    subroutine generate_matrix(A, b, N)
        integer, intent(in) :: N
        double precision, intent(out) :: A(N,N), b(N)
        integer :: i, j
        double precision :: sum, CF, rand_val

        CF = 1.0d1
        call random_seed()

        do i = 1, N
            sum = 0.0d0
            do j = 1, N
                if (i /= j) then
                    call random_number(A(i,j))
                    A(i,j) = A(i,j) * CF
                    sum = sum + abs(A(i,j))
                end if
            end do
            call random_number(rand_val)
            A(i,i) = sum + rand_val * CF + 1.0d0
        end do

        do i = 1, N
            call random_number(b(i))
            b(i) = b(i) * CF
        end do
    end subroutine

    subroutine print_matrix(A, b, N)
        integer, intent(in) :: N
        double precision, intent(in) :: A(N,N), b(N)
        integer :: i, j
        do i = 1, N
            write(*, '(10F9.3, " | ", F10.3)') (A(i,j), j=1,N), b(i)
        end do
    end subroutine

    subroutine simple_iteration(A, b, x, epsilon, N)
        integer, intent(in) :: N
        double precision, intent(in) :: A(N,N), b(N), epsilon
        double precision, intent(out) :: x(N)
        double precision :: x_new(N), error
        integer :: i, j, iter

        x = 0.0d0
        error = epsilon + 1
        iter = 0

        do while (error > epsilon)
            error = 0.0d0
            do i = 1, N
                x_new(i) = b(i)
                do j = 1, N
                    if (i /= j) x_new(i) = x_new(i) - A(i,j) * x(j)
                end do
                x_new(i) = x_new(i) / A(i,i)
            end do

            do i = 1, N
                error = max(error, abs(x_new(i) - x(i)))
                x(i) = x_new(i)
            end do

            iter = iter + 1
            write(*, '(A, I4, A, F12.6)') "Iteration ", iter, ", error: ", error
        end do
    end subroutine

    subroutine gauss_seidel(A, b, x, epsilon, N)
        integer, intent(in) :: N
        double precision, intent(in) :: A(N,N), b(N), epsilon
        double precision, intent(out) :: x(N)
        double precision :: error, x_old
        integer :: i, j, iter

        x = 0.0d0
        error = epsilon + 1
        iter = 0

        do while (error > epsilon)
            error = 0.0d0

            do i = 1, N
                x_old = x(i)
                x(i) = b(i)
                do j = 1, N
                    if (i /= j) x(i) = x(i) - A(i,j) * x(j)
                end do
                x(i) = x(i) / A(i,i)
                error = max(error, abs(x(i) - x_old))
            end do

            iter = iter + 1
            write(*, '(A, I4, A, F12.6)') "Iteration ", iter, ", error: ", error
        end do
    end subroutine

    subroutine sor(A, b, x, epsilon, omega, N)
        integer, intent(in) :: N
        double precision, intent(in) :: A(N,N), b(N), epsilon, omega
        double precision, intent(out) :: x(N)
        double precision :: error, x_old, gs_update
        integer :: i, j, iter

        x = 0.0d0
        error = epsilon + 1
        iter = 0

        do while (error > epsilon)
            error = 0.0d0

            do i = 1, N
                x_old = x(i)
                gs_update = b(i)
                do j = 1, N
                    if (i /= j) gs_update = gs_update - A(i,j) * x(j)
                end do
                gs_update = gs_update / A(i,i)
                x(i) = (1 - omega) * x_old + omega * gs_update
                error = max(error, abs(x(i) - x_old))
            end do

            iter = iter + 1
            write(*, '(A, I4, A, F12.6)') "Iteration ", iter, ", error: ", error
        end do
    end subroutine

    subroutine output_solutions(A, b, x, method_name, N)
        integer, intent(in) :: N
        double precision, intent(in) :: A(N,N), b(N), x(N)
        character(len=*), intent(in) :: method_name
        double precision :: Ax(N), residual
        integer :: i, j

        do i = 1, N
            Ax(i) = 0.0d0
            do j = 1, N
                Ax(i) = Ax(i) + A(i,j) * x(j)
            end do
        end do

        residual = 0.0d0
        do i = 1, N
            residual = max(residual, abs(Ax(i) - b(i)))
        end do

        write(*, *) "Solution (", method_name, "):"
        do i = 1, N
            write(*, '(A, I4, A, F12.6)') "x[", i, "] = ", x(i)
        end do
        write(*, '(A, F12.6)') "Residual: ", residual
    end subroutine

