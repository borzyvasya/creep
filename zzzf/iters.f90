program IterativeMethods
    use IterativeMethodsModule
    implicit none

    integer, parameter :: N = 3
    double precision, parameter :: EPS = 0.0002, OMEGA = 1.15
    double precision :: A(N,N), b(N), x(N)

    call generate_matrix(A, b, N)
    write(*, *) "Matrix A and b:"
    call print_matrix(A, b, N)

    write(*, *) "--------------------------------------------------"
    write(*, *) "Simple Iterations (Jacobi Method):"
    call simple_iteration(A, b, x, EPS, N)
    call output_solutions(A, b, x, "Simple Iterations", N)

    write(*, *) "--------------------------------------------------"
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
    
    A = reshape([ & 
        3.00, 1.00, -1.00, & 
        2.00, 6.00, 3.00, & 
        -1.00, 1.00, 4.00], [N, N])

    b = [7.00, -2.00, 4.00]
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
    double precision :: x_new(N), Ax(N), error
    integer :: i, iter

    x = 0.0d0
    error = epsilon + 1.0d0
    iter = 0

    do while (error > epsilon .and. iter < 1000)
        iter = iter + 1

        Ax = matmul(A, x)

        do i = 1, N
            if (abs(A(i,i)) < 1.0d-10) then
                write(*,*) "Error: Zero or near-zero diagonal element at ", i
                return
            end if
            x_new(i) = x(i) + (b(i) - Ax(i)) / A(i,i)
        end do

        error = maxval(abs(x_new - x))
        x = x_new

        write(*, '(A, I4, A, F12.6)') "Iteration ", iter, ", error: ", error
    end do

end subroutine

subroutine gauss_seidel(A, b, x, epsilon, N)
    integer, intent(in) :: N
    double precision, intent(in) :: A(N,N), b(N), epsilon
    double precision, intent(out) :: x(N)
    double precision :: error
    integer :: i, j, iter, max_iter
    double precision :: Ax(N)  

    x = 0.0d0
    error = epsilon + 1.0d0
    iter = 0
    max_iter = 1000

    do while (error > epsilon .and. iter < max_iter)
        iter = iter + 1
        error = 0.0d0

        do i = 1, N
            x(i) = (b(i) - sum(A(i, 1:i-1) * x(1:i-1)) - &
             sum(A(i, i+1:N) * x(i+1:N))) / A(i, i)
        end do

        
        Ax = matmul(A, x)
        error = maxval(abs(Ax - b))

        write(*, '(A, I4, A, F12.6)') "Iteration ", iter, ", error: ", error
    end do

    if (iter == max_iter) then
        write(*, *) "Warning: Reached maximum iterations in Gauss-Seidel method"
    end if
end subroutine

subroutine sor(A, b, x, epsilon, omega, N)
    integer, intent(in) :: N
    double precision, intent(in) :: A(N,N), b(N), epsilon, omega
    double precision, intent(out) :: x(N)
    double precision :: error, x_old, gs_update
    integer :: i, j, iter, max_iter
    double precision :: Ax(N)  ! Declare Ax to be a vector of length N

    x = 0.0d0
    error = epsilon + 1.0d0
    iter = 0
    max_iter = 1000

    do while (error > epsilon .and. iter < max_iter)
        iter = iter + 1
        error = 0.0d0

        do i = 1, N
            x_old = x(i)
            gs_update = b(i)
            do j = 1, N
                if (i /= j) gs_update = gs_update - A(i,j) * x(j)
            end do
            if (abs(A(i,i)) < 1.0d-10) then
                write(*,*) "Error: Zero or near-zero diagonal element at ", i
                return
            end if
            gs_update = gs_update / A(i,i)
            x(i) = (1 - omega) * x_old + omega * gs_update
            error = max(error, abs(x(i) - x_old))
        end do

        write(*, '(A, I4, A, F12.6)') "Iteration ", iter, ", error: ", error
    end do
end subroutine


subroutine output_solutions(A, b, x, method_name, N)
    integer, intent(in) :: N
    double precision, intent(in) :: A(N,N), b(N), x(N)
    character(len=*), intent(in) :: method_name
    double precision :: Ax(N), residual
    integer :: i

    Ax = matmul(A, x)
    
    residual = maxval(abs(Ax - b))

    write(*, *) "Solution (", method_name, "):"
    do i = 1, N
        write(*, '(A, I2, A, F12.6)') "x[", i, "] = ", x(i)
    end do
    write(*, '(A, F12.6)') "Residual: ", residual
end subroutine
