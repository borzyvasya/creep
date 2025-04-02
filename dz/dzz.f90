program main
    use dzz_m
    implicit none
    double precision :: a, b, root, tolerance, h
    integer :: N, iteration_count

    a = -1.0d0
    b = 1.0d0
    tolerance = 1.0e-4
    N = 1000
    iteration_count = 0

    call loop(a, b, tolerance, N, root, h, iteration_count)
    
    write(*, "(A, F14.10)") "Root with precision(10): ", root
    write(*, "(A, I10)")    "Steps taken: ", iteration_count
    write(*, "(A, F14.10)") "Step length (h): ", h

    call check_convergence(iteration_count, N, h, tolerance)
    
end program main
     subroutine check_convergence(iter, max_iter, step, tol)
        integer, intent(in) :: iter, max_iter
        double precision, intent(in) :: step, tol
        
        if (iter >= max_iter) then
            write(*,*) "The method did not converge within the maximum number of iterations."
        else if (step > tol) then
            write(*,*) "The final interval is larger than tolerance."
        else
            write(*,*) "Success check!"
        end if
    end subroutine check_convergence
    
    function f(x, a, b) result(res)
        double precision, intent(in) :: x, a, b
        double precision :: res

        res = atan(1+2*x)-cos(x)
    end function f

    subroutine loop(a, b, tol, max_iter, r, step, iter)
    double precision, intent(inout) :: a, b
    double precision, intent(in) :: tol
    integer, intent(in) :: max_iter
    double precision, intent(out) :: r, step
    integer, intent(out) :: iter
    double precision :: fr

    iter = 0
    do while (iter < max_iter)
        r = (a + b) / 2.0d0    
        step = b - a           
        fr = f(r, a, b)     

        if (step <= tol) then
            exit           
        end if

        if (f(a, a, b) * fr < 0.0d0) then
            b = r              
        else
            a = r              
        end if
        iter = iter + 1
    end do
end subroutine loop
   