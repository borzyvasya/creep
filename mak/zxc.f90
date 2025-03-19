program main
    implicit none
    double precision :: start_p, pi, root, &
         tolerance, f_lower_bound, f_upper_bound, f_root, h
    integer :: N, iteration_count

    start_p = 0.0d0                    
    pi = atan(1.0d0) * 4.0d0      
    tolerance = 1.0e-4  ! Epsilon            
    N = 1000                  
    iteration_count = 0                    

    f_lower_bound = f(start_p)
    f_upper_bound = f(pi)
    
    call loop(start_p, pi, tolerance, N, root, f_root, h, iteration_count, &
                       f_lower_bound, f_upper_bound)
    
    call adjust_function_value(f_root, tolerance)
    
    write(*, "(A, F14.10)") "Root with precision(10): ", root
    write(*, "(A, F14.10)") "Function value at root: ", f_root
    write(*, "(A, I10)")    "Steps taken: ", iteration_count
    write(*, "(A, F14.10)") "Step length (h): ", h

    call check_convergence(iteration_count, N, h, tolerance)

contains
    ! Функция f(x) = sin(x) - cos(x)
    double precision function f(x)
        implicit none
        double precision, intent(in) :: x
        f = sin(x) - cos(x)
    end function f

    subroutine loop(a, b, tol, max_iter, r, fr, step, iter, fa, fb)
        implicit none
        double precision, intent(inout) :: a, b, fa, fb
        double precision, intent(in) :: tol
        integer, intent(in) :: max_iter
        double precision, intent(out) :: r, fr, step
        integer, intent(out) :: iter
        
        iter = 0
        do while (iter < max_iter)
            r = (a + b) / 2.0d0    ! Середина интервала
            fr = f(r)
            step = b - a           ! Длина интервала

            ! Проверка условия точности
            if (step <= tol) then
                exit  
            end if

            ! Обновление границ
            if (fa * fr < 0.0d0) then
                b = r
                fb = fr
            else
                a = r
                fa = fr
            end if
            iter = iter + 1
        end do
    end subroutine loop

    subroutine adjust_function_value(fr, tol)
        implicit none
        double precision, intent(inout) :: fr
        double precision, intent(in) :: tol
        
        if (abs(fr) < tol) then
            fr = 0.0d0
        end if
    end subroutine adjust_function_value

    ! Для проверки сходимости
    subroutine check_convergence(iter, max_iter, step, tol)
        implicit none
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

end program main