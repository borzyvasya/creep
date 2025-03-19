module zxc_m
    interface
        double precision function f(x)
            implicit none
            double precision, intent(in) :: x
        end function f
        
        subroutine loop(a, b, tol, max_iter, r, fr, step, iter, fa, fb)
            implicit none
            double precision, intent(inout) :: a, b, fa, fb
            double precision, intent(in) :: tol
            integer, intent(in) :: max_iter
            double precision, intent(out) :: r, fr, step
            integer, intent(out) :: iter
        end subroutine loop
        
        subroutine adjust_function_value(fr, tol)
            implicit none
            double precision, intent(inout) :: fr
            double precision, intent(in) :: tol
        end subroutine adjust_function_value
        
        subroutine check_convergence(iter, max_iter, step, tol)
            implicit none
            integer, intent(in) :: iter, max_iter
            double precision, intent(in) :: step, tol
        end subroutine check_convergence
    end interface
end module