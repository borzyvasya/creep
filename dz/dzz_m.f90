module dzz_m
    implicit none
    
interface
    function f(x, a, b) result(res)
        double precision, intent(in) :: x, a, b
        double precision :: res
    end function f

    subroutine check_convergence(iter, max_iter, step, tol)
        integer, intent(in) :: iter, max_iter
        double precision, intent(in) :: step, tol
    end subroutine check_convergence

    subroutine loop(a, b, tol, max_iter, r, step, iter)
        double precision, intent(inout) :: a, b
        double precision, intent(in) :: tol
        integer, intent(in) :: max_iter
        double precision, intent(out) :: r, step
        integer, intent(out) :: iter
    end subroutine loop
    
end interface
end module dzz_m