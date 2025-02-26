program gauss_solver
    implicit none
    integer, parameter :: Q = 10, Q1 = Q + 1, QE = Q * 2
    real :: a(Q, Q1), x(Q), cn

    call rand_mat(a)
    write(*,*) "Original matrix with b:"
    call print_mat(a)
    call gauss_solve(a, x)
    write(*,*) "Checking solution:"
    call check_sol(a, x)
    cn = mat_norm(a) * inv_norm(a)
    write(*,*) "Condition number: ", cn

contains
    subroutine rand_mat(m)
        real, intent(out) :: m(Q, Q1)
        call random_seed()
        call random_number(m)
        m = m / 1.0e4
    end subroutine

    subroutine print_mat(m)
        real, intent(in) :: m(Q, Q1)
        integer :: i
        do i = 1, Q
            write(*, "(10F8.2, A, F8.2)") m(i, 1:Q), " | ", m(i, Q1)
        end do
        write(*,*)
    end subroutine

    subroutine gauss_solve(m, s)
        real, intent(inout) :: m(Q, Q1)
        real, intent(out) :: s(Q)
        integer :: i, k, j, mr
        real :: me, f, t
        
        do i = 1, Q
            me = abs(m(i, i)); mr = i
            do k = i + 1, Q
                if (abs(m(k, i)) > me) me = abs(m(k, i)); mr = k
            end do
            if (mr /= i) then
                do j = 1, Q1
                    t = m(mr, j); m(mr, j) = m(i, j); m(i, j) = t
                end do
            end if
            if (abs(m(i, i)) < 1.0e-6) then
                write(*,*) "Matrix is singular!"; return
            end if
            do k = 1, Q
                if (k /= i) then
                    f = m(k, i) / m(i, i)
                    m(k, i:Q1) = m(k, i:Q1) - f * m(i, i:Q1)
                    where (abs(m(k, :)) < 1.0e-4) m(k, :) = 0.0
                end if
            end do
        end do
        do i = 1, Q
            if (abs(m(i, i)) < 1.0e-7) then
                write(*,*) "Cannot solve!"; return
            end if
            s(i) = m(i, Q1) / m(i, i)
        end do
    end subroutine

    real function mat_norm(m)
        real, intent(in) :: m(Q, Q1)
        mat_norm = maxval(sum(abs(m(:, 1:Q)), dim=1))
    end function

    real function inv_norm(m)
        real, intent(in) :: m(Q, Q1)
        real :: aug(Q, QE)
        integer :: i, k, j
        real :: f, p
        
        aug(:, 1:Q) = m(:, 1:Q); aug(:, Q+1:QE) = 0.0
        forall(i=1:Q) aug(i, Q+i) = 1.0
        do i = 1, Q
            if (abs(aug(i, i)) < 1.0e-9) then
                inv_norm = -1.0; return
            end if
            
            do k = i + 1, Q
                f = aug(k, i) / aug(i, i)
                aug(k, :) = aug(k, :) - f * aug(i, :)
            end do
            
        end do
        do i = Q, 1, -1
            p = aug(i, i)
            if (abs(p) < 1.0e-4) then
                inv_norm = -1.0; return
            end if
            aug(i, :) = aug(i, :) / p
            do k = 1, i - 1
                aug(k, :) = aug(k, :) - aug(k, i) * aug(i, :)
            end do
        end do
        inv_norm = maxval(sum(abs(aug(:, Q+1:QE)), dim=1))
    end function

    subroutine check_sol(m, s)
        real, intent(in) :: m(Q, Q1), s(Q)
        integer :: i
        real :: sum
        do i = 1, Q
            sum = dot_product(m(i, 1:Q), s)
            if (abs(sum - m(i, Q1)) > 1.0e-4) then
                write(*, "(A, I2, A, F10.4, A, F10.4)") "Equation ", i, " fails: ", sum, " != ", m(i, Q1)
                write(*,*) "Solution is incorrect"
                return
            end if
        end do
        write(*,*) "Solution is correct"
    end subroutine

end program gauss_solver