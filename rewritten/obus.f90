program obuslov
    implicit none
    integer, parameter :: Q = 10, Q1 = Q + 1, QE = Q * 2
    real :: a(Q, Q1), x(Q), cn

    ! Генерация случайной матрицы
    call rand_mat(a)
    
    write(*,*) "Original matrix:"
    call print_mat(a)

    ! Решение системы и вывод верхнетреугольной матрицы
    call forwardSub(a, x)
    
    write(*,*) "Upper triangular matrix:"
    call print_mat(a)
    
    write(*,*) "Checking solution:"
    call check_sol(a, x)
    
    ! Вычисление числа обусловленности 
    cn = mat_norm(a) * inv_norm(a)
    write(*, '(A, F8.2)') "Condition number: ", cn


contains
    
    ! Генерация случайной матрицы
    subroutine rand_mat(m)
        real, intent(out) :: m(Q, Q1)
        call random_seed()
        call random_number(m)
        m = m * 15
    end subroutine rand_mat

    ! Печать матрицы
    subroutine print_mat(m)
        real, intent(in) :: m(Q, Q1)
        integer :: i
        do i = 1, Q
            write(*, "(10F8.2, A, F8.2)") m(i, 1:Q), " | ", m(i, Q1)
        end do
        write(*,*)
    end subroutine print_mat

    ! Решение системы методом Гаусса без обратного хода
    subroutine forwardSub(m, s)
        real, intent(inout) :: m(Q, Q1)
        real, intent(out) :: s(Q)
        integer :: i, k, j, mr
        real :: me, f, t
    
        ! Прямой ход 
        do i = 1, Q
            me = abs(m(i, i)); mr = i
            do k = i + 1, Q
                if (abs(m(k, i)) > me) then
                    me = abs(m(k, i)); mr = k
                end if
            end do
            if (mr /= i) then
                do j = 1, Q1
                    t = m(mr, j); m(mr, j) = m(i, j); m(i, j) = t
                end do
            end if
            if (abs(m(i, i)) < 1.0e-6) then
                write(*,*) "Matrix is singular!"; return
            end if
            do k = i + 1, Q  ! Обнуляем только ниже диагонали
                f = m(k, i) / m(i, i)
                do j = i, Q1  ! Начинаем с диагонального элемента
                    m(k, j) = m(k, j) - f * m(i, j)
                end do
                where (abs(m(k, :)) < 1.0e-6) m(k, :) = 0.0
            end do
        end do
    
        ! Прямое вычисление решения из верхнетреугольной матрицы
        do i = Q, 1, -1  ! Теперь идем снизу вверх
            if (abs(m(i, i)) < 1.0e-7) then
                write(*,*) "Cannot solve!"; return
            end if
            s(i) = m(i, Q1)
            do j = i + 1, Q
                s(i) = s(i) - m(i, j) * s(j)
            end do
            s(i) = s(i) / m(i, i)
        end do
    end subroutine forwardSub

    ! Вычисление нормы ||A||
    real function mat_norm(m)
        real, intent(in) :: m(Q, Q1)
        mat_norm = maxval(sum(abs(m(:, 1:Q)), dim=1))
    end function mat_norm

    ! Вычисление нормы ||A⁻¹||
    real function inv_norm(m)
        real, intent(in) :: m(Q, Q1)
        real :: aug(Q, QE)  ! Расширенная матрица [A | E]
        integer :: i, k, j
        real :: f, p
    
        aug(:, 1:Q) = m(:, 1:Q)
        aug(:, Q+1:QE) = 0.0
        forall(i=1:Q) aug(i, Q+i) = 1.0

        do i = 1, Q
            p = abs(aug(i, i))
            k = i
            do j = i + 1, Q
                if (abs(aug(j, i)) > p) then
                    p = abs(aug(j, i))
                    k = j
                end if
            end do
        
            if (p < 1.0e-9) then
                inv_norm = -1.0  ! Матрица сингулярная
                return
            end if
        
            if (k /= i) then
                aug([i, k], :) = aug([k, i], :)
            end if

            aug(i, :) = aug(i, :) / aug(i, i)
            do k = 1, Q
                if (k /= i) then
                    f = aug(k, i)
                    aug(k, :) = aug(k, :) - f * aug(i, :)
                end if
            end do
        end do
        inv_norm = maxval(sum(abs(aug(:, Q+1:QE)), dim=1))
    end function inv_norm
    
    ! Проверка решения
    subroutine check_sol(m, s)
        real, intent(in) :: m(Q, Q1), s(Q)
        integer :: i, j
        real :: sum, sol
        ! Проверка для верхнетреугольной матрицы
        do i = 1, Q
            sum = 0.0
            do j = i, Q  ! Учитываем только верхний треугольник
                sum = sum + m(i, j) * s(j)
            end do
            sol = m(i, Q1)
            if (abs(sum - sol) > 1.0e-4) then
                write(*, "(A, I2, A, F10.4, A, F10.4)") & 
                    "Equation ", i, " fails: ", sum, " != ", sol
                write(*,*) "Solution is incorrect"
                return
            end if
        end do
        write(*,*) "Solution is correct"
    end subroutine check_sol

end program obuslov