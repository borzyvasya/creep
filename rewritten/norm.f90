program norms
    implicit none
    integer, parameter :: SIZE = 4
    real :: x(SIZE) = [12.7, 4.4, 2.3, 5.5]
    
    ! Вывод результатов с помощью write
    write(*,'(A, F8.1)') "Norm of the array (L1): ", normL1(x, SIZE)
    write(*,'(A, F8.1)') "Norm of the array (L2): ", normL2(x, SIZE)
    write(*,'(A, F6.1)') "Norm of the array (Linf): ", l_inf_norm(x, SIZE)

contains
    ! Функция для L1 нормы (сумма абсолютных значений)
    real function normL1(arr, n)
        integer, intent(in) :: n
        real, intent(in) :: arr(n)  ! Размер массива передается как аргумент
        real :: sum
        integer :: i
    
        sum = 0.0
        do i = 1, n
            sum = sum + abs(arr(i))
        end do
        normL1 = sum
    end function normL1

    ! Функция для L2 нормы 
    real function normL2(arr, n)
        integer, intent(in) :: n
        real, intent(in) :: arr(n)  ! Размер массива передается как аргумент
        real :: sum
        integer :: i
    
        sum = 0.0
        do i = 1, n
            sum = sum + sqrt(arr(i) * arr(i))  
        end do
        normL2 = sum
    end function normL2 

    ! Функция для L∞ нормы (максимальное абсолютное значение)
    real function l_inf_norm(arr, n)
        integer, intent(in) :: n
        real, intent(in) :: arr(n)  ! Размер массива передается как аргумент
        real :: max_value, abs_value
        integer :: i
    
        max_value = abs(arr(1))
        do i = 1, n
            abs_value = abs(arr(i))
            if (abs_value > max_value) max_value = abs_value
        end do
        l_inf_norm = max_value
    end function l_inf_norm

end program norms