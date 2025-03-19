program main
    use mat_mod
    implicit none
    integer :: I, J, N, M
    parameter (N = 4, M = 5)
    real :: A(N, N), B(N, N), C(N, N), D(M, M)
    
    do I = 1, N
        do J = 1, N
            A(I, J) = 1.0
            B(I, J) = 2.0
        end do
    end do

    call MATSUMM(A, B, C, N)
    call OUTPUT(C, N, N)
    
    do I = 1, M
        do J = 1, N
            D(I, J) = I + J
        end do
    end do

    call OUTPUT(D, M, N)
    end program main

    subroutine MATSUMM(M1, M2, M3, SZ)
    implicit none
    integer I, J, SZ
    real M1(SZ, SZ), M2(SZ, SZ), M3(SZ, SZ)
    
    do I = 1, SZ
        do J = 1, SZ
            M3(I, J) = M1(I, J) + M2(I, J)
        end do
    end do
    return 
end

subroutine OUTPUT(AR, MSZ, SZ)
    implicit none
    integer I, J, SZ, MSZ
    real AR(MSZ, SZ)
    
    write(*, 60) "Matrix size: ", MSZ, "X", SZ
    do I = 1, MSZ
        write(*,50) (AR(I, J), J = 1, SZ)
    end do
    
50  format(5f5.1)
60  format(/A13, I1, A1, I1)
    return
end