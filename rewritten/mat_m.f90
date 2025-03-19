module mat_mod
    

    interface
        subroutine MATSUMM(M1, M2, M3, SZ)
            implicit none
            integer, intent(in) :: SZ
            real, intent(in):: M1(SZ, SZ), M2(SZ, SZ) 
            real, intent(out) :: M3(SZ, SZ)
        end subroutine MATSUMM

        subroutine OUTPUT(AR, MSZ, SZ)
            implicit none
            integer, intent(in) :: MSZ, SZ
            real, intent(in) :: AR(MSZ, SZ)
        end subroutine OUTPUT
    end interface

end module mat_mod
