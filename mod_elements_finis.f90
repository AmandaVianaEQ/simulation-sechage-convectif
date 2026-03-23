MODULE elements_finis
    USE donnees
    IMPLICIT NONE

CONTAINS

! -------------------- Calcul de la matrice K --------------------
    SUBROUTINE calcul_K(K,D)
        IMPLICIT NONE
        INTEGER :: i
        DOUBLE PRECISION, DIMENSION(NN)    :: D
        DOUBLE PRECISION, DIMENSION(NN,NN) :: K

        K=0

        K(1,1) = (deltat/(2*deltaz))*(D(1)+D(2))
        K(1,2) = -(deltat/(2*deltaz))*(D(1)+D(2))

        DO i=2,NN-1
            K(i,i-1) = -(deltat/(2*deltaz))*(D(i-1)+D(i))
            K(i,i) = (deltat/(2*deltaz))*(D(i-1)+D(i)) + (deltat/(2*deltaz))*(D(i+1)+D(i))
            K(i,i+1) = -(deltat/(2*deltaz))*(D(i+1)+D(i))
        END DO

        K(NN,NN-1) = -(deltat/(2*(deltaz)))*(D(NN-1)+D(NN))
        K(NN,NN) = (deltat/(2*deltaz))*(D(NN-1)+D(NN))

    END SUBROUTINE calcul_K

! -------------------- Calcul de la matrice M --------------------
    SUBROUTINE calcul_M(M)
        IMPLICIT NONE
        INTEGER :: i
        DOUBLE PRECISION, DIMENSION(NN,NN) :: M

        M=0

        M(1,1) = (1.0/3.0)*deltaz
        M(1,2) = (1.0/6.0)*deltaz

        DO i=2,NN-1
            M(i,i-1)=(1.0/6.0)*deltaz
            M(i,i)=(1.0/3.0)*deltaz + (1.0/3.0)*deltaz
            M(i,i+1)=(1.0/6.0)*deltaz
        END DO

        M(NN,NN-1)=(1.0/6.0)*deltaz
        M(NN,NN)=(1.0/3.0)*deltaz

    END SUBROUTINE calcul_M

! -------------------- Calcul de la matrice B - Methode I --------------------
    SUBROUTINE calcul_B(B,Fm,FmNN,A)
        IMPLICIT NONE
        INTEGER :: i
        DOUBLE PRECISION :: Fm, FmNN
        DOUBLE PRECISION, DIMENSION(NN) :: A
        DOUBLE PRECISION, DIMENSION(NN) :: B

        B(1) = -deltat*((g*(A(1)+A(2)))+((Fm/psec)))

        DO i=2,NN-1
            B(i)= deltat*g*((A(i-1)-A(i+1))/2)
        END DO

        B(NN) = deltat*(((g*(A(NN-1)+A(NN)))/2)-((FmNN/psec)))

    END SUBROUTINE calcul_B

! -------------------- Calcul de la matrice B - Methode II --------------------
    SUBROUTINE calcul_BII(B,Fm,FmNN)
        IMPLICIT NONE
        INTEGER :: i
        DOUBLE PRECISION :: Fm, FmNN
        DOUBLE PRECISION, DIMENSION(NN) :: B

        B = 0
        B(1) = -deltat*((Fm/psec))
        B(NN) = deltat*(-((FmNN/psec)))

    END SUBROUTINE calcul_BII

! -------------------- Resolution du Systeme - Algorithme de Thomas --------------------
    SUBROUTINE resolution_systeme(A_global,C_global,W_new)
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(NN,NN) :: A_global
        DOUBLE PRECISION, DIMENSION(NN) :: C_global
        DOUBLE PRECISION, DIMENSION(NN) :: W_new
        DOUBLE PRECISION, DIMENSION(NN) :: c_prime, d_prime
        DOUBLE PRECISION :: m, GV
        INTEGER :: i

        GV = 1.0d40

        IF(Vinf>5.0d0) THEN
            A_global(1,1) = GV
            A_global(NN,NN) = GV
            C_global(1) = GV*Weq
            C_global(NN) = GV*Weq

        END IF

        c_prime(1) = A_global(1,2) / A_global(1,1)
        d_prime(1) = C_global(1) / A_global(1,1)

        DO i = 2, NN - 1
            m = A_global(i,i) - A_global(i,i-1)*c_prime(i-1)
            c_prime(i) = A_global(i,i+1)/m
            d_prime(i) = (C_global(i) - A_global(i,i-1)*d_prime(i-1))/m
        END DO

        m = A_global(NN,NN) - A_global(NN,NN-1)*c_prime(NN-1)
        d_prime(NN) = (C_global(NN) - A_global(NN,NN-1) * d_prime(NN-1))/m

        W_new(NN) = d_prime(NN)

        DO i = NN-1,1,-1
            W_new(i) = d_prime(i) - c_prime(i) * W_new(i+1)
        END DO

    END SUBROUTINE

END MODULE
