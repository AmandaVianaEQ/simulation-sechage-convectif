MODULE physique
    USE donnees
    IMPLICIT NONE

CONTAINS

! -------------------- Calcul de la teneur en eau d'équilibre --------------------
    SUBROUTINE Calc_Weq()
        IMPLICIT NONE
        DOUBLE PRECISION A_aw,B_aw,C_aw

        A_aw = -2.86d-5*(Tinf+273)**2 - 1.07d-2*(Tinf+273) + 10.24
        B_aw = -5.41d-4*(Tinf+273) + 1.01
        C_aw = 4.97d-6*(Tinf+273)**2 - 2.67d-3*(Tinf+273) + 0.35

        Weq = (DLOG((C_aw-DLOG(HR))/A_aw))/(100*DLOG(B_aw))

    END SUBROUTINE Calc_Weq

! -------------------- Calcul Temperature Humide --------------------

     SUBROUTINE Calc_Thumide()
        DOUBLE PRECISION :: Th_a, Th_b, Th_moy
        DOUBLE PRECISION :: f_a, f_moy
        DOUBLE PRECISION :: PsatTh, PsatTinf
        DOUBLE PRECISION :: DeltaHv_a, DeltaHv_moy
        DOUBLE PRECISION :: Hv, Cpv, Cpl
        DOUBLE PRECISION :: terme1
        INTEGER :: i
        INTEGER, PARAMETER :: imax = 1000
        DOUBLE PRECISION, PARAMETER :: tol = 1.0d-6

        Hv = 2501000.0d0
        Cpv = 2000.0d0
        Cpl = 4185.0d0

        !hm/ht * Mv/R
        terme1 = (1.0d0/1000.0d0) * (mvsec/r)

        PsatTinf = 100 * EXP(75.051-(7293.3/(Tinf+273))-8.593*DLOG(Tinf+273)+0.00617*(Tinf+273))

        Th_a = 0.0d0
        Th_b = Tinf

        !METHODE DICHOTOMIE
        DO i = 1, imax
            Th_moy = (Th_a + Th_b) / 2.0d0

            DeltaHv_a = Hv - (Cpl-Cpv)*Th_a
            PsatTh = 100 * EXP(75.051-(7293.3/(Th_a+273))-8.593*DLOG(Th_a+273)+0.00617*(Th_a+273))
            f_a = (Tinf - Th_a) - terme1 * ((PsatTh/(Th_a+273))- (HR*PsatTinf/(Tinf+273)))*DeltaHv_a

            DeltaHv_moy = Hv - (Cpl-Cpv)*Th_moy
            PsatTh = 100 * EXP(75.051-(7293.3/(Th_moy+273))-8.593*DLOG(Th_moy+273)+0.00617*(Th_moy+273))
            f_moy = (Tinf - Th_moy) - terme1 * ((PsatTh/(Th_moy+273))- (HR*PsatTinf/(Tinf+273)))*DeltaHv_moy

            IF (ABS(Th_b - Th_a) < tol) EXIT

            IF (f_a * f_moy < 0.0d0) THEN
                Th_b = Th_moy
            ELSE
                Th_a = Th_moy
            END IF
        END DO

        Thumide = Th_moy

    END SUBROUTINE Calc_Thumide

! -------------------- Calcul Flux de Masse --------------------

    SUBROUTINE Calcul_Fm(Fm, W)
        IMPLICIT NONE
        DOUBLE PRECISION Fm, W
        DOUBLE PRECISION T, pv, pvinf, PressionVsat, PressionVapeur, PVsatTinf
        DOUBLE PRECISION A_aw, B_aw, C_aw, a_temp, b_temp, hm, aw

    ! CALCUL hm
        hm = 9.454d-3 * Vinf**0.5003

    ! CALCUL TEMPERATURE
        IF(W>Whyg) then
            T = Thumide
        ELSE
            a_temp = (Tinf - (Thumide)) / ((Tinf)*(Thumide)*(Whyg-Weq))
            b_temp = 1.0d0/(Thumide) - a_temp*Whyg
            T = 1.0d0/(a_temp*W + b_temp)
        END IF

    ! CALCUL aw
        A_aw = -2.86d-5*(T+273)**2 - 1.07d-2*(T+273) + 10.24
        B_aw = -5.41d-4*(T+273) + 1.01
        C_aw = 4.97d-6*(T+273)**2 - 2.67d-3*(T+273) + 0.35

        aw = EXP(-A_aw*B_aw**(100*W)+C_aw)

    ! CALCUL PRESSION
        PressionVsat = 100 * EXP(75.051-(7293.3/(T+273))-8.593*DLOG(T+273)+0.00617*(T+273))
        PVsatTinf = 100 * EXP(75.051-(7293.3/(Tinf+273))-8.593*DLOG(Tinf+273)+0.00617*(Tinf+273))

        IF(W>Whyg) THEN
            PressionVapeur = PressionVsat
        ELSE
            PressionVapeur = aw*PressionVsat
        END IF

    ! CALCUL FM
        pv = (mvsec*PressionVapeur)/(R*(T+273))
        pvinf = (mvsec*PVsatTinf*HR)/(R*(Tinf+273))

        Fm = hm*(pv-pvinf)

    END SUBROUTINE Calcul_Fm

! -------------------- Calcul DA --------------------
    SUBROUTINE calcul_DA(W,D,A)
        DOUBLE PRECISION, DIMENSION(NN) :: W
        DOUBLE PRECISION, DIMENSION(NN) :: D, A
        INTEGER :: i

        DO i = 1,NN
            ! CALCUL D
            IF(W(i)<Whyg) then
                D(i) = B1*5d-07
            ELSEIF(W(i)>=Whyg) then
                IF(W(i)<(Whyg+0.15)) then
                    D(i) = B1*(5d-07 + 5d-06*(W(i)-Whyg))
                ELSEIF(W(i)>=(Whyg+0.15)) then
                    D(i) =B1*(2d-06 - 1d-06*(W(i)-Whyg-0.15))
                END IF
            END IF

            ! CALCUL A
                A(i) = GAMMA1*EXP(-20.31+1.40*W(i)**2.0 * DLOG(W(i))-(0.87/W(i)**2.0))
        END DO

    END SUBROUTINE calcul_DA

END MODULE
