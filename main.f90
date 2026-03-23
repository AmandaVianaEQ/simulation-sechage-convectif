PROGRAM projetCSII
    USE donnees
    USE physique
    USE elements_finis

    IMPLICIT NONE
    INTEGER :: i, j, choix, choix_th
    DOUBLE PRECISION :: t, Fm, FmNN, t_next_stockage

    DOUBLE PRECISION, ALLOCATABLE :: WI(:), WII(:), W_newI(:), W_newII(:)
    DOUBLE PRECISION, ALLOCATABLE :: Winter(:), W_moy(:), k1(:), k2(:), A_moy(:)
    DOUBLE PRECISION, ALLOCATABLE :: A(:), D(:), B(:)
    DOUBLE PRECISION, ALLOCATABLE :: M(:,:), K(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: A_global(:,:), C_global(:)

    CALL lecture_donnees()

    WRITE(*,*) "Calculer Thumide automatiquement?"
    WRITE(*,*) "(1)Oui   (2)Non"
    READ(*,*) choix_th
    IF (choix_th == 1) THEN
        CALL Calc_Thumide()
    END IF

    WRITE(*,'(A,F8.3,A)') " Thumide utilisee: ", Thumide, " C"

    CALL Calc_Weq()

    ALLOCATE(WI(NN), WII(NN), W_newI(NN), W_newII(NN))
    ALLOCATE(Winter(NN), W_moy(NN), k1(NN), k2(NN), A_moy(NN))
    ALLOCATE(A(NN), D(NN), B(NN))
    ALLOCATE(M(NN,NN), K(NN,NN))
    ALLOCATE(A_global(NN,NN), C_global(NN))

    WI = Wini
    WII = Wini

    WRITE(*,*) "Quelle methode souhaitez-vous utiliser?"
    WRITE(*,*) "(1)Methode I (2)Methode II (3)Methodes I et II"
    READ(*,*)choix

    t = 0.0

    IF (choix==1.OR.choix==3) OPEN(20, FILE='results.txt',STATUS='REPLACE')
    IF (choix==2.OR.choix==3) OPEN(30, FILE='resultsII.txt',STATUS='REPLACE')
    IF (choix==3) OPEN(40, FILE='results_comparasion.txt', STATUS='REPLACE')

    t=0
    t_next_stockage = t_stockage

    IF (choix==1.OR.choix==3) WRITE(20,'(f10.2,500(1x,f10.6))') 0.0d0, (WI(j),j=1,NN)
    IF (choix==2.OR.choix==3) WRITE(30,'(f10.2,500(1x,f10.6))') 0.0d0, (WII(j),j=1,NN)
    IF (choix==3) WRITE(40,'(f10.2,500(1x,f10.6))') 0.0d0, (0.0d0,j=1,NN)

    DO WHILE (t < t_fin)

! -------------------------- METHODE I --------------------------

        IF (choix==1.OR.choix==3) THEN
            CALL calcul_DA(WI,D,A)

            CALL Calcul_Fm(Fm,WI(1))
            CALL Calcul_Fm(FmNN,WI(NN))

            CALL calcul_K(K, D)
            CALL calcul_M(M)
            CALL calcul_B(B, Fm, FmNN, A)

            A_global = M + K
            C_global = B + MATMUL(M, WI)

            CALL resolution_systeme(A_global, C_global, W_newI)
            WI = W_newI
        END IF

! -------------------------- METHODE II --------------------------

        IF (choix==2.OR.choix==3) THEN
            CALL calcul_DA(WII, D, A)

            CALL Calcul_Fm(Fm,WII(1))
            CALL Calcul_Fm(FmNN,WII(NN))

            !RUNGE-KUTTA 2
            DO i=1, NN-1
                k1(i)= -(deltat/deltaz)*g*(A(i+1) - A(i))
            END DO

            k1(NN) = -(deltat/deltaz)*g*(A(NN) - A(NN-1))

            W_moy = WII + k1

            CALL calcul_DA(W_moy, D, A_moy)

            DO i = 1, NN-1
                k2(i) = -(deltat/deltaz)*g*(A_moy(i+1) - A_moy(i))
            END DO

            k2(NN) = -(deltat/deltaz)*g*(A_moy(NN) - A_moy(NN-1))

            Winter = WII + 0.5d0*(k1+k2)

            IF (Vinf > 5.0d0) THEN
                Winter(1)  = Weq
                Winter(NN) = Weq
            END IF

            CALL calcul_DA(WII, D, A)
            CALL calcul_K(K, D)
            CALL calcul_M(M)
            CALL calcul_BII(B, Fm, FmNN)

            A_global = M + K
            C_global = B + MATMUL(M, Winter)

            CALL resolution_systeme(A_global, C_global, W_newII)
            WII = W_newII

        END IF

        t = t + deltat

         IF (t >= t_next_stockage) THEN
            IF (choix==1.OR.choix==3) WRITE(20,'(f10.0,500(1x,f10.6))') t, (WI(j),j=1,NN)
            IF (choix==2.OR.choix==3) WRITE(30,'(f10.0,500(1x,f10.6))') t, (WII(j),j=1,NN)
            IF (choix==3) WRITE(40,'(f10.0,500(1x,f10.6))') t, (WI(j)-WII(j),j=1,NN)
            t_next_stockage = t_next_stockage + t_stockage
         END IF
    END DO

    IF (choix==1.OR.choix==3) CLOSE(20)
    IF (choix==2.OR.choix==3) CLOSE(30)
    IF (choix==3) CLOSE(40)

    DEALLOCATE(WI, WII, W_newI, W_newII, Winter, W_moy, k1, k2, A_moy,A, D, B, M, K, A_global, C_global)

END PROGRAM



