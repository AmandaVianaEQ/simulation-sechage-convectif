MODULE donnees
    IMPLICIT NONE

    ! Paramètres géométriques et numériques
    DOUBLE PRECISION :: e, deltaz, deltat, t_fin, t_stockage
    INTEGER :: NN

    ! Paramètres physiques
    DOUBLE PRECISION :: Wini, Whyg, psec, g
    DOUBLE PRECISION :: mvsec, r, Vinf, Tinf, HR, Thumide
    DOUBLE PRECISION :: B1, GAMMA1
    DOUBLE PRECISION :: Weq

CONTAINS

    SUBROUTINE lecture_donnees()
        IMPLICIT NONE

        OPEN(10, FILE='donnees.txt', STATUS='OLD')
            READ(10,*) e
            READ(10,*) Wini
            READ(10,*) Whyg
            READ(10,*) psec
            READ(10,*) g
            READ(10,*) mvsec
            READ(10,*) r
            READ(10,*) Vinf
            READ(10,*) Tinf
            READ(10,*) HR
            READ(10,*) Thumide
            READ(10,*) NN
            READ(10,*) B1
            READ(10,*) GAMMA1
            READ(10,*) deltat
            READ(10,*) t_fin
            READ(10,*) t_stockage
        CLOSE(10)

        deltaz = e / (NN - 1)

    END SUBROUTINE

END MODULE
