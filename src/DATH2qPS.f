      SUBROUTINE DATH2qPS
*     -------------------
*
      IMPLICIT NONE
      INTEGER LIN,LOUT,NX,IX,ISET,ITE
      REAL*8 XXV,XV,A(501),B(501),C(501),D(501),H2qPS
      CHARACTER*80 STAR
      COMMON/H2qPS/ H2qPS(22,501,4)
      COMMON/XV/    XV(501)
      COMMON/XXV/   XXV(501)
*
*
*    >>> READ IN THE GRID DATA AND COMPUTE THE SPLINE PARAMETERS
*
      LIN = 10
      LOUT= 11
      NX  = 501
*
      OPEN(UNIT=LIN,FILE='BLUM/xvalues.dat',STATUS='OLD')
      DO IX = 1,NX
          READ(LIN,*) XV(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS20000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
          DO IX = 1,NX
          READ(LIN,*) H2qPS(1,IX,1)
          A(IX)=H2qPS(1,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(1,IX,2)=B(IX)
          H2qPS(1,IX,3)=C(IX)
          H2qPS(1,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS20100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(2,IX,1)
          A(IX)=H2qPS(2,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(2,IX,2)=B(IX)
          H2qPS(2,IX,3)=C(IX)
          H2qPS(2,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS20200.dat',STATUS='OLD')
      READ(LIN,1000) STAR
          DO IX = 1,NX
          READ(LIN,*) H2qPS(3,IX,1)
          A(IX)=H2qPS(3,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(3,IX,2)=B(IX)
          H2qPS(3,IX,3)=C(IX)
          H2qPS(3,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS21000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(4,IX,1)
          A(IX)=H2qPS(4,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(4,IX,2)=B(IX)
          H2qPS(4,IX,3)=C(IX)
          H2qPS(4,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS22000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(5,IX,1)
          A(IX)=H2qPS(5,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(5,IX,2)=B(IX)
          H2qPS(5,IX,3)=C(IX)
          H2qPS(5,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(6,IX,1)
          A(IX)=H2qPS(6,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(6,IX,2)=B(IX)
          H2qPS(6,IX,3)=C(IX)
          H2qPS(6,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30001.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(7,IX,1)
          A(IX)=H2qPS(7,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(7,IX,2)=B(IX)
          H2qPS(7,IX,3)=C(IX)
          H2qPS(7,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(8,IX,1)
          A(IX)=H2qPS(8,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
          DO IX = 1,NX
          H2qPS(8,IX,2)=B(IX)
          H2qPS(8,IX,3)=C(IX)
          H2qPS(8,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(9,IX,1)
          A(IX)=H2qPS(9,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(9,IX,2)=B(IX)
          H2qPS(9,IX,3)=C(IX)
          H2qPS(9,IX,4)=D(IX)
          END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30110.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(10,IX,1)
          A(IX)=H2qPS(10,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(10,IX,2)=B(IX)
          H2qPS(10,IX,3)=C(IX)
          H2qPS(10,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30200.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(11,IX,1)
          A(IX)=H2qPS(11,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(11,IX,2)=B(IX)
          H2qPS(11,IX,3)=C(IX)
          H2qPS(11,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30210.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(12,IX,1)
          A(IX)=H2qPS(12,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(12,IX,2)=B(IX)
          H2qPS(12,IX,3)=C(IX)
          H2qPS(12,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30300.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(13,IX,1)
          A(IX)=H2qPS(13,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(13,IX,2)=B(IX)
          H2qPS(13,IX,3)=C(IX)
          H2qPS(13,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS30310.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(14,IX,1)
          A(IX)=H2qPS(14,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(14,IX,2)=B(IX)
          H2qPS(14,IX,3)=C(IX)
          H2qPS(14,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS31000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(15,IX,1)
          A(IX)=H2qPS(15,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(15,IX,2)=B(IX)
          H2qPS(15,IX,3)=C(IX)
          H2qPS(15,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS31010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(16,IX,1)
          A(IX)=H2qPS(16,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(16,IX,2)=B(IX)
          H2qPS(16,IX,3)=C(IX)
          H2qPS(16,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS31100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(17,IX,1)
          A(IX)=H2qPS(17,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(17,IX,2)=B(IX)
          H2qPS(17,IX,3)=C(IX)
          H2qPS(17,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS32000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(18,IX,1)
          A(IX)=H2qPS(18,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(18,IX,2)=B(IX)
          H2qPS(18,IX,3)=C(IX)
          H2qPS(18,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS32010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(19,IX,1)
          A(IX)=H2qPS(19,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(19,IX,2)=B(IX)
          H2qPS(19,IX,3)=C(IX)
          H2qPS(19,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS32100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(20,IX,1)
          A(IX)=H2qPS(20,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(20,IX,2)=B(IX)
          H2qPS(20,IX,3)=C(IX)
          H2qPS(20,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS33000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(21,IX,1)
          A(IX)=H2qPS(21,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(21,IX,2)=B(IX)
          H2qPS(21,IX,3)=C(IX)
          H2qPS(21,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2qPS33010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2qPS(22,IX,1)
          A(IX)=H2qPS(22,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2qPS(22,IX,2)=B(IX)
          H2qPS(22,IX,3)=C(IX)
          H2qPS(22,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
*     >> WRITE TO DISK
*
      OPEN(UNIT=LOUT,FILE='BLUM/XXV.dat',STATUS='REPLACE')
      DO IX=1,501
      WRITE(LOUT,*) XV(IX)
      END DO
*
      OPEN(UNIT=LOUT,FILE='BLUM/H2qPS.dat',STATUS='REPLACE')
      DO ISET=1,22
      DO IX=1,501
      DO ITE=1,4
      WRITE(LOUT,*) H2qPS(ISET,IX,ITE)
      END DO
      END DO
      END DO
*
      RETURN
1000  FORMAT(A80)
      END
