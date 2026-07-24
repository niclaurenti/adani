      SUBROUTINE DATH2gS
*     ------------------
*
      IMPLICIT NONE
      INTEGER LIN,LOUT,NX,IX,ISET,ITE
      REAL*8 XXV,XV,A(501),B(501),C(501),D(501),H2gS
      CHARACTER*80 STAR
      COMMON/H2gS/ H2gS(29,501,4)
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
      OPEN(UNIT=LIN,FILE='BLUM/H2gS10000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
          DO IX = 1,NX
          READ(LIN,*) H2gS(1,IX,1)
          A(IX)=H2gS(1,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(1,IX,2)=B(IX)
          H2gS(1,IX,3)=C(IX)
          H2gS(1,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS10100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(2,IX,1)
          A(IX)=H2gS(2,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(2,IX,2)=B(IX)
          H2gS(2,IX,3)=C(IX)
          H2gS(2,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS11000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
          DO IX = 1,NX
          READ(LIN,*) H2gS(3,IX,1)
          A(IX)=H2gS(3,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(3,IX,2)=B(IX)
          H2gS(3,IX,3)=C(IX)
          H2gS(3,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS20000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(4,IX,1)
          A(IX)=H2gS(4,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(4,IX,2)=B(IX)
          H2gS(4,IX,3)=C(IX)
          H2gS(4,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS20100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(5,IX,1)
          A(IX)=H2gS(5,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(5,IX,2)=B(IX)
          H2gS(5,IX,3)=C(IX)
          H2gS(5,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS20200.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(6,IX,1)
          A(IX)=H2gS(6,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(6,IX,2)=B(IX)
          H2gS(6,IX,3)=C(IX)
          H2gS(6,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS21000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(7,IX,1)
          A(IX)=H2gS(7,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(7,IX,2)=B(IX)
          H2gS(7,IX,3)=C(IX)
          H2gS(7,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS21100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(8,IX,1)
          A(IX)=H2gS(8,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
          DO IX = 1,NX
          H2gS(8,IX,2)=B(IX)
          H2gS(8,IX,3)=C(IX)
          H2gS(8,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS22000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(9,IX,1)
          A(IX)=H2gS(9,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(9,IX,2)=B(IX)
          H2gS(9,IX,3)=C(IX)
          H2gS(9,IX,4)=D(IX)
          END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(10,IX,1)
          A(IX)=H2gS(10,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(10,IX,2)=B(IX)
          H2gS(10,IX,3)=C(IX)
          H2gS(10,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30001.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(11,IX,1)
          A(IX)=H2gS(11,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(11,IX,2)=B(IX)
          H2gS(11,IX,3)=C(IX)
          H2gS(11,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(12,IX,1)
          A(IX)=H2gS(12,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(12,IX,2)=B(IX)
          H2gS(12,IX,3)=C(IX)
          H2gS(12,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(13,IX,1)
          A(IX)=H2gS(13,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(13,IX,2)=B(IX)
          H2gS(13,IX,3)=C(IX)
          H2gS(13,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30110.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(14,IX,1)
          A(IX)=H2gS(14,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(14,IX,2)=B(IX)
          H2gS(14,IX,3)=C(IX)
          H2gS(14,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30200.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(15,IX,1)
          A(IX)=H2gS(15,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(15,IX,2)=B(IX)
          H2gS(15,IX,3)=C(IX)
          H2gS(15,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30210.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(16,IX,1)
          A(IX)=H2gS(16,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(16,IX,2)=B(IX)
          H2gS(16,IX,3)=C(IX)
          H2gS(16,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30300.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(17,IX,1)
          A(IX)=H2gS(17,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(17,IX,2)=B(IX)
          H2gS(17,IX,3)=C(IX)
          H2gS(17,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS30310.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(18,IX,1)
          A(IX)=H2gS(18,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(18,IX,2)=B(IX)
          H2gS(18,IX,3)=C(IX)
          H2gS(18,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS31000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(19,IX,1)
          A(IX)=H2gS(19,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(19,IX,2)=B(IX)
          H2gS(19,IX,3)=C(IX)
          H2gS(19,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*

      OPEN(UNIT=LIN,FILE='BLUM/H2gS31010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(20,IX,1)
          A(IX)=H2gS(20,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(20,IX,2)=B(IX)
          H2gS(20,IX,3)=C(IX)
          H2gS(20,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS31100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(21,IX,1)
          A(IX)=H2gS(21,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(21,IX,2)=B(IX)
          H2gS(21,IX,3)=C(IX)
          H2gS(21,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS31110.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(22,IX,1)
          A(IX)=H2gS(22,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(22,IX,2)=B(IX)
          H2gS(22,IX,3)=C(IX)
          H2gS(22,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS31200.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(23,IX,1)
          A(IX)=H2gS(23,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(23,IX,2)=B(IX)
          H2gS(23,IX,3)=C(IX)
          H2gS(23,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS31210.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(24,IX,1)
          A(IX)=H2gS(24,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(24,IX,2)=B(IX)
          H2gS(24,IX,3)=C(IX)
          H2gS(24,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS32000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(25,IX,1)
          A(IX)=H2gS(25,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(25,IX,2)=B(IX)
          H2gS(25,IX,3)=C(IX)
          H2gS(25,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS32010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(26,IX,1)
          A(IX)=H2gS(26,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(26,IX,2)=B(IX)
          H2gS(26,IX,3)=C(IX)
          H2gS(26,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS32100.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(27,IX,1)
          A(IX)=H2gS(27,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(27,IX,2)=B(IX)
          H2gS(27,IX,3)=C(IX)
          H2gS(27,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS33000.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(28,IX,1)
          A(IX)=H2gS(28,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(28,IX,2)=B(IX)
          H2gS(28,IX,3)=C(IX)
          H2gS(28,IX,4)=D(IX)
      END DO
      CLOSE(LIN)
*
      OPEN(UNIT=LIN,FILE='BLUM/H2gS33010.dat',STATUS='OLD')
      READ(LIN,1000) STAR
      DO IX = 1,NX
          READ(LIN,*) H2gS(29,IX,1)
          A(IX)=H2gS(29,IX,1)
      END DO
      CALL SPLINE(NX,XV,A,B,C,D)
      DO IX = 1,NX
          H2gS(29,IX,2)=B(IX)
          H2gS(29,IX,3)=C(IX)
          H2gS(29,IX,4)=D(IX)
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
      OPEN(UNIT=LOUT,FILE='BLUM/H2gS.dat',STATUS='REPLACE')
      DO ISET=1,29
      DO IX=1,501
      DO ITE=1,4
      WRITE(LOUT,*) H2gS(ISET,IX,ITE)
      END DO
      END DO
      END DO
*
      RETURN
1000  FORMAT(A80)
      END
