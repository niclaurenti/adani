      REAL*8 FUNCTION GRH2qPS(IO,X,AS,LM,LQ,NF,FLAV)
*     -----------------------------------------------
*
      IMPLICIT NONE
      INTEGER IO,I1,ISEARCHX,NX,K
      REAL*8 X,AS,LM,LQ,NF,FLAV,G,DX
      REAL*8 H2qPS
      REAL*8 XXV
      REAL*8 XV
      REAL*8 T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13
      REAL*8 T14,T15,T16,T17,T18,T19,T20,T21,T22
      COMMON/H2qPS/  H2qPS(22,501,4)
      COMMON/XXV/   XXV(501)
      COMMON/XV/    XV(501)
*
      NX=501
      IF(IO.LT.2) WRITE(6,*) ' >>> IO=', IO, 'STOP: GRH2qPS'
      IF(IO.GT.3) WRITE(6,*) ' >>> IO=', IO, 'STOP: GRH2qPS'
      IF(IO.EQ.2) GOTO 101
      IF(IO.EQ.3) GOTO 102
      STOP
101   CONTINUE
      I1=ISEARCHX(NX,XXV,X)
      DX=X-XXV(I1)
      WRITE(6,*) 'DX=',DX
*
      T1  =0.0D0
      T2  =0.0D0
      T3  =0.0D0
      T4  =0.0D0
      T5  =0.0D0
      DO K=1,4
      T1  = T1   + H2qPS(1,I1,K)*DX**(K-1)
      T2  = T2   + H2qPS(2,I1,K)*DX**(K-1)
      T3  = T3   + H2qPS(3,I1,K)*DX**(K-1)
      T4  = T4   + H2qPS(4,I1,K)*DX**(K-1)
      T5  = T5   + H2qPS(5,I1,K)*DX**(K-1)
      END DO
*
      G = AS**2*(T1+T2*LQ+T3*LQ**2+T4*LM+T5*LM**2)
      GOTO 200
102   CONTINUE
      I1=ISEARCHX(NX,XXV,X)
      DX=X-XXV(I1)
      WRITE(6,*) 'DX=',DX
*
      T6  =0.0D0
      T7  =0.0D0
      T8  =0.0D0
      T9  =0.0D0
      T10 =0.0D0
      T11 =0.0D0
      T12 =0.0D0
      T13 =0.0D0
      T14 =0.0D0
      T15 =0.0D0
      T16 =0.0D0
      T17 =0.0D0
      T18 =0.0D0
      T19 =0.0D0
      T20 =0.0D0
      T21 =0.0D0
      T22 =0.0D0
      WRITE(6,*) 'I1=',I1
      DO K=1,4
      T6  = T6   + H2qPS(6,I1,K)*DX**(K-1)
      T7  = T7   + H2qPS(7,I1,K)*DX**(K-1)
      T8  = T8   + H2qPS(8,I1,K)*DX**(K-1)
      T9  = T9   + H2qPS(9,I1,K)*DX**(K-1)
      T10 = T10  + H2qPS(10,I1,K)*DX**(K-1)
      T11 = T11  + H2qPS(11,I1,K)*DX**(K-1)
      T12 = T12  + H2qPS(12,I1,K)*DX**(K-1)
      T13 = T13  + H2qPS(13,I1,K)*DX**(K-1)
      T14 = T14  + H2qPS(14,I1,K)*DX**(K-1)
      T15 = T15  + H2qPS(15,I1,K)*DX**(K-1)
      T16 = T16  + H2qPS(16,I1,K)*DX**(K-1)
      T17 = T17  + H2qPS(17,I1,K)*DX**(K-1)
      T18 = T18  + H2qPS(18,I1,K)*DX**(K-1)
      T19 = T19  + H2qPS(19,I1,K)*DX**(K-1)
      T20 = T20  + H2qPS(20,I1,K)*DX**(K-1)
      T21 = T21  + H2qPS(21,I1,K)*DX**(K-1)
      T22 = T22  + H2qPS(22,I1,K)*DX**(K-1)
*      WRITE(6,*) 'K,T21=',K,T21
      END DO
*      WRITE(6,*) 'T1 =',T1
*      WRITE(6,*) 'T2 =',T2
*      WRITE(6,*) 'T3 =',T3
*      WRITE(6,*) 'T4 =',T4
*      WRITE(6,*) 'T5 =',T5
*      WRITE(6,*) 'T6 =',T6
*      WRITE(6,*) 'T7 =',T7
*      WRITE(6,*) 'T8 =',T8
*      WRITE(6,*) 'T9 =',T9
*      WRITE(6,*) 'T10=',T10
*      WRITE(6,*) 'T11=',T11
*      WRITE(6,*) 'T12=',T12
*      WRITE(6,*) 'T13=',T13
*
      G = as**3*(T6 + FLAV*T7 + NF*T8 + LQ*T9 + LQ*NF*T10 + LQ**2*T11
     &   + LQ**2*NF*T12 + LQ**3*T13 + LQ**3*NF*T14 + LM*T15 + LM*NF*T16
     &   + LM*LQ*T17 + LM**2*T18 + LM**2*NF*T19 + LM**2*LQ*T20
     &   + LM**3*T21 + LM**3*NF*T22)
*      WRITE(6,*) 'G=',G
*
200   CONTINUE
*
      GRH2qPS= G
*
      RETURN
      END
