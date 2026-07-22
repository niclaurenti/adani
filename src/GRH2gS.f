      REAL*8 FUNCTION GRH2gS(IO,X,AS,LM,LQ,NF,FLAV)
*     ---------------------------------------------
*
      IMPLICIT NONE
      INTEGER IO,I1,ISEARCHX,NX,K
      REAL*8 X,AS,LM,LQ,NF,FLAV,G,DX
      REAL*8 H2gS
      REAL*8 XXV
      REAL*8 XV
      REAL*8 T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13
      REAL*8 T14,T15,T16,T17,T18,T19,T20,T21,T22,T23,T24
      REAL*8 T25,T26,T27,T28,T29
      COMMON/H2gS/  H2gS(29,501,4)
      COMMON/XXV/   XXV(501)
      COMMON/XV/    XV(501)
*
      NX=501
      IF(IO.LT.1) WRITE(6,*) ' >>> IO=', IO, 'STOP: GRH2gS'
      IF(IO.GT.3) WRITE(6,*) ' >>> IO=', IO, 'STOP: GRH2gS'
      IF(IO.EQ.1) GOTO 101
      IF(IO.EQ.2) GOTO 102
      IF(IO.EQ.3) GOTO 103
      STOP
101   CONTINUE
      I1=ISEARCHX(NX,XXV,X)
      DX=X-XXV(I1)
      WRITE(6,*) 'DX=',DX
*
      T1  =0.0D0
      T2  =0.0D0
      T3  =0.0D0
      DO K=1,4
      T1  = T1   + H2gS(1,I1,K)*DX**(K-1)
      T2  = T2   + H2gS(2,I1,K)*DX**(K-1)
      T3  = T3   + H2gS(3,I1,K)*DX**(K-1)
      END DO
*
      G = AS*(T1+T2*LQ+T3*LM)
      GOTO 200
102   CONTINUE
      I1=ISEARCHX(NX,XXV,X)
      DX=X-XXV(I1)
      WRITE(6,*) 'DX=',DX
*
      T4  =0.0D0
      T5  =0.0D0
      T6  =0.0D0
      T7  =0.0D0
      T8  =0.0D0
      T9  =0.0D0
      DO K=1,4
      T4  = T4   + H2gS(4,I1,K)*DX**(K-1)
*      WRITE(6,*) 'K,T4=',K,T4
      T5  = T5   + H2gS(5,I1,K)*DX**(K-1)
      T6  = T6   + H2gS(6,I1,K)*DX**(K-1)
      T7  = T7   + H2gS(7,I1,K)*DX**(K-1)
      T8  = T8   + H2gS(8,I1,K)*DX**(K-1)
      T9  = T9   + H2gS(9,I1,K)*DX**(K-1)
      END DO
*
      G = AS**2*(T4+T5*LQ+T6*LQ**2+T7*LM+T8*LM*LQ+T9*LM**2)
      GOTO 200
103   CONTINUE
      I1=ISEARCHX(NX,XXV,X)
      DX=X-XXV(I1)
      WRITE(6,*) 'DX=',DX
*
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
      T23 =0.0D0
      T24 =0.0D0
      T25 =0.0D0
      T26 =0.0D0
      T27 =0.0D0
      T28 =0.0D0
      T29 =0.0D0
      WRITE(6,*) 'I1=',I1
      DO K=1,4
      T10 = T10  + H2gS(10,I1,K)*DX**(K-1)
      T11 = T11  + H2gS(11,I1,K)*DX**(K-1)
      T12 = T12  + H2gS(12,I1,K)*DX**(K-1)
      T13 = T13  + H2gS(13,I1,K)*DX**(K-1)
      T14 = T14  + H2gS(14,I1,K)*DX**(K-1)
      T15 = T15  + H2gS(15,I1,K)*DX**(K-1)
      T16 = T16  + H2gS(16,I1,K)*DX**(K-1)
      T17 = T17  + H2gS(17,I1,K)*DX**(K-1)
      T18 = T18  + H2gS(18,I1,K)*DX**(K-1)
      T19 = T19  + H2gS(19,I1,K)*DX**(K-1)
      T20 = T20  + H2gS(20,I1,K)*DX**(K-1)
      T21 = T21  + H2gS(21,I1,K)*DX**(K-1)
      T22 = T22  + H2gS(22,I1,K)*DX**(K-1)
      T23 = T23  + H2gS(23,I1,K)*DX**(K-1)
      T24 = T24  + H2gS(24,I1,K)*DX**(K-1)
      T25 = T25  + H2gS(25,I1,K)*DX**(K-1)
      T26 = T26  + H2gS(26,I1,K)*DX**(K-1)
      T27 = T27  + H2gS(27,I1,K)*DX**(K-1)
      T28 = T28  + H2gS(28,I1,K)*DX**(K-1)
      T29 = T29  + H2gS(29,I1,K)*DX**(K-1)
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
      G = as**3*(T10 + FLAV*T11 + NF*T12 + LQ*T13 + LQ*NF*T14
     &   + LQ**2*T15 + LQ**2*NF*T16 + LQ**3*T17 + LQ**3*NF*T18
     &   + LM*T19 + LM*NF*T20 + LM*LQ*T21 + LM*LQ*NF*T22
     &   + LM*LQ**2*T23 + LM*LQ**2*NF*T24 + LM**2*T25
     &   + LM**2*NF*T26+ LM**2*LQ*T27 + LM**3*T28
     &   + LM**3*NF*T29)
*      WRITE(6,*) 'G=',G
*
200   CONTINUE
*
      GRH2gS= G
*
      RETURN
      END
