      REAL*8 FUNCTION REGH2gS(IO,X,AS,LM,LQ,NF,FLAV)
*     ----------------------------------------------
*
      IMPLICIT NONE
      INTEGER IO
      REAL*8 X,AS,LM,LQ,NF,FLAV
      REAL*8 GRH2gS,SXH2gS,LXH2gS,T1,T2,T3
*
      T1=SXH2gS(IO,X,AS,LM,LQ,NF,FLAV)
      T2=LXH2gS(IO,X,AS,LM,LQ,NF,FLAV)
      T3=GRH2gS(IO,X,AS,LM,LQ,NF,FLAV)
*
      REGH2gS = T1+T2+T3
*
      RETURN
      END
