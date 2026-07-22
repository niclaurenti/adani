      REAL*8 FUNCTION REGH2qPS(IO,X,AS,LM,LQ,NF,FLAV)
*     -----------------------------------------------
*
      IMPLICIT NONE
      INTEGER IO
      REAL*8 X,AS,LM,LQ,NF,FLAV
      REAL*8 GRH2qPS,SXH2qPS,T1,T3
*
      T1=SXH2qPS(IO,X,AS,LM,LQ,NF,FLAV)
      T3=GRH2qPS(IO,X,AS,LM,LQ,NF,FLAV)
*
      REGH2qPS = T1+T3
*
      RETURN
      END
