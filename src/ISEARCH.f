* ---------------------------------------------------------------------
      INTEGER FUNCTION ISEARCHX(NX,XV,Y)
* ---------------------------------------------------------------------
* THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
* X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
*
* FUNCTION TAKEN FROM AAC GROUP (KUMANO et al.)
*
      IMPLICIT NONE
      REAL*8   XV
      REAL*8   Y
      INTEGER  NX,MIN,MAX,MID
*
      DIMENSION XV(NX)
*
      MIN=1
      MAX=NX+1
*
   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.XV(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10
*
      ISEARCHX=MIN
*
      RETURN
      END
