      SUBROUTINE SELFD(I,J,R2,IOVPAR,EM,NE)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL IOVPAR
      DIMENSION PARM(13),EM(NE,5)
C
      ID=IOVPAR(I,J,R2,PARM)
      DO 11 L=1,5
      DO 10 M=1,5
10    EM(L,M)=0.0D0
11    EM(L,L)=PARM(11)
      RETURN
      END
