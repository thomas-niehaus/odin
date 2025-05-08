      SUBROUTINE SKPP(X,X2,I,J,R2,IOVPAR,EM,NE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(6),DM(6),X2(6),PARM(13),EM(NE,3),EPP(6)
      EXTERNAL IOVPAR
C
      ID=IOVPAR(I,J,R2,PARM)
      DO 10 L=1,3
      EPP(L)=X2(L)
   10 EPP(L+3)=X(L)*X(L+1)
      DO 11 L=1,3
      HP=EPP(L)
   11 DM(L)=HP*PARM(6)+(1.0D0-HP)*PARM(7)
      DO 12 L=4,6
   12 DM(L)=EPP(L)*(PARM(6)-PARM(7))
      DO 13 IR=1,3
      DO 13 IS=1,IR
      II=IR-IS
      K=3*II-(II*(II-1))/2+IS
      EM(IS,IR)=DM(K)
   13 EM(IR,IS)=DM(K)
      RETURN
      END
