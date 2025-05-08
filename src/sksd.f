      SUBROUTINE SKSD(X,X2,I,J,R2,IOVPAR,EM,EMT,NE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(6),X2(6),PARM(13),EM(NE,5),EMT(NE,5),ES(5)
      EXTERNAL IOVPAR
C
      R3=SQRT(3.0D0)
      D4=X2(3)-0.5D0*(X2(1)+X2(2))
      D5=X2(1)-X2(2)
      ID=IOVPAR(I,J,R2,PARM)
C
      DO 10 L=1,3
   10 ES(L)=R3*X(L)*X(L+1)
      ES(4)=0.5D0*R3*D5
      ES(5)=D4
      DO 11 L=1,5
      EM(1,L)=ES(L)*PARM(8)
   11 EMT(L,1)=EM(1,L)
      RETURN
      END
