      SUBROUTINE SKPD(X,X2,I,J,R2,IOVPAR,EM,EMT,NE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(6),X2(6),PARM(13),EM(NE,5),EMT(NE,5),EPD(13,2),DM(15)
      EXTERNAL IOVPAR
C
      R3=SQRT(3.0D0)
      D3=X2(1)+X2(2)
      D4=X2(3)-0.5D0*D3
      D5=X2(1)-X2(2)
      D6=X(1)*X(2)*X(3)
      ID=IOVPAR(I,J,R2,PARM)
      DO 10 L=1,3
      EPD(L,1)=R3*X2(L)*X(L+1)
      EPD(L,2)=X(L+1)*(1.0D0-2*X2(L))
      EPD(L+4,1)=R3*X2(L)*X(L+2)
      EPD(L+4,2)=X(L+2)*(1.0D0-2*X2(L))
      EPD(L+7,1)=0.5D0*R3*X(L)*D5
   10 EPD(L+10,1)=X(L)*D4
      EPD(4,1)=R3*D6
      EPD(4,2)=-2*D6
      EPD(8,2)=X(1)*(1.0D0-D5)
      EPD(9,2)=-X(2)*(1.0D0+D5)
      EPD(10,2)=-X(3)*D5
      EPD(11,2)=-R3*X(1)*X2(3)
      EPD(12,2)=-R3*X(2)*X2(3)
      EPD(13,2)=R3*X(3)*D3
      DO 11 L=1,15
   11 DM(L)=0.0D0
      DO 12 M=1,2
      DM(1)=DM(1)+EPD(1,M)*PARM(M+3)
      DM(2)=DM(2)+EPD(6,M)*PARM(M+3)
      DM(3)=DM(3)+EPD(4,M)*PARM(M+3)
      DM(5)=DM(5)+EPD(2,M)*PARM(M+3)
      DM(6)=DM(6)+EPD(7,M)*PARM(M+3)
      DM(7)=DM(7)+EPD(5,M)*PARM(M+3)
      DM(9)=DM(9)+EPD(3,M)*PARM(M+3)
      DO 12 L=8,13
   12 DM(L+2)=DM(L+2)+EPD(L,M)*PARM(M+3)
      DM(4)=DM(3)
      DM(8)=DM(3)
      DO 13 IR=1,5
      DO 13 IS=1,3
      K=3*(IR-1)+IS
      EMT(IR,IS)=-DM(K)
   13 EM(IS,IR)=DM(K)
      RETURN
      END
