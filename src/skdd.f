      SUBROUTINE SKDD(X,X2,I,J,R2,IOVPAR,EM,NE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(6),X2(6),PARM(13),EM(NE,5),E(15,3),DM(15),DD(3)
      EXTERNAL IOVPAR
C
      R3=SQRT(3.0D0)
      D3=X2(1)+X2(2)
      D4=X2(3)-0.5D0*D3
      D5=X2(1)-X2(2)
      ID=IOVPAR(I,J,R2,PARM)
      DO 3 L=1,3
      E(L,1)=X2(L)*X2(L+1)
      E(L,2)=X2(L)+X2(L+1)-4*E(L,1)
      E(L,3)=X2(L+2)+E(L,1)
3     E(L,1)=3*E(L,1)
      E(4,1)=D5*D5
      E(4,2)=D3-E(4,1)
      E(4,3)=X2(3)+0.25D0*E(4,1)
      E(4,1)=0.75D0*E(4,1)
      E(5,1)=D4*D4
      E(5,2)=3*X2(3)*D3
      E(5,3)=D3*D3*0.75D0
      DD(1)=X(1)*X(3)
      DD(2)=X(2)*X(1)
      DD(3)=X(3)*X(2)
      DO 4 L=1,2
      E(L+5,1)=3*X2(L+1)*DD(L)
      E(L+5,2)=DD(L)*(1.0D0-4*X2(L+1))
4     E(L+5,3)=DD(L)*(X2(L+1)-1.0D0)
      E(8,1)=DD(1)*D5*1.5D0
      E(8,2)=DD(1)*(1.0D0-2*D5)
      E(8,3)=DD(1)*(0.5D0*D5-1.0D0)
      E(9,1)=D5*0.5D0*D4*R3
      E(9,2)=-D5*X2(3)*R3
      E(9,3)=D5*0.25D0*(1.0D0+X2(3))*R3
      E(10,1)=X2(1)*DD(3)*3
      E(10,2)=(0.25D0-X2(1))*DD(3)*4
      E(10,3)=DD(3)*(X2(1)-1.0D0)
      E(11,1)=1.5D0*DD(3)*D5
      E(11,2)=-DD(3)*(1.0D0+2*D5)
      E(11,3)=DD(3)*(1.0D0+0.5D0*D5)
      E(13,3)=0.5D0*D5*DD(2)
      E(13,2)=-2*DD(2)*D5
      E(13,1)=E(13,3)*3
      E(12,1)=D4*DD(1)*R3
      E(14,1)=D4*DD(3)*R3
      E(15,1)=D4*DD(2)*R3
      E(15,2)=-2*R3*DD(2)*X2(3)
      E(15,3)=0.5D0*R3*(1.0D0+X2(3))*DD(2)
      E(14,2)=R3*DD(3)*(D3-X2(3))
      E(14,3)=-R3*0.5D0*DD(3)*D3
      E(12,2)=R3*DD(1)*(D3-X2(3))
      E(12,3)=-R3*0.5D0*DD(1)*D3
      DO 11 L=1,15
      DM(L)=0.0D0
      DO 11 M=1,3
   11 DM(L)=DM(L)+E(L,M)*PARM(M)
      DO 12 IR=1,5
      DO 12 IS=1,IR
      II=IR-IS
      K=5*II-(II*(II-1))/2+IS
      EM(IR,IS)=DM(K)
   12 EM(IS,IR)=DM(K)
      RETURN
      END
