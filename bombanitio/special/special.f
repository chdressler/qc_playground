      subroutine SPECIALFUNC(K,X,F1,F2,F3,F4,F5)
c
c Evaluates special functions Fk(X) involved in nuclear-attraction and
c electron-repulsion integrals. 
c Only the first K functions (1 <= K <= 5) are evaluated.
c
      external ERROR
c
c      , intent(out)
      integer K
      double precision X,F1,F2,F3,F4,F5,
     $xthresh(5),cost,xpower(9),erf
cf2py intent(in)  K
cf2py intent(in)  X
cf2py intent(out)  F1,F2,F3,F4,F5
c
      data xthresh /.0001d0,.1d0,.25d0,.65d0,.75d0/
      cost=(2.d0/dasin(1.d0))**.5d0
      xpower(2)=X**2.d0
      if (X .ge. xthresh(1)) call ERROR(X,erf)
      if (K .gt. 1 .and. X .lt. xthresh(K)) then
       xpower(4)=X**4.d0
       xpower(6)=X**6.d0
       xpower(8)=X**8.d0
      endif
      do 10, i=2,K
       xpower(2*i-1)=X**dble(2*i-1)
10     continue
c
      if (X .lt. xthresh(1)) then
       F1=cost*(1.d0-xpower(2)/3.d0)
      else
       F1=erf/X
      endif
      if (K .eq. 1) goto 20
c
      if (X .lt. xthresh(2)) then
       F2=cost*(-2.d0/3.d0+.4d0*xpower(2)-xpower(4)/7.d0+
     $xpower(6)/27.d0-xpower(8)/132.d0)
      else
       F2=(cost*X*exp(-xpower(2))-erf)/xpower(3)
      endif
      if (K .eq. 2) goto 20
c
      if (X .lt. xthresh(3)) then
       F3=cost*(.8d0-4.d0/7.d0*xpower(2)+2.d0/9.d0*xpower(4)-
     $2.d0/33.d0*xpower(6)+xpower(8)/78.d0)
      else
       F3=(3.d0*erf-cost*exp(-xpower(2))*(3.d0*X+2.d0*xpower(3)))/
     $xpower(5)
      endif
      if (K .eq. 3) goto 20
c
      if (X .lt. xthresh(4)) then
       F4=cost*(-8.d0/7.d0+8.d0/9.d0*xpower(2)-4.d0/11.d0*xpower(4)+
     $4.d0/39.d0*xpower(6)-xpower(8)/45.d0)
      else
       F4=(cost*exp(-xpower(2))*(15.d0*X+10.d0*xpower(3)+4.d0*
     $xpower(5))-15.d0*erf)/xpower(7)
      endif
      if (K .eq. 4) goto 20
c
      if (X .lt. xthresh(5)) then
       F5=cost*(16.d0/9.d0-16.d0/11.d0*xpower(2)+8.d0/13.d0*xpower(4)-
     $8.d0/45.d0*xpower(6)+2.d0/51.d0*xpower(8))
      else
       F5=(105.d0*erf-cost*exp(-xpower(2))*(105.d0*X+70.d0*xpower(3)+
     $28.d0*xpower(5)+8.d0*xpower(7)))/xpower(9)
      endif
c
20    return
      end
c


      subroutine ERROR(X,ERR)
C
C     =========================================
C     Purpose: Compute error function erf(x)
C     Input:   x   --- Argument of erf(x)
C     Output:  ERR --- erf(x)
C     =========================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EPS=1.0D-15
      PI=3.141592653589793D0
      X2=X*X
      IF (DABS(X).LT.3.5D0) THEN
         ER=1.0D0
         R=1.0D0
         DO 10 K=1,50
            R=R*X2/(K+0.5D0)
            ER=ER+R
            IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
10       CONTINUE
15       C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
         ERR=C0*ER
      ELSE
         ER=1.0D0
         R=1.0D0
         DO 20 K=1,12
            R=-R*(K-0.5D0)/X2
20          ER=ER+R
         C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
         ERR=1.0D0-C0*ER
         IF (X.LT.0.0) ERR=-ERR
      ENDIF
      RETURN
      END

