*
*-----------------------------------------------------------------------
* On the next-to-next-to-leading order QCD corrections
* to heavy-quark production in deep-inelastic
*
* H. Kawamura, N. Lo Presti, S. Moch and A. Vogt
*-----------------------------------------------------------------------
*
* The threshold approximations for the gluon coefficient function c_{2,g}
* at two loops:
* cgt2 in Eq. (3.18)
*
*-----------------------------------------------------------------------
* ..The 2-loop threshold approximation to the gluon c_{2,g}^2
*   see Eq. (3.18)
*   power ln(mu)^0
*
      real*8 function cgt2(eta,xi,nf)
      implicit none
      real*8 eta,xi,beta,lnbeta,pi
      real*8 nlotconst,born_t,born_l
      integer nf
      parameter(pi = 3.14159265359d0)

      beta = dsqrt(eta/(1d0 + eta))
      lnbeta = dlog(beta)

      cgt2 = 1d0 * (
     ~  + 3.607744112d+0/beta**2
     ~  + 2.894748951d+2/beta
     ~  - 1.579136704d+2/beta*lnbeta**2
     ~  + 4.187651872d+1/beta*lnbeta
     ~  + 6.149252757d-1/beta*nf
     ~  - 4.386490844d+0/beta*nf*lnbeta
     ~  + 1.152d+3*lnbeta**4
     ~  - 1.672966687d+3*lnbeta**3
     ~  - 3.328893061d+3*lnbeta**2
     ~  + 2.262291573d+3*lnbeta
     ~  + 4.266666666d+1*nf*lnbeta**3
     ~  - 8.024907466d+1*nf*lnbeta**2
     ~  + 4.706482425d+1*nf*lnbeta
     ~  - 3.289868133d+0*nlotconst(xi)/beta
     ~  + 4.8d+1*nlotconst(xi)*lnbeta**2
     ~  - 2.018680599d+1*nlotconst(xi)*lnbeta
     ~ )

      return
      end
*
*-----------------------------------------------------------------------
* ..The 2-loop threshold approximation to the gluon c_{2,g}^2
*   see Eq. (3.18)
*   power ln(mu)^1
*
      real*8 function cgt2br1(eta,xi)
      implicit none
      real*8 eta,xi,beta,lnbeta,pi
      real*8 nlotconst,nlotbarconst,born_t,born_l
      parameter(pi = 3.14159265359d0)

      beta = dsqrt(eta/(1d0 + eta))
      lnbeta = dlog(beta)

      cgt2br1 = 1d0*(
     ~  + 4.934802200d+1/beta
     ~  + 7.895683520d+1/beta*lnbeta
     ~  - 1.152d+3*lnbeta**3
     ~  + 1.179777919d+2*lnbeta**2
     ~  + 2.222515190d+3*lnbeta
     ~  - 2.4d+1*nlotconst(xi)*lnbeta
     ~  - 3.289868133d+0*nlotbarconst(xi)/beta
     ~  + 4.8d+1*nlotbarconst(xi)*lnbeta**2
     ~  - 2.018680599d+1*nlotbarconst(xi)*lnbeta
     ~ )

      return
      end
*
*-----------------------------------------------------------------------
* ..The 2-loop threshold approximation to the gluon c_{2,g}^2
*   see Eq. (3.18)
*   power ln(mu)^2
*
      real*8 function cgt2br2(eta,xi)
      implicit none
      real*8 eta,xi,beta,lnbeta,pi
      real*8 nlotconst,nlotbarconst,born_t,born_l
      parameter(pi = 3.14159265359d0)

      beta = dsqrt(eta/(1d0 + eta))
      lnbeta = dlog(beta)

      cgt2br2 = 1d0*(
     ~  + 2.88d+2*lnbeta**2
     ~  + 2.912527760d+2*lnbeta
     ~  - 2.4d+1*nlotbarconst(xi)*lnbeta
     ~ )

      return
      end

*
*
* ------------------------------------------------------------------------
* ..The exact expression for the constant term O(beta^0) of the
*   gluon coefficient function at 1 loop
*   see Eq.(3.10)
*   power ln(mu)^0
*   extracted from:
*   E. Laenen, S. Riemersma, J. Smith, W.L. van Neerven, Nucl.Phys. B392 (1993) 162-228
*
      real*8 function nlotconst(xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision ln2,LOGA,LOGB,LOGC,LOGD
      double precision inv2pxi,inv4pxi
      parameter ( z2 = 1.6449 34066 84822 64365 D0)
      parameter(pi = 3.14159265359d0)

      real*8 dilog

      ca = 3.d0
      cf = 4.d0/3.d0

      ln2 = dlog(2d0)
      YSQ = DSQRT(1d0+4d0/xi)
      LOGA = DLOG(YSQ-1d0)
      LOGB = DLOG(XI)
      LOGC = DLOG(4d0+XI)
      LOGD = DLOG(2d0+XI)

      inv2pxi = (2d0 + xi)**(-1)
      inv4pxi = (4d0 + xi)**(-1)

      consta = 56d0 - pi**2
     #  + 12d0*LOGA*YSQ**(-1)
     #  + 4d0*LOGA**2
     #  + LOGD**2
     #  + 4d0*LOGC
     #  - 4d0*LOGC**2
     #  + 6d0*LOGB*YSQ**(-1)
     #  + 4d0*LOGB*LOGA
     #  + LOGB**2
     #  - 4d0*LN2
     #  - 12d0*LN2*YSQ**(-1)
     #  - 8d0*LN2*LOGA
     #  - 2d0*LN2*LOGD
     #  + 8d0*LN2*LOGC
     #  - 4d0*LN2*LOGB
     #  + LN2**2
     #  + 2d0*DILOG( - 2d0*inv2pxi)
     #  + 36d0*inv2pxi
     #  - 6d0*inv2pxi*XI
     #  + 24d0*inv2pxi*LOGC
     #  - 4d0*inv2pxi*LOGC*XI
     #  - 24d0*inv2pxi*LN2
     #  + 4d0*inv2pxi*LN2*XI
     #  + 16d0*inv2pxi**(+2)*LOGC
     #  - 12d0*inv2pxi**(+2)*LOGC*XI
     #  + 6d0*inv2pxi**(+2)*LOGC*XI**2
     #  - 16d0*inv2pxi**(+2)*LN2
     #  + 12d0*inv2pxi**(+2)*LN2*XI
     #  - 6d0*inv2pxi**(+2)*LN2*XI**2
     #  - 128d0*inv2pxi**(+2)*inv4pxi*LOGC
     #  + 128d0*inv2pxi**(+2)*inv4pxi*LN2
     #  - 160d0*inv2pxi*inv4pxi
     #  - 64d0*inv2pxi*inv4pxi*LOGC
     #  + 64d0*inv2pxi*inv4pxi*LN2
     #  + 128d0*inv2pxi*inv4pxi**(+2)
     #  - 48d0*inv4pxi
     #  - 2d0*inv4pxi*pi**2
     #  - 8d0*inv4pxi*LOGA**2
     #  - 8d0*inv4pxi*LOGB*LOGA
     #  - 2d0*inv4pxi*LOGB**2
     #  + 16d0*inv4pxi*LN2*LOGA
     #  + 8d0*inv4pxi*LN2*LOGB
     #  - 8d0*inv4pxi*LN2**2
     #  + 64d0*inv4pxi**(+2)

      constf =        - 12d0
     #  - 4d0*Z2
     #  - 24d0*LOGA*YSQ**(-1)
     #  - 8d0*LOGA**2
     #  + 2d0*LOGD**2
     #  - 8d0*LOGC
     #  - 12d0*LOGB*YSQ**(-1)
     #  - 8d0*LOGB*LOGA
     #  - 2d0*LOGB**2
     #  + 8d0*LN2
     #  + 24d0*LN2*YSQ**(-1)
     #  + 16d0*LN2*LOGA
     #  - 4d0*LN2*LOGD
     #  + 8d0*LN2*LOGB
     #  - 6d0*LN2**2
     #  + 4d0*DILOG( - 2d0*inv2pxi)
     #  - 60d0*inv2pxi
     #  - 6d0*inv2pxi*XI
     #  - 16d0*inv2pxi*LOGC
     #  - 12d0*inv2pxi*LOGC*XI
     #  + 16d0*inv2pxi*LN2
     #  + 12d0*inv2pxi*LN2*XI
     #  - 64d0*inv2pxi**(+2)*LOGC
     #  + 56d0*inv2pxi**(+2)*LOGC*XI
     #  + 14d0*inv2pxi**(+2)*LOGC*XI**2
     #  + 64d0*inv2pxi**(+2)*LN2
     #  - 56d0*inv2pxi**(+2)*LN2*XI
     #  - 14d0*inv2pxi**(+2)*LN2*XI**2
     #  + 256d0*inv2pxi**(+2)*inv4pxi*LOGC
     #  - 256d0*inv2pxi**(+2)*inv4pxi*LN2
     #  + 360d0*inv2pxi*inv4pxi
     #  + 128d0*inv2pxi*inv4pxi*LOGC
     #  - 128d0*inv2pxi*inv4pxi*LN2
     #  - 544d0*inv2pxi*inv4pxi**(+2)
     #  + 44d0*inv4pxi
     #  + 40d0*inv4pxi*Z2
     #  + 48d0*inv4pxi*LOGA**2
     #  - 8d0*inv4pxi*LOGD**2
     #  + 48d0*inv4pxi*LOGB*LOGA
     #  + 12d0*inv4pxi*LOGB**2
     #  - 96d0*inv4pxi*LN2*LOGA
     #  + 16d0*inv4pxi*LN2*LOGD
     #  - 48d0*inv4pxi*LN2*LOGB
     #  + 40d0*inv4pxi*LN2**2
     #  - 16d0*inv4pxi*DILOG( - 2d0*inv2pxi)
     #  - 272d0*inv4pxi**(+2)

      nlotconst=consta*ca+constf*cf

      return
      end
*
* ------------------------------------------------------------------------
* ..The exact expression for the constant term O(beta^0) of the
*   gluon coefficient function at 1 loop
*   see Eq.(3.11)
*   power ln(mu)^1
*   extracted from:
*   E. Laenen, S. Riemersma, J. Smith, W.L. van Neerven, Nucl.Phys. B392 (1993) 162-228
*
      real*8 function nlotbarconst(xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision ln2,LOGC
      parameter(pi = 3.14159265359d0)

      ca = 3.0d0
      tf = 0.5d0

      ln2 = dlog(2d0)
      LOGC = DLOG(4d0+XI)

      nlotbarconst = 4d0*ca*(2d0+LOGC-2d0*ln2)
      nlotbarconst = nlotbarconst - 4.d0/3.d0*tf

      return
      end
*
* ------------------------------------------------------------------------
* ..The Pade estimate for the constant term O(beta^0) of the
*   gluon coefficient function at 2 loops
*   power ln(mu)^0
*
      real*8 function cgt2pade(eta,xi)
      implicit none
      real*8 ca,ln2,eta,xi,pi
      real*8 nlotconst,born_t,born_l
      parameter(pi = 3.14159265359d0)

      ca = 3.d0
      ln2 = dlog(2d0)

      cgt2pade = 1d0*(
     ~  nlotconst(xi)+36d0*ca*ln2**2-60d0*ca*ln2
     ~ )**2

      return
      end
*
* ------------------------------------------------------------------------
* ..The transverse Born gluon coefficient function
*
      double precision function born_t(eta,xi)
      implicit none
      double precision eta, xi, pi,tf
      parameter(pi = 3.14159265359d0)

      tf = 0.5d0

      born_t = 0.5d0*pi*tf*(1.d0 + eta + 0.25d0*xi)**(-3)*
     #         (-2.d0*((1.d0 + eta - 0.25d0*xi)**2 + eta + 1.d0)*
     #         dsqrt(eta/(1.d0 + eta)) + (2.d0*(1.d0 + eta)**2 +
     #         0.125d0*xi**2 + 2.d0*eta + 1.d0)*
     #         dlog((dsqrt(1.d0 + eta) + dsqrt(eta))/
     #              (dsqrt(1.d0 + eta) - dsqrt(eta))))
      return
      end
*
* ------------------------------------------------------------------------
* ..The longitudinal Born gluon coefficient function
*
      double precision function born_l(eta,xi)
      implicit none
      double precision eta, xi, pi, tf
      parameter(pi = 3.14159265359d0)

      tf = 0.5d0

      born_l = 0.5d0*pi*tf*xi*(1.d0 + eta + 0.25d0*xi)**(-3)*
     #         (2.d0*dsqrt(eta*(1.d0 + eta)) -
     #         dlog((dsqrt(1.d0 + eta) + dsqrt(eta))/
     #               (dsqrt(1.d0 + eta) - dsqrt(eta))))
      return
      end
*
*
* File: dilog.f
*
* ..A routine for the dilogarithm, by Jos Vermaseren
* ------------------------------------------------------------------------
*
       double precision function dilog(x)
       implicit double precision  (a-z)
       dimension b(8)
       integer ncall
       data ncall/0/,pi6/1.644934066848226d+00/,een,vier/1.d+00,.25d+00/
       ncall = 0
       if(ncall.eq.0)go to 2
1      if(x.lt.0)go to 3
       if(x.gt.0.5)go to 4
       z=-dlog(1.-x)
7      z2=z*z
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       if(x.gt.een)dilog=-dilog-.5*u*u+2.*pi6
       return
2      b(1)=een
       b(2)=een/36.
       b(3)=-een/3600.
       b(4)=een/211680.
       b(5)=-een/(30.*362880.d+00)
       b(6)=5./(66.*39916800.d+00)
       b(7)=-691./(2730.*39916800.d+00*156.)
       b(8)=een/(39916800.d+00*28080.)
       ncall=1
       go to 1
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier-u*(z+.5*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een-u)+z2*vier+pi6
       if(x.gt.een)dilog=-dilog-.5*z*z+pi6*2.
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1./x
       if(x.gt.2.)go to 11
       z=dlog(x)
       y=1.-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.-xx)
       go to 7
20     dilog=pi6
       return
       end
