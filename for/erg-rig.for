       program ERG_RIG

c........+.........+.........+.........+.........+.........+.........+..
c      Energy to Rigidity conversion example
c........+.........+.........+.........+.........+.........+.........+..
c      \/ Example of the use of the Energy to Rigidity conversion program
c         Convert geomagnetic cutoff rigidity to proton cutoff energy
c                 need to define proton atomic number, charge and mass
c                 rigidity in often used in units of GV;
c       Subroutine azrgeg needs   rigidity  in units of MV
c       Subroutine azrgeg returns  energy   in units of MeV per nucleon
c........+.........+.........+.........+.........+.........+.........+..
c     Programmed by Don F. Smart  (sssrc@msn.com)
c     Note - programming adheres to convention that variables beginning
c            with i, j, k, l, m, n   are integer values,
c            variables beginning with c are character variables
c            all other variables are real*8
c........+.........+.........+.........+.........+.........+.........+..
c
      implicit integer (i-n)
      implicit REAl*8 (a-b)
      implicit REAl*8 (d-h)
      implicit REAl*8 (o-z)
c
c........+.........+.........+.........+.........+.........+.........+..

      open (7,file='erg-rig.txt', status='unknown')

c........+.........+.........+.........+.........+.........+.........+..
c
c     Define atomic number, atomic change and rest mass for oxygen16
c     (Note, any element or isotope can be specified)
c........+.........+.........+.........+.........+.........+.........+..
c
      na = 16
      nz = 8
      pamu = 16.00
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Fist demonstration  does Energy to Rigidity from 1 to 100 MeV
c........+.........+.........+.........+.........+.........+.........+..
c
      Do 100 i=1,100
         rpmv = 0.0
         epnmev = float (i)
         call azrgeg (na,nz,pamu,rpmv,epnmev,beta)

         write (7, 1070)  rpmv,epnmev

  100 continue
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Second demonstration  does Energy to Rigidity 
c                              from 100 to 1000 MeV in 10 MeV increments
c........+.........+.........+.........+.........+.........+.........+..
c
      Do 110 i=100,1000,10
         rpmv = 0.0
         epnmev = float (i)
         call azrgeg (na,nz,pamu,rpmv,epnmev,beta)

         write (7, 1070)  rpmv,epnmev
  110 continue
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Second demonstration  does Energy to Rigidity 
c                              from 1000 to 10000 MeV in 100 MeV increments
c........+.........+.........+.........+.........+.........+.........+..
c
      Do 120 i=1000,10000,100
         rpmv = 0.0
         epnmev = float (i)
         call azrgeg (na,nz,pamu,rpmv,epnmev,beta)

         write (7, 1070)  rpmv,epnmev
  120 continue

 1070  format (1f10.1, 1f10.3)

      stop
      end
      subroutine azrgeg (na,nz,pamu,rigin,epn,beta)
c
c........+.........+.........+.........+.........+.........+.........+..
c     subroutine to convert rigidity to energy and visa versa
cLast mod 17 March 94
c........+.........+.........+.........+.........+.........+.........+..
c     Programmed by Don F. Smart  (sssrc@msn.com)
c     Note - programming adheres to convention that variables beginning
c            with i, j, k, l, m, n   are integer values,
c            variables beginning with c are character variables
c            all other variables are real*4
c........+.........+.........+.........+.........+.........+.........+..
c
      implicit integer (i-n)
      implicit REAl*8 (a-b)
      implicit REAl*8 (d-h)
      implicit REAl*8 (o-z)
c
c........+.........+.........+.........+.........+.........+.........+..
c
c     write (*,6010) na,nz,pamu,rigin,epn                               ! diag
c6010 format (' azrgeg1 ',2i5,3f15.3)                                   ! diag
c
c........+.........+.........+.........+.........+.........+.........+..
c     Check, if na, nz, or pamu not specified, put in default for protons
c            epamu is rest mass energy per atomic mass unit
c........+.........+.........+.........+.........+.........+.........+..
c
      if (pamu.le.0.0)  pamu = 1.0081451
      if (na.le.0)      na = 1
      if (nz.le.0)      nz = 1
c
      epamu = 931.141
c
      anuc = na
      zcharg = nz
      rsmspn = (pamu/anuc)*epamu
c
      trig = rigin
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ if trig .le. 0.0    do energy   to rigidity conversion
c        if trig .gt. 0.0    do rigidity to energy conversion
c........+.........+.........+.........+.........+.........+.........+..
c
      if (trig.le.0.0)  then
c                                          Energy to Rigidity conversion
         gmaeg = (anuc*epn+anuc*rsmspn)/(anuc*rsmspn)
         gmaegg = (epn+rsmspn)/rsmspn
         rigin = dsqrt(gmaeg*gmaeg-1.0)*rsmspn*anuc/zcharg
         relgma = gmaeg
      else
c                                          Rigidity to Energy conversion
         gmarg = dsqrt(((rigin*zcharg)/(rsmspn*anuc))**2+1.0)
         epn = (gmarg-1.0)*rsmspn
         relgma = gmarg
      endif
c
c     write (*,6020) anuc,zcharg,rsmspn,rigin,epn,gmaeg,gmaegg,gmarg    ! diag
c     write (*,6030) na,nz,pamu,rigin,epn                               ! diag
c6020 format (' azrgeg2 ',8f10.3)                                       ! diag
c6030 format (' azrgeg3 ',2i5,3f15.5)                                   ! diag
c
      beta = dsqrt(1.0-1.0/(relgma*relgma))
c
      return
c
c          beta is v/c (speed as fraction of light speed)
c          energy in MeV
c          epamu is mass-energy conversion = 931.141 MeV per amu
c          epn is kinetic energy per nucleon
c          na is atomic number
c          nz is charge
c          pamu is rest mass in physical atomic mass units
c          relgam is relativistic factor 'gamma'
c          rigidity in mv
c          rsmspn is rest mass per nucleon in MeV
c
      end
