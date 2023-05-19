      PROGRAM TJALLMAG
C
C........+.........+.........+.........+.........+.........+.........+..
c     Mod 08 Dec 00  Use ALLMAG for magnetic field calculations
C........+.........+.........+.........+.........+.........+.........+..
C     Multi-platform COSMIC-RAY TRAJECTORY PROGRAM
C     FORTRAN 77  transportable version (no time functions)
C     Read in control card; LAT, LON, RIG, ZENITH, AZIMUTH, DELPC, INDO
C          Then calculate        INDO    trajectories
C               Starting at      PC
C               Incrementing at  DELPC   intervals
C     Includes conversion from Geodetic to Geocentric coordinates
C     Includes re-entrant albedo calculations
C     Uses subroutine SINGLE to do trajectory calculations
C     Magnetic field - as determined by the NASA ALLMAG routine
C........+.........+.........+.........+.........+.........+.........+..
C     Restrictions: Cannot run over N or S pole; will get BETA blowup
C........+.........+.........+.........+.........+.........+.........+..
C     Mod History
CLast Mod 22 Dec 00  Use NASA ALLMAG for magnetic field calculations
C     Mod 21 Dec 00  Make all intrinsic function double precision for PC
C     Mod 20 Dec 00  Insert 8 character format 1000 with AZ & ZE
C     Mod 17 Feb 99  if (ymax.lt.6.6) IFATE = 3
C     Mod 17 Feb 99  set limit to 600000
C     Mod    Aug 97  Adjust step size to minimize beta problems
C     Mod    Jan 97  High latitude step size adjust, introduce AHLT
C     Mod    Jun 96  EDIF limit set to 1.0e-5
C     Mod    Jun 96  IERRPT formats, Boundary and look ahead 
C     Mod    May 96  Make all FORTRAN coding upper case
C     Mod    Feb 96  Standard reference TJ1V line check
C     **************************************************************
C          Timing estimates base on COMPAQ Digital FORTRAN
C     Will run on PIII PC at 850 MHZ        55000 steps/sec (Real*8)
C     Will run on PIII PC at 700 MHZ        39000 steps/sec (Real*8)
C     Will run on PIII PC at 550 MHZ        32000 steps/sec (Real*8)
C     Will run on PII  PC at 400 MHZ        23000 steps/sec (Real*8)
**************************************************************
C     *  TAPE*       Monitor program operation
C     *  TAPE1       Trajectory control cards
C     *  TAPE7       80  character card image   output
C     *  TAPE8       132 character line printer output
C     *  TAPE16      Diagnostic output for trouble shooting
C     *              Normally turned off (open statement commented out)
C     **************************************************************
C
C........+.........+.........+.........+.........+.........+.........+..
C     Programmer - Don F. Smart; FORTRAN77
C     Note - The programming adheres to the conventional FORTRAN
C            default standard that variables beginning with
C            'i','j','k','l','m',or 'n' are integer variables
C            Variables beginning with "c" are character variables
C            All other variables are real
C........+........+.........+.........+.........+.........+.........+..
C            Do not mix different type variables in same common block
C            Some computers do not allow this
C........+.........+.........+.........+.........+.........+.........+..
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL * 8(A-B)
      IMPLICIT REAL * 8(D-H)
      IMPLICIT REAL * 8(O-Z)
C
C........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /WRKVLU/ F(6),Y(6),ERAD,EOMC,VEL,BR,BT,BP,B
      COMMON /WRKTSC/ TSY2,TCY2,TSY3,TCY3
      COMMON /TRIG/   PI,RAD,PIO2
      COMMON /GEOID/  ERADPL, ERECSQ
      COMMON /SNGLR/  SALT,DISOUT,GCLATD,GDLATD,GLOND,GDAZD,GDZED,
     *                RY1,RY2,RY3,RHT,TSTEP
      COMMON /SNGLI/  LIMIT,NTRAJC,IERRPT
C
C........+.........+.........+.........+.........+.........+.........+..
C
      OPEN (1, FILE='TAPE1', STATUS='OLD')
      OPEN (7, FILE='TAPE7', STATUS='UNKNOWN')
      OPEN (8, FILE='TAPE8', STATUS='UNKNOWN')
      OPEN (16,FILE='TAPE16',STATUS='UNKNOWN')
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ User defined program control
C........+.........+.........+.........+.........+.........+.........+..
C
      FSTEP = 4.0E08
      LIMIT  = 600000
C
C........+.........+.........+.........+.........+.........+.........+..
C     /\ FSTEP  is total number of steps before run is terminated
C        LIMIT  is       max number of steps before trajectory declared F
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Define program constants
C........+.........+.........+.........+.........+.........+.........+..
C        DISOUT is       radial distance for trajectory termination
C        ERAD   is       average earth radius
C        NTRAJC is       number of trajectory computed in this run
C        RHT    is       top of atmosphere for re-entrant trajectory
C        TSTEP  is       number of steps executed in this run
C........+.........+.........+.........+.........+.........+.........+..
C
      NTRAJC = 0
      TSTEP = 0.0
C
      DISOUT = 25.0
      ERAD   = 6371.2
      RHT    = 20.0
      VEL    = 2.99792458E5/ERAD
C
C........+.........+.........+.........+.........+.........+.........+..
C      "VEL" is light velocity in earth radii per second
C      Light speed defined as    299792458 m/s
C       Ref: E. R. Cohne AND B. N. Taylor, "The Fundamental Physical
C            Constants, Physics Today P11, August 1987.
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Define essential trigonometric values
C........+.........+.........+.........+.........+.........+.........+..
C
C     PI   = ACOS(-1.0)                                                 !Sngl
      PI   = DACOS(-1.0D0)                                              !Dbl
      RAD  = 180.0/PI
      PIO2 = PI/2.0
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ TAPE1 must contain trajectory control cards
C        Terminate program if no data on TAPE1 file
C        Terminate if EOF encountered
C        Terminate if negative data found on input file
C        Terminate if bad      data found on input file
C........+.........+.........+.........+.........+.........+.........+..
C
  100 READ (1,1010,IOSTAT=IOSTAT,ERR=120,END=110) GDLATD,GLOND,PC,
     *     GDZED,GDAZD,DELPC,INDO,IERRPT,INDEX
 1010 FORMAT (BZ,6F8.2,3I8)
C
  110 CONTINUE
      IF (IOSTAT.LT.0) THEN
         WRITE (*,1020)
         GO TO 150
      ENDIF
 1020 FORMAT (' END OF FILE ON TAPE 1 (DATA INPUT)')

  120 IF (IOSTAT.GT.0) THEN
         WRITE (*,1030) IOSTAT,GDLATD,GLOND,PC,DELPC,
     *                  INDO,IERRPT,INDEX
         GO TO 150
      ENDIF
 1030 FORMAT (' ERROR ON DATA INPUT FILE (TAPE1), IOSTAT =',I5/
     *   4F8.3,3I8)
C
      IF (PC.LE.0) THEN
         WRITE (*,1040)
         GO TO 150
      ENDIF
 1040 FORMAT (' END OF DATA INPUT (NEGATIVE VALUE READ IN)')
C
      WRITE (*,1050) GDLATD,GLOND,PC,GDZED,GDAZD,DELPC,INDO,IERRPT,INDEX
 1050 FORMAT (' TAPE 1  ',6F7.2,3I6)
C
C.......+.........+.........+.........+.........+.........+.........+..
C     \/ Start at top of atmosphere (20 km above surface of oblate earth)
C        Coding is relic of past when ISALT was read in
C.......+.........+.........+.........+.........+.........+.........+..
C
      ISALT = 0
      IF (ISALT.LE.0) SALT = 20.0
      IF (ISALT.GT.0) SALT = ISALT
C
      KNT = 0
      IDELPC = DELPC*1000.0+0.0001
      INDXPC = PC*1000.0+0.0001
C
C.......+.........+.........+.........+.........+.........+.........+..
C     For trajectories from Earth
C         convert from Geodetic coordinates to Geocentric coordinates
C                Geodetic   coordinates used for input
C                GEOCENTRIC coordinates used for output
C        coordinates are placed in common block /MNSINGL/
C        All calculation are done in  Geocentric coordinates!
C     \/ Conversion from Geodetic to Geocentric coordinates
C.......+.........+.........+.........+.........+.........+.........+..
C
      CALL GDGC (TCD, TSD)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Remember positron of initial point on trajectory
C                 in Geocentric coordinates
C
C        Y(1) is distance in earth radii from geocenter
C             Start with height above geoid and convert to earth radii
C                   The initial values of Y(1), Y(2) AND Y(3) are
C                   calculated in subroutine GDGC
C
C        Coordinate reference system
C            Y(1) = R      = vertical
C            Y(2) = THETA  = south
C            Y(3) = PHI    = east
C........+.........+.........+.........+.........+.........+.........+..
C
      RY2 = Y(2)
      RY3 = Y(3)
      RY1 = Y(1)
C
      GDAZ = GDAZD/RAD
      GDZE = GDZED/RAD
C     TSGDZE = SIN(GDZE)                                                !Sngl
C     TCGDZE = COS(GDZE)                                                !Sngl
C     TSGDAZ = SIN(GDAZ)                                                !Sngl
C     TCGDAZ = COS(GDAZ)                                                !Sngl
      TSGDZE = DSIN(GDZE)                                               !Dbll
      TCGDZE = DCOS(GDZE)                                               !Dbl
      TSGDAZ = DSIN(GDAZ)                                               !Dbll
      TCGDAZ = DCOS(GDAZ)                                               !Dbl
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Get Y1, Y2, Y3 components in Geodetic coordinates
C         Azimuth is measured clockwise from the north
C         in R, THETA, PHI coordinates, in the THETA-PHI plane
C         The angle is 180 - AZD
C........+.........+.........+.........+.........+.........+.........+..
C
      Y1GD =  TCGDZE
      Y2GD = -TSGDZE*TCGDAZ
      Y3GD =  TSGDZE*TSGDAZ
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ The small angle delta at the point in space between the
C        downward Geodetic   direction and the
C        downward Geocentric direction is given by
C        DELTA = Geocentric co-latitude + Geodetic latitude  - 90 (deg)
C
C        We are looking up
C           The rotation from Geodetic vertical to Geocentric Vertical
C               Is always rotation toward the equator
C
C     \/ Convert from Geodetic to Geocentric Components for Y1, Y2,
C........+.........+.........+.........+.........+.........+.........+..
C
      Y1GC =  Y1GD*TCD+Y2GD*TSD
      Y2GC = -Y1GD*TSD+Y2GD*TCD
      Y3GC =  Y3GD
C
C     WRITE (*,1060) GDZED,GDZE,GDAZD,GDAZ,TSGDZE,TCGDZE,TSGDAZ,TCGDAZ
C     WRITE (*,1060) Y1GD,Y2GD,Y3GD,Y1GC,Y2GC,Y3GC
C1060 FORMAT (' 1050',8F15.5)
C
C........+.........+.........+.........+.........+.........+.........+..
C     ***************************************************
C     Main control of trajectory calculations begins here
C     Trajectories are calculated in subroutine SINGLE
C     ***************************************************
C
C     PC     =  rigidity IN GV
C     INDXPC =  index of rigidity in MV (integer)
C     IRSLT  =  trajectory result
C               IRSLT     +1     allowed
C               IRSLT      0     failed
C               IRSLT     -1     re-entrant 
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 130 NDO = 1, INDO
C
         IF (IERRPT.GE.1) WRITE (16,1070) GDLATD,GLOND,KNT,INDO,NDO,
     *                    IDELPC,INDXPC,DELPC, PC
C
         CALL SINGLTJ (PC,IRSLT,INDXPC,Y1GC,Y2GC,Y3GC)
C
         KNT = KNT+1
         INDXPC = INDXPC-IDELPC
         PC = FLOAT(INDXPC)/1000.0
C
C        +.........+.........+.........+.........+.........+.........+..
C        \/ Check termination conditions
C        +.........+.........+.........+.........+.........+.........+..
C
         IF (PC    .LE. 0.0)    GO TO 140
         IF (TSTEP .GE. FSTEP)  GO TO 150
C
  130 CONTINUE
  140 CONTINUE
 1070 FORMAT (' 1070 ',2F7.2,5I6,2F6.2) 
C
C........+.........+.........+.........+.........+.........+.........+..
C        ************************
C        End of main control loop
C        ************************
C     /\ Go read in next control card
C........+.........+.........+.........+.........+.........+.........+..
C
      GO TO 100
C
C........+.........+.........+.........+.........+.........+.........+..
C     ******************************
C     End of trajectory calculations
C     ******************************
C........+.........+.........+.........+.........+.........+.........+..
C
  150 CONTINUE
C
      WRITE (*, 1120) TSTEP,NTRAJC
      WRITE (8, 1120) TSTEP,NTRAJC
 1120 FORMAT (//' TOTAL NUMBER OF STEPS        ',F15.0///
     *          ' TOTAL NUMBER OF TRAJECTORIES',I15///)
      Write (*,1130)
 1130 format (' End program TJALLMAG')
C
      STOP
C
C........+.........+.........+.........+.........+.........+.........+..
C     Y(1) is R coordinate         Y(2) is THETA coordinate
C     Y(3) is PHI coordinate       Y(4) is V(R)
C     Y(5) is V(THETA)             Y(6) is V(PHI)
C     F(1) is R dot                F(2) is THETA dot
C     F(3) is PHI dot              F(4) is R dot dot
C     F(5) is THETA dot dot        F(6) is PHI dot dot
C     BR   is B(R)                 BT   is B(THETA)
C     BP   is B(PHI)               B    is magnitude of magnetic field
C........+.........+.........+.........+.........+.........+.........+..
C
C     ierrpt vlu  Program Format Variables printed out
C     IERRPT = 1  "MAIN"   1070  Input to SINGLTJ 
C     IERRPT = 1  SINGLTJ  2000  Input to SINGLTJ 
C     IERRPT = 2  SINGLTJ  2070  PC,BETA,KBF,RCKBETA,NSTEP,TBETA,Y,H
C     IERRPT = 4  SINGLTJ  2090  Y,F,ACCER,H,NSTEP
C     IERRPT = 3  SINGLTJ  2100  H,HCK,Y(1),DELACC,PC,NSTEP
C     IERRPT = 3  SINGLTJ  2110  H,HCK,Y(1),RFA,   PC,NSTEP
C     IERRPT = 3  SINGLTJ  2120  H,HCK,Y(1),NAMX,F(ICK),ICK,FOLD(ICK),
C                                ICK,PC,STEP
C     IERRPT = 4  SINGLTJ  2130  Y(1),DISCK,PVEL,H,HSNEK,HOLD,NSTEP
C     IERRPT = 4  SINGLTJ  2140  Y(1),DISCK,PVEL,H,      HOLD,NSTEP
C
      END
      SUBROUTINE GDGC (TCD, TSD)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Convert from Geodetic to Geocentric coordinates
C        Adopted from NASA ALLMAG
C         GDLATD = Geodetic   latitude (in degrees)
C         GCLATD = Geocentric latitude (in degrees)
C         GDCLT  = Geodetic co-latitude
C         ERPLSQ is earth radius AT poles   squared = 40408585 (km sq)
C         EREQSQ is earth radius AT equator squared = 40680925 (km sq)
C         ERADPR is earth polar      radius = 6356.774733 (km)
C         ERADER is earth equatorial radius = 6378.160001 (km)
C         ERAD   is earth average    radius = 6371.25     (km)
C         ERADFL is flattening factor = 1.0/298.25
C         ERADFL =  (ERADEQ - factor)/ERADEQ
C         ERECSQ is eccentricity squared = 0.00673966
C         ERECSQ =  EREQSQ/ERPLSQ - 1.0
C........+.........+.........+.........+.........+.........+.........+..
C
CLast Mod 15 Jan 97  Common block SNGLR & SNGLI
C     Mod    Feb 96  Standard reference TJ1V line check
C
C........+.........+.........+.........+.........+.........+.........+..
C     Programmer - Don F. Smart; FORTRAN77
C     Note - The programming adheres to the conventional FORTRAN
C            default standard that variables beginning with
C            'i','j','k','l','m',or 'n' are integer variables
C            Variables beginning with "c" are character variables
C            All other variables are real
C........+........+.........+.........+.........+.........+.........+..
C            Do not mix different type variables in same common block
C            Some computers do not allow this
C........+.........+.........+.........+.........+.........+.........+..
C
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL * 8(A-B)
      IMPLICIT REAL * 8(D-H)
      IMPLICIT REAL * 8(O-Z)
C
C........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /WRKVLU/ F(6),Y(6),ERAD,EOMC,VEL,BR,BT,BP,B
      COMMON /WRKTSC/ TSY2,TCY2,TSY3,TCY3
      COMMON /TRIG/   PI,RAD,PIO2
      COMMON /GEOID/  ERADPL, ERECSQ
      COMMON /SNGLR/  SALT,DISOUT,GCLATD,GDLATD,GLOND,GDAZD,GDZED,
     *                RY1,RY2,RY3,RHT,TSTEP
C
C........+.........+.........+.........+.........+.........+.........+..
C
      ERPLSQ = 40408585.0
      EREQSQ = 40680925.0
C     ERADPL = SQRT(ERPLSQ)                                             !Sngl
      ERADPL = DSQRT(ERPLSQ)                                            !Dbl
      ERECSQ = EREQSQ/ERPLSQ - 1.0
C
      GDCLT = PIO2-GDLATD/RAD
C     TSGDCLT = SIN(GDCLT)                                              !Sngl
C     TCGDCLT = COS(GDCLT)                                              !Sngl
      TSGDCLT = DSIN(GDCLT)                                             !Dbll
      TCGDCLT = DCOS(GDCLT)                                             !Dbll

      ONE = EREQSQ*TSGDCLT*TSGDCLT
      TWO = ERPLSQ*TCGDCLT*TCGDCLT
      THREE = ONE+TWO
C     RHO = SQRT(THREE)                                                 !Sngl
      RHO = DSQRT(THREE)                                                !Dbll
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Get geocentric distance from geocenter in kilometers
C........+.........+.........+.........+.........+.........+.........+..
C
C     DISTKM = SQRT(SALT*(SALT+2.0*RHO)+(EREQSQ*ONE+ERPLSQ*TWO)/THREE)  !Sngl
      DISTKM = DSQRT(SALT*(SALT+2.0*RHO)+(EREQSQ*ONE+ERPLSQ*TWO)/THREE) !Dbll
C
C........+.........+.........+.........+.........+.........+.........+..
C     TCD and TSD are sine and cosine of the angle the Geodetic vertical
C         must be rotated to form the Geocentric Vertical
C........+.........+.........+.........+.........+.........+.........+..
C
      TCD = (SALT+RHO)/DISTKM
      TSD = (EREQSQ-ERPLSQ)/RHO*TCGDCLT*TSGDCLT/DISTKM
      TCY2 = TCGDCLT*TCD-TSGDCLT*TSD
      TSY2 = TSGDCLT*TCD+TCGDCLT*TSD
C
C     Y(2) = ACOS(TCY2)                                                 !Sngl
      Y(2) = DACOS(TCY2)                                                !Dbll
      Y(3) = GLOND/RAD
      Y(1) = DISTKM/ERAD
C
      GCLATD = (PIO2-Y(2))*RAD
C
C     WRITE (*,1200) GDLATD,GDCLT,TSGDCLT,TCGDCLT,ONE,TWO,THREE,RHO
C1200 FORMAT (' 1200',8F15.5)
C     WRITE (*,1200) DISTKM,TCD,TSD,TCY2,TSY2,GCLATD
C
      RETURN
      END
      SUBROUTINE SINGLTJ (PC,IRSLT,INDXPC,Y1GC,Y2GC,Y3GC)
C
C........+.........+.........+.........+.........+.........+.........+..
C     Cosmic-ray trajectory calculations subroutine
C                calculates cosmic ray trajectory at rigidity  PC
C........+.........+.........+.........+.........+.........+.........+..
C        PC     = rigidity in GV
C        IRSLT  = trajectory result
C        INDXPC = index of rigidity in mv (integer)
C                 Y1GC,Y2GC,Y3GC are initial geocentric coorinates
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Step size optimization & look ahead for potential BETA problems
C             monitor accelerating terms and reduce step length 
C                     if large increase occurs
C        Restart at smaller step size if BETA error occurs
C........+.........+.........+.........+.........+.........+.........+..
C     Restrictions: Cannot run over N or S pole; will get BETA blowup
C........+.........+.........+.........+.........+.........+.........+..
CLast Mod 17 Feb 99  if (ymax.lt.6.6) IFATE = 3
C     Mod 18 Jan 97  Patch high latitude beta problem
C     Mod    Jan 97  High latitude step size adjust, introduce AHLT
C     Mod    Jun 96  EDIF limit set to 1.0e-5
C     Mod    Jun 96  IERRPT formats, Boundary and look ahead 
C     Mod    FEB 96  standard reference TJ1V (line check 17 Feb)

C........+.........+.........+.........+.........+.........+.........+..
C     Programmer - Don F. Smart; FORTRAN77
C     Note - The programming adheres to the conventional FORTRAN
C            default standard that variables beginning with
C            'i','j','k','l','m',or 'n' are integer variables
C            Variables beginning with "c" are character variables
C            All other variables are real
C........+........+.........+.........+.........+.........+.........+..
C            Do not mix different type variables in same common block
C            Some computers do not allow this
C........+.........+.........+.........+.........+.........+.........+..
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL * 8(A-B)
      IMPLICIT REAL * 8(D-H)
      IMPLICIT REAL * 8(O-Z)
C
C........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /WRKVLU/ F(6),Y(6),ERAD,EOMC,VEL,BR,BT,BP,B
      COMMON /WRKTSC/ TSY2,TCY2,TSY3,TCY3
      COMMON /TRIG/   PI,RAD,PIO2
      COMMON /GEOID/  ERADPL, ERECSQ
      COMMON /SNGLR/  SALT,DISOUT,GCLATD,GDLATD,GLOND,GDAZD,GDZED,
     *                RY1,RY2,RY3,RHT,TSTEP
      COMMON /SNGLI/  LIMIT,NTRAJC,IERRPT
C
C........+.........+.........+.........+.........+.........+.........+..
C
      DIMENSION P(6),Q(6),R(6),S(6),YB(6),FOLD(6),YOLD(6)
C
C........+.........+.........+.........+.........+.........+.........+..
C
      CHARACTER*1 CF,CR
      CHARACTER*6 CNAME
C
      DATA CF,CR / 'F','R'/
      DATA CNAME / '  I95 '/                                            ###
C
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (IERRPT.GT.0) WRITE (16,2000) PC,INDXPC,RY1,RY2,RY3
 2000 FORMAT (' SINGLTJ ',F8.3,I8,3F8.4)
C
      BETAST = 2.0
      LSTEP = 0
      KBF = 0
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Runge-Kutta constants
C........+.........+.........+.........+.........+.........+.........+..
C
      RC1O6 = 1.0D0/6.0D0                                               !Dbl
C     SR2 = SQRT(2.0)                                                   !Sngl
      SR2 = DSQRT(2.0D0)                                                !Dbl
      TMS2O2 = (2.0-SR2)/2.0
      TPS2O2 = (2.0+SR2)/2.0
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Initialize Runge-Kutta variables to zero
C........+.........+.........+.........+.........+.........+.........+..
C
  100 DO 110 I = 1, 6
         YB(I) = 0.0
         S(I) = 0.0
         R(I) = 0.0
         Q(I) = 0.0
         P(I) = 0.0
         F(I) = 0.0
  110 CONTINUE
C
      NMAX = 0
      NMIN = 0
      NSTEP = 0
      NSTEPT = 0
C
      TAU = 0.0
      TU100 = 0.0
      YMAX = RY1
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Define initial point at start of trajectory
C        Y(1), Y(2), Y(3) are the position vectors
C........+.........+.........+.........+.........+.........+.........+..
C
      Y(1) = RY1
      Y(2) = RY2
      Y(3) = RY3
C     GRNDKM = (ERADPL/SQRT(1.0-ERECSQ*TSY2SQ))                         !Sngl
      GRNDKM = (ERADPL/DSQRT(1.0-ERECSQ*TSY2SQ))                        !Dbl
      Y10 = (RHT+GRNDKM)/ERAD
      R120KM = (ERAD+120.0)/ERAD
C
C........+.........+.........+.........+.........+.........+.........+..
C     Rigidity = momentum/charge
C                use oxygen 16 as reference isotope
C     Constants used from Handbook of Physics (7-170)
C                         1 amu = 0.931141  GeV
C........+.........+.........+.........+.........+.........+.........+..
C
      ANUC = 16.0
      ZCHARGE = 8.0
C
      EMCSQ = 0.931141
C     TENG = SQRT((PC*ZCHARGE)**2+(ANUC*EMCSQ)**2)                      !Sngl
      TENG = DSQRT((PC*ZCHARGE)**2+(ANUC*EMCSQ)**2)                     !Dbl
      EOMC = -8987.566297*ZCHARGE/TENG
C     GMA  = SQRT(((PC*ZCHARGE)/(EMCSQ*ANUC))**2+1.0)                   !Sngl
C     BETA = SQRT(1.0-1.0/(GMA*GMA))                                    !Sngl
      GMA  = DSQRT(((PC*ZCHARGE)/(EMCSQ*ANUC))**2 + 1.0D0)              !Dbl
      BETA = DSQRT(1.0D0 - 1.0D0/(GMA*GMA))                             !Dbl
      PVEL = VEL*BETA
      HMAX = 1.0/PVEL
      DISCK = DISOUT - 1.1*HMAX*PVEL
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Set max step length ("HMAX") to 1 earth radii
C        PVEL  is particle velocity in earth radii per second
C        DISCK is check for approaching termination boundary
C                 (within 1.1 steps)
C........+.........+.........+.........+.........+.........+.........+..
C
      EDIF = BETA*1.0E-4
      IF (EDIF.LT.1.0-5)  EDIF = 1.0E-5
      IF (BETA.LT.0.1)    EDIF = 1.0E-4
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Y(4), Y(5), Y(6) are the velocity vectors
C........+.........+.........+.........+.........+.........+.........+..
C
      Y(4) = BETA*Y1GC
      Y(5) = BETA*Y2GC
      Y(6) = BETA*Y3GC
C
      AZD = GDAZD
      ZED = GDZED
      IAZ = AZD+0.01
      IZE = ZED+0.01
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Set HSTART to about 1 % of the time to complete one gyro-radius
C                            in a 1 Gauss field
C        H = [(2.0*PI*33.333*PC)/(BETA*C)]/0.01
C            if restart after BETA error, set HCK to small value
C        Introduce   AHLT  to control step size at high lat (beta problem)
C                    HCK   - reduce step size when large acceleration
C                    HOLD  - last step size used
C                    HCNG  - only allow 20% max growth in step size
C                    HSNEK - attempt to approach boundary quickly
C        Problem at z=90 at high lat
C                add zen angle in deg to reduce first step
C........+.........+.........+.........+.........+.........+.........+..
C
      PTCY2 = ABS(TCY2)
      AHLT = (1.0 + PTCY2)**2
      HSTART = 6.0E-6*PC/(BETA*AHLT + ZED*PTCY2)
      IF (HSTART.LT.1.0E-6) HSTART = 1.0E-6
      HOLD = HSTART
      HCK  = HSTART
      HCNG = HSTART
C
C     WRITE (16, 2010)  HMAX,HOLD,HCK,HCNG,Y(4),Y(5),Y(6),PVEL, NSTEP
C2010 FORMAT (' 2010 ',18X, 4F9.6, 3F9.4, F9.4,9X,15X,I6)
C
C........+.........+.........+.........+.........+.........+.........+..
C     Start Runge-Kutta
C      \/\/\/\/\/\/\/\/
C        \/\/\/\/\/\/
C          \/\/\/\/
C            \/\/
C             \/
C........+.........+.........+.........+.........+.........+.........+..
C     Change in step size criteria, Aug 97
C            remove cos VxB step size, causes problems in tight loops
C            step size is now only a function of B and BETA
C........+.........+.........+.........+.........+.........+.........+..
C
  130 IF (HCK.LT.1.0E-6) HCK = 1.0E-6
      CALL FGRADA
      HB = 1.6E-5*PC/(B*BETA)
      H = HB/BETAST
C
      IF (KBF.GT.0)  H=H/(FLOAT(KBF*2))
      IF (H.GT.HMAX) H = HMAX
      IF (H.GT.HCNG) H = HCNG
      IF (H.GT.HCK)  H = HCK
C
      DO 140 I = 1, 6
         S(I) = H*F(I)
         P(I) = 0.5*S(I)-Q(I)
         YB(I) = Y(I)
         Y(I) = Y(I)+P(I)
         R(I) = Y(I)-YB(I)
         Q(I) = Q(I)+3.0*R(I)-0.5*S(I)
  140 CONTINUE
C
      CALL FGRADA
C
      DO 150 I = 1, 6
         S(I) = H*F(I)
         P(I) = TMS2O2*(S(I)-Q(I))
         YB(I) = Y(I)
         Y(I) = Y(I)+P(I)
         R(I) = Y(I)-YB(I)
         Q(I) = Q(I)+3.0*R(I)-TMS2O2*S(I)
  150 CONTINUE
C
      CALL FGRADA
C
      DO 160 I = 1, 6
         S(I) = H*F(I)
         P(I) = TPS2O2*(S(I)-Q(I))
         YB(I) = Y(I)
         Y(I) = Y(I)+P(I)
         R(I) = Y(I)-YB(I)
         Q(I) = Q(I)+3.0*R(I)-TPS2O2*S(I)
  160 CONTINUE
C
      CALL FGRADA
C
      DO 170 I = 1, 6
         S(I) = H*F(I)
         P(I) = RC1O6*(S(I)-2.0*Q(I))
         YB(I) = Y(I)
         Y(I) = Y(I)+P(I)
         R(I) = Y(I)-YB(I)
         Q(I) = Q(I)+3.0*R(I)-0.5*S(I)
  170 CONTINUE
C
C........+.........+.........+.........+.........+.........+.........+..
C             /\
C            /\/\
C          /\/\/\/\
C        /\/\/\/\/\/\
C      /\/\/\/\/\/\/\/\
c      One Runge-Kutta
c      step completed
C........+.........+.........+.........+.........+.........+.........+..
C
      NSTEP = NSTEP+1
      NSTEPT = NSTEPT + 1
      TAU = TAU+H
      HOLD = H
      HCNG = H*1.2
      HCK  = HCNG
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Emergency diagnostic printout if desired
C........+.........+.........+.........+.........+.........+.........+..
C     WRITE (16, 2030)     H,  Y(1),Y(2),Y(3),   PVEL,B, NSTEP
C     WRITE (16, 2040)  HB,H,HMAX,HOLD,HCK,HCNG,Y(4),Y(5),Y(6),
C    *                  PVEL,B,NSTEP
C2030 FORMAT (' 2030 ',  9X,    F9.6, 36X,    3F9.5,F9.4,F9.5,18X,I6)
C2040 FORMAT (' 2040 ',        6F9.6,         3F9.5,F9.4,F9.5,18X,I6)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Check for altitude less than 100 km
C        if less than 120 km, compute exact altitude above oblate earth
C           and sum time trajectory is below 100 km altitude.
C        set re-entrant altitude at  RHT  km above oblate earth
C               computed from international reference ellipsoid
C     ...+.........+.........+.........+.........+.........+.........+..
C
      IF (Y(1).LT.R120KM) THEN
C        TSY2SQ = SIN(Y(2))**2                                          !Sngl
         TSY2SQ = DSIN(Y(2))**2                                         !Dbl
C        GRNDKM = (ERADPL/SQRT(1.0-ERECSQ*TSY2SQ))                      !Sngl
         GRNDKM = (ERADPL/DSQRT(1.0-ERECSQ*TSY2SQ))                     !Dbl
         R100KM = (100.0+GRNDKM)/ERAD
         R120KM = (120.0+GRNDKM)/ERAD
         IF (Y(1).LT.R100KM) TU100 = TU100+H
         PSALT = Y(1)*ERAD-GRNDKM
         Y10 = (RHT+GRNDKM)/ERAD
C
         IF (NSTEP.GT.5)  THEN
            IF (Y(1).LT.Y10.OR.PSALT.LE.0.0) THEN
               IF (IERRPT.GT.2)  WRITE (16, 2045) PSALT, Y(1), Y10 
               IRT = -1
               GO TO 260
            ENDIF
         ENDIF
      ENDIF
 2045 FORMAT (' 2045 PSALT,Y(1),Y10',F10.6,1PE14.6,E14.6)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Begin error checks
C        (1) Check for unacceptable changes in BETA
C........+.........+.........+.........+.........+.........+.........+..
C
C     RCKBETA = SQRT(Y(4)*Y(4)+Y(5)*Y(5)+Y(6)*Y(6))                     !Sngl
      RCKBETA = DSQRT(Y(4)*Y(4)+Y(5)*Y(5)+Y(6)*Y(6))                    !Dbl
      TBETA = BETA-RCKBETA
      IF (ABS(TBETA).GT.EDIF) THEN
         KBF = KBF+1
         BETAST = BETAST + AHLT
         EDIF = 2.0*EDIF
         IF (RCKBETA.GT.(1.0+EDIF)) THEN
             BETAST = BETAST+FLOAT(KBF)*(1.0+AHLT)
             WRITE (*,2050) KBF,BETA,RCKBETA
             WRITE (*,2060) Y,H,PC,NSTEP
             WRITE (16,2050) KBF,BETA,RCKBETA
             WRITE (16,2060) Y,H,PC,NSTEP
         ENDIF
C
         WRITE (16,2070)  PC,BETA,KBF,RCKBETA,NSTEP,TBETA,Y,H
         WRITE ( *,2070)  PC,BETA,KBF,RCKBETA,NSTEP,TBETA,Y,H
C
C        +.........+.........+.........+.........+.........+.........+..
C        \/ Check for irrecoverable beta error
C           if KBF > 4, set fate to failed and start next rigidity
C        +.........+.........+.........+.........+.........+.........+..
C
         IF (KBF.lt.5) THEN
            GO TO 100
         ELSE
            IRT = 0
            PATH = -PVEL*TAU
            ISALT = SALT+0.0001
            WRITE (*,2080)
            GO TO 280
         ENDIF
      ENDIF
C
 2050 FORMAT (' 2050 ','KBF=',i5,5x,' error, the velocity of the ',
     *   'particle (BETA) has exceeded the velocity of light'/
     *   '          BETA at start of trajectory was ',F10.7/
     *   '          BETA now equals                 ',F10.7/
     *   '          reduce step size and try again    ')
 2060 FORMAT (' 2060 ','Y,H,PC,NSTEP=',8F12.6,I10)
 2070 FORMAT (' 2070 ','ERROR, Particle BETA changed;',
     *        ' PC = ',F10.6,' GV'/
     *     6x,' Beta at start was ',F10.7,17X,'KBF=',I6/
     *     6x,' Beta new equals   ',F10.7,10X,'step number',I6/
     *     6x,' Beta difference   ',F10.7/
     *     6x,' step length reduced and trajectory recalculated '/,
     *     6x,'Y=',6F12.6,6X,' H=',F15.7)
 2080    FORMAT (' 2080 ','Irrecoverable BETA error problem'/6x,
     *      ' Set fate to Failed & make path length negative'/6x,
     *      ' Terminate this trajectory & continue with program')
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ Continue status checks
C       (2) Compute composite acceleration
C     ...+.........+.........+.........+.........+.........+.........+..
C
C     ACCER = SQRT(F(4)*F(4)+F(5)*F(5)+F(6)*F(6))                       !Sngl
      ACCER = DSQRT(F(4)*F(4)+F(5)*F(5)+F(6)*F(6))                      !Dbl
C
      IF (IERRPT.GT.3) WRITE (16,2090) Y, F, ACCER, H, NSTEP
 2090 FORMAT (' Y,F,A,H ',f7.4,5f7.3,1x,1pe8.1,5e8.1,1x,e8.1, e9.2,I6)
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ Continue status checks, make adjustment latitude dependant
C        (3) Monitor change in composite acceleration
C            If composite acceleration (new-old) change > 5
C            If composite acceleration (new/old) ratio  > 2
C               change step size to a smaller value
C    ....+.........+.........+.........+.........+.........+.........+..
C
      IF (NSTEP.GE.2) THEN
         IF (ACCER.GT.ACCOLD) THEN
            DELACC = ACCER-ACCOLD
            IF (DELACC.GT.5.0)  THEN
               HCK = HCK/(1.0+AHLT)
               IF (IERRPT.GT.2) WRITE (16,2100) 
     *                          H,HCK,y(1),DELACC,PC,NSTEP
               RFA = ACCER/ACCOLD
               IF (RFA.GT.2.0) THEN
                  HCK = HCK/(1.0+AHLT)
                  IF (IERRPT.GT.2) WRITE (16,2110)
     *                             H,HCK,Y(1),RFA,PC,NSTEP
               ENDIF
            ENDIF
         ENDIF
C
 2100    FORMAT (' 2100 ','H-REDUCE',2x,'H=',F8.6,2x,'HCK=',F8.6,2x,
     *    'Y(1)=',f7.4,2X,'DELACC=',F6.2,4X,'PC=',F8.3,4x,'NSTEP=',I8)
 2110    FORMAT (' 2110 ','H-REDUCE',2x,'H=',F8.6,2x,'HCK=',F8.6,2x,
     *    'Y(1)=',f7.4,4X,' RFA=',F6.2,4X,'PC=',F8.3,4x,'NSTEP=',I8)
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ Continue status checks, make adjustment latitude dependant
C        (4) Monitor change in acceleration components
C            If change in any acceleration component is more than
C               a factor of 3, reduce step length
C     ...+.........+.........+.........+.........+.........+.........+..
C
         DO 200 ICK = 4, 6
            AFOLD = ABS(FOLD(ICK))
            IF (AFOLD.GT.3.0) THEN
               RFCK = ABS(F(ICK)/AFOLD)
               IF (RFCK.GT.3.0) THEN
                  HCK = HCK/(1.0+AHLT)
                  IF (IERRPT.GT.2)  THEN
                      WRITE (16,2120)  H,HCK,Y(1),NMAX,ICK,F(ICK),
     &                       ICK,FOLD(ICK),PC,NSTEP
                  ENDIF
               ENDIF
            ENDIF
  200    CONTINUE
      ENDIF
C
 2120 FORMAT (' 2120 ','H-reduce',2X,'H=',F8.6,2X,'HCK=',F8.6,2X,
     *        'Y(1)=',F7.4,2X,'NAMX=',I4,2X,'F(',I1,')=',F6.2,2X, 
     *        'FOLD(',I1,')=',F6.2,2X,'PC=',F6.3,2X,'NSTEP=',I6)
C
      ACCOLD = ACCER
C
C........+.........+.........+.........+.........+.........+.........+..
C     /\ Error checks complete
C
C     \/ Find if a max or a min has occurred
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (NSTEP.GT.1) THEN
         IF (YOLD(4).LE.0.0.AND.Y(4).GT.0.0) NMIN = NMIN+1
         IF (YOLD(4).GE.0.0.AND.Y(4).LT.0.0) NMAX = NMAX+1
      ENDIF
C
      IF (Y(1).GT.YMAX) YMAX = Y(1)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Check for termination conditions
C         Allowed    - radial distance exceeded disout
C         Failed     - number of steps exceeded
C         Re-entrant - trajectory is below "top" of atmosphere
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ (1) Check for step limit exceeded
C     ...+.........+.........+.........+.........+.........+.........+..
C
      IF (NSTEP.GE.LIMIT) THEN
         IRT = 0
         GO TO 260
      ENDIF
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/  (2) Check if y(1) within 1.1 max step lengths of disout.
C                   if so, reduce step size and
C                          approach boundary at smaller step
C     ...+.........+.........+.........+.........+.........+.........+..
C
      IF (Y(1).GT.DISCK) THEN
         DISTR = ABS(DISOUT - Y(1))
         HSNEK = DISTR/PVEL
         HCNG = HCNG/2.0
         HCK  = HCK/2.0
         IF (HSNEK .LT. HCNG)  HCNG = HSNEK 
         IF (HSNEK .LT. HCK)   HCK  = HSNEK 
         DISCK = DISOUT - DISTR/2.0
         IF (DISCK.GE.DISOUT)  THEN
            DISCK = 24.999
             GO TO  210
         ENDIF
         IF (H.LT.1.0E-5 .OR. HCK.LT.1.0E-5 .OR. HCNG.LT.1.0E-5)  THEN
            H    = 1.0E-5
            HCK  = 1.0E-5
            HCNG = 1.0E-5
         ENDIF
C
         IF (IERRPT.GT.3) WRITE (16,2130) Y(1),DISCK,PVEL,H,HSNEK,NSTEP
C
  210    IF (Y(1).GT.DISOUT) THEN
            IF (H.LE.1.0E-5)  THEN
               IRT = 1
               GO TO 260
            ENDIF
            TAU = TAU - H
C
C           .......+.........+.........+.........+.........+.........+..
C           \/ Backup option invoked if you are here
C           .......+.........+.........+.........+.........+.........+..
C
            DO 220 I = 1, 6
               Y(I) = YOLD(I)
               F(I) = FOLD(I)
  220       CONTINUE
         ENDIF
         GO TO 130
      ENDIF
 2130 FORMAT (' 2130 ',2X,'Y(1),DISCK,PVEL,H,HSNEK', 
     *        4X,1PE12.6,4X,E12.6,4X,E12.6,4X,2E9.2,22X,I6)
C
C     +.........+.........+.........+.........+.........+.........+..
c     \/ Have penetrated boundary if you are here.
c          if large step size, go back one step and
c             reduce step length (and adjust "TAU")
C     +.........+.........+.........+.........+.........+.........+..
C
  230    IF (Y(1).GT.DISOUT) THEN
C
         IF (IERRPT.GT.3)  WRITE (16, 2140) Y(1),DISCK,PVEL,H,NSTEP
C
         IF (H.LT.1.0E-5 .OR. HCK.LT.1.0E-5 .OR. HCNG.LT.1.0E-5)  THEN
            IRT = 1
            GO TO 260
         ELSE
            HCK  = HCK/2.0
            HCNG = HCNG/2.0
            TAU = TAU - H
            DO 240 I = 1, 6
               Y(I) = YOLD(I)
               F(I) = FOLD(I)
  240       CONTINUE
            GO TO 130
         ENDIF
      ENDIF
C
 2140 FORMAT (' 2140 ',2X,'Y(1),DISCK,PVEL,H', 
     *        4X,1PE12.6,4X,E12.6,4X,E12.6,4X,E9.2,27X,I6)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ STORE VALUES OF  Y  AND  F  AS FOLD & YOLD
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 250 I = 1, 6
         YOLD(I) = Y(I)
         FOLD(I) = F(I)
  250 CONTINUE
C
      GO TO 130
C
C........+.........+.........+.........+.........+.........+.........+..
C*********************************************************************
C        **********          **********          **********
C                  **********          **********
C                           **********
C                 TRAJECTORY COMPLETE IF YOU ARE HERE
C........+.........+.........+.........+.........+.........+.........+..
C
  260 CONTINUE
C
      IF (Y(1).GE.DISOUT)  IRT = 1
      PATH  = PVEL*TAU
      ISALT = SALT+0.0001
      LSTEP = BETAST - 1.9
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ WRITE OUT RESULTS
C        IRT     +1     ALLOWED          (FATE = 0)
C        IRT      0     FAILED           (FATE = 2)
C        IRT     -1     RE-ENTRANT       (FATE = 1)
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (IRT.GT.0) THEN
C        TCY2 = COS(Y(2))                                               !Sngl
C        TSY2 = SIN(Y(2))                                               !Sngl
         TCY2 = DCOS(Y(2))                                              !Dbll
         TSY2 = DSIN(Y(2))                                              !Dbll
         YDA5 = Y(5)*TCY2+Y(4)*TSY2
         ATRG1 = Y(4)*TCY2-Y(5)*TSY2
C        ATRG2 = SQRT(Y(6)*Y(6)+YDA5*YDA5)                              !Sngl
         ATRG2 = DSQRT(Y(6)*Y(6)+YDA5*YDA5)                             !Dbll
         FASLAT = 0.0
C        IF (ATRG1.NE.0.0.AND.ATRG2.NE.0.0)                             !Sngl
C    *       FASLAT = ATAN2(ATRG1,ATRG2)*RAD                            !Sngl
         IF (ATRG1.NE.0.0.AND.ATRG2.NE.0.0)                             !Dbll
     *       FASLAT = DATAN2(ATRG1,ATRG2)*RAD                           !Dbll
         FASLON = Y(3)*RAD
C        IF (Y(6).NE.0.0.AND.YDA5.NE.0.0)                               !Sngl
C    *       FASLON = (Y(3)+ ATAN2(Y(6),YDA5))*RAD                      !Sngl
         IF (Y(6).NE.0.0.AND.YDA5.NE.0.0)                               !Dbl
     *       FASLON = (Y(3)+ DATAN2(Y(6),YDA5))*RAD                     !Dbl
         IF (FASLON.LT.0.0)    FASLON = FASLON+360.0
         IF (FASLON.GT.360.0)  FASLON = FASLON-360.0
C
         WRITE (8,2150) GDLATD,GCLATD,GLOND,IZE,IAZ,PC,FASLAT,FASLON,
     *                  PATH,NMAX,NSTEP,TU100,YMAX,LSTEP,SALT,CNAME
C
         IFATE = 0
         WRITE (7,2160) GDLATD,GLOND,PC,ZED,AZD,ISALT,FASLAT,FASLON,
     *                  NSTEP,IFATE,CNAME
      ENDIF
 2150 FORMAT (2F7.2,F9.2,I5,I4,F10.3,2F8.2,F11.5,I4,I7,F9.5,F9.4,
     *        I4,F11.1,1X,A6,13X)
 2160 FORMAT (F7.2,F8.2,F9.3,2F6.1,I7,F7.2,F8.2,I7,3X,I3,3X,A6)
C
      IF (IRT.LT.0) THEN
         RENLAT = (PIO2-Y(2))*RAD
         RENLON = Y(3)*RAD
C
         WRITE (8,2170) GDLATD,GCLATD,GLOND,IZE,IAZ,PC,CR,CR,PATH,NMAX
     *      ,NSTEP,TU100,YMAX,LSTEP,SALT,CNAME,RENLAT,RENLON
C
         IFATE = 1
         WRITE (7,2180) GDLATD,GLOND,PC,ZED,AZD,ISALT,NSTEP,IFATE,CNAME
      ENDIF
 2170 FORMAT (2F7.2,F9.2,I5,I4,F10.3,5X,A1,2X,5X,A1,2X,F11.5,I4,I7,
     *        F9.5,F9.4,I4,F11.1,1X,A6,F6.1,F7.1)
 2180 FORMAT (F7.2,F8.2,F9.3,2F6.1,I7,4X,'R',7X,'R',I9,3X,I3,3X,A6)
C
  280 IF (IRT.EQ.0) THEN
C
         WRITE (8,2190) GDLATD,GCLATD,GLOND,IZE,IAZ,PC,CF,CF,PATH,
     *                  NMAX,NSTEP,TU100,YMAX,LSTEP,SALT,CNAME
C
         IFATE = 2
         IF (YMAX.LT.6.6) IFATE = 3
         WRITE (7,2200) GDLATD,GLOND,PC,ZED,AZD,ISALT,NSTEP,IFATE,CNAME
      ENDIF
 2190 FORMAT (2F7.2,F9.2,I5,I4,F10.3,5X,A1,2X,5X,A1,2X,F11.5,I4,I7,
     *        F9.5,F9.5,I4,F11.1,1X,A6,13X)
 2200 FORMAT (F7.2,F8.2,F9.3,2F6.1,I7,4X,'F',7X,'F',I9,3X,I3,3X,A6)
C
      NTRAJC = NTRAJC+1
      TSTEP = TSTEP+FLOAT(NSTEP)
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ Comment out to reduce IO
C     ...+.........+.........+.........+.........+.........+.........+..
C
C     WRITE (*,2210) PC, ZED, AZD, NSTEP, IFATE
C2210 FORMAT (1H+, 22X, 3F7.2,7x,2I6)
C
      IRSLT = IRT
      RETURN
      END
      SUBROUTINE FGRADA
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Version of FGRAD to run with NASA NSSDC ALLMAG subroutine   
C........+.........+.........+.........+.........+.........+.........+..
CLast Mod 22 Dec 00  Use NASA ALLMAG for magnetic field calculations
C     Mod    Feb 96  standard reference TJ1V (line check 17 Feb)
c........+.........+.........+.........+.........+.........+.........+..
c     Programmer   -  Don F. Smart; FORTRAN77
c     Note -  The programming adheres to the conventional FORTRAN
c             default standard that variables beginning with
c            'i','j','k','l','m',or 'n' are integer variables
c             Variables beginning with "c" are character variables
c             All other variables are real
c........+.........+.........+.........+.........+.........+.........+..
c        Do not mix different type variables in same common block
c           Some computers do not allow this
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL * 8(A-B)
      IMPLICIT REAL * 8(D-H)
      IMPLICIT REAL * 8(O-Z)
c
c........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /WRKVLU/ F(6),Y(6),ERAD,EOMC,VEL,BR,BT,BP,B
      COMMON /WRKTSC/ TSY2,TCY2,TSY3,TCY3
C
c........+.........+.........+.........+.........+.........+.........+..
C
      F(1) = VEL*Y(4)
      F(2) = VEL*Y(5)/Y(1)
C     TSY2 = SIN(Y(2))                                                  !Sngl
C     TCY2 = COS(Y(2                                                    !Sngl
      TSY2 = DSIN(Y(2))                                                 !Dbl
      TCY2 = DCOS(Y(2))                                                 !Dbl
      F(3) = VEL*Y(6)/(Y(1)*TSY2)
      SQY6 = Y(6)*Y(6)/Y(1)
      Y5OY1 = Y(5)/Y(1)
      TAY2 = TSY2/TCY2

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Use NSSDC routine ALLMAG for magnetic field calculations
c            define MODEL and epoch year (TM)
c            need sine and cosine of phi (Y(3) for ALLMAG
C            need radial distance from earth center in kilometers (RKM)
c        Remember, ALMAG returns magnetic field in units of Gauss
c........+.........+.........+.........+.........+.........+.........+..

      MODEL = 14
      TM = 1995.0
C     TSY3 = SIN(Y(3))                                                  !Sngl
C     TCY3 = COS(Y(3))                                                  !Sngl
      TSY3 = DSIN(Y(3))                                                 !Dbl
      TCY3 = DCOS(Y(3))                                                 !Dbl
      RKM = Y(1)*6371.2

      CALL ALLMAG1 (MODEL,TM,RKM,TSY2,TCY2,TSY3,TCY3,BR,BT,BP,B)

c     write (16,1818)  MODEL,TM, Y(1),Y(3),Y(3), RKM,  
c    *                 TSY2,TCY2,TSY3,TCY3, BR,BT,BP,B
c1818 format( i5, f8.1, 2x, 3F10.5,2X, f12.3, 2x, 4f10.5,2x, 4f10.5)

      F(4) = EOMC*(Y(5)*BP-Y(6)*BT)+VEL*(Y(5)*Y5OY1+SQY6)
      F(5) = EOMC*(Y(6)*BR-Y(4)*BP)+VEL*(SQY6/TAY2-Y5OY1*Y(4))
      F(6) = EOMC*(Y(4)*BT-Y(5)*BR)-VEL*((Y5OY1*Y(6))/TAY2+Y(4)*Y(6)/
     *       Y(1))
      RETURN
C
C........+.........+.........+.........+.........+.........+.........+..
C     Y(1) is R coordinate         Y(2) is THETA coordinate
C     Y(3) is PHI coordinate       Y(4) is V(R)
C     Y(5) is V(THETA)             Y(6) is V(PHI)
C     F(1) is R dot                F(2) is THETA dot
C     F(3) is PHI dot              F(4) is R dot dot
C     F(5) is THETA dot dot        F(6) is PHI dot dot
C     BR   is B(R)                 BT   is B(THETA)
C     BP   is B(PHI)               B    is magnitude of magnetic field
C........+.........+.........+.........+.........+.........+.........+..
C
      END
C  Program: allmag_sub.f     Version: 1.1     Last Updated: 12/30/97 13:24:51
C     This subroutine contains all the geomagnetic field coefficients.
C                         Written in Digital UNIX FORTRAN      (12/01/1997)  
!***************************************************************************! 
!                           SUBROUTINE ALLMAG1                              !
!***************************************************************************!
      SUBROUTINE ALLMAG1 (MODEL,TM,RKM,ST,CT,SPH,CPH,BR,BT,BP,B)        ALMGL001

C     NOTE ADDITION OF NEXT STATEMENT
C ****  GEOCENTRIC VERSION OF GEOMAGNETIC FIELD ROUTINE                 ALMGL002
C ****  LONG DECK, THROUGH NMAX=13, FIXED INDICES WITHOUT DO LOOPS      ALMGL004
C ****  EXECUTION TIME PER CALL FACTOR OF THREE LESS THAN SHORT DECK    ALMGL005
C ****  PROGRAM DESIGNED AND TESTED BY E G STASSINOPOULOS AND G D MEAD, ALMGL006
C ****  CODE 641, NASA GODDARD SPACE FLT CTR, GREENBELT, MD 20771       ALMGL007
C *****   INPUT  MODEL    CHOICE OF 14 MODELS - SEE BELOW               ALMGL008
C *****          RKM      GEOCENTRIC DISTANCE IN KILOMETERS             ALMGL009
C *****          TM       TIME IN YEARS FOR DESIRED FIELD               ALMGL010
C *****          ST,CT    SIN + COS OF GEOCENTRIC COLATITUDE            ALMGL011
C *****          SPH,CPH  SIN + COS OF EAST LONGITUDE                   ALMGL012
C *****  OUTPUT  BR,BT,BP GEOCENTRIC FIELD COMPONENTS IN GAUSS          ALMGL013
C *****          B        FIELD MAGNITUDE IN GAUSS                      ALMGL014
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /TRAJAC/ CONSTEM,T,FILENAM
      COMMON /DIPOLE/ WLONG,COLAT,EM

      DIMENSION T0(14),NMX(14),ISUM(14,3),G(13,13)                      ALMGL024

      DATA T0 /4*1960.,2*1965.,1970.,1980.,3*1975.,1985.,1990.,1995./,  ALMGL025
     $ NMX /10,11,12,11,9,9,13,11,9,13,13,13,13,13/

      INTEGER LSUM(14,3)/-1646106,-1795169,-1865298,-1777057,-158472,   ALMGL026
     A-156856,-2191704,-1996220,-168051,-2252599,-2445733,-2369772,     ALMGL027
     B-2409473,-246795,-62661,-96778,-181519,-83555,-9569,-9599,-8593,
     C-5412,-7351,-11947,-10278,-6777,-8938,-12380,1,-10618,7*1,-2698,
     D1,1,1,1/

      INTEGER*4 G1(13,13),GT1(13,13),GTT1(13,13),G2(13,13),GT2(13,13),  ALMGL029
     1 GTT2(13,13),G3(13,13),GT3(13,13),GTT3(13,13),G4(13,13),          ALMGL030
     2 GT4(13,13),GTT4(13,13),G5(13,13),GT5(13,13),GTT5(13,13),         ALMGL031
     3 G6(13,13),GT6(13,13),GTT6(13,13),G7(13,13),GT7(13,13),GTT7(13,13)ALMGL032
     4 ,G8(13,13),GT8(13,13),GTT8(13,13),
     5 G9(13,13),GT9(13,13),GTT9(13,13),
     6 G10(13,13),GT10(13,13),GTT10(13,13),
     7 G11(13,13),GT11(13,13),GTT11(13,13),
     8 G12(13,13),GT12(13,13),GTT12(13,13),
     8 G13(13,13),GT13(13,13),GTT13(13,13),
     8 G14(13,13),GT14(13,13),GTT14(13,13)
     9 ,LG(13,13,14),LGT(13,13,14),LGTT(13,13,14)                       ALMGL033

      REAL*4  GG(13,13,14),GGT(13,13,14),GGTT(13,13,14),SHMIT(13,13)    ALMGL034
      EQUIVALENCE (G1(1),GG(1),LG(1)), (GT1(1),GGT(1),LGT(1)),          ALMGL035
     A     (GTT1(1),GGTT(1),LGTT(1)),                                   ALMGL036
     B (G2(1),LG(1,1,2)), (GT2(1),LGT(1,1,2)), (GTT2(1),LGTT(1,1,2)),   ALMGL037
     C (G3(1),LG(1,1,3)), (GT3(1),LGT(1,1,3)), (GTT3(1),LGTT(1,1,3)),   ALMGL038
     D (G4(1),LG(1,1,4)), (GT4(1),LGT(1,1,4)), (GTT4(1),LGTT(1,1,4)),   ALMGL039
     E (G5(1),LG(1,1,5)), (GT5(1),LGT(1,1,5)), (GTT5(1),LGTT(1,1,5)),   ALMGL040
     F (G6(1),LG(1,1,6)), (GT6(1),LGT(1,1,6)), (GTT6(1),LGTT(1,1,6)),   ALMGL041
     G (G7(1),LG(1,1,7)), (GT7(1),LGT(1,1,7)), (GTT7(1),LGTT(1,1,7)),   ALMGL042
     H (G8(1),LG(1,1,8)), (GT8(1),LGT(1,1,8)), (GTT8(1),LGTT(1,1,8)),
     I (G9(1),LG(1,1,9)), (GT9(1),LGT(1,1,9)), (GTT9(1),LGTT(1,1,9)),
     J (G10(1),LG(1,1,10)),(GT10(1),LGT(1,1,10)),(GTT10(1),LGTT(1,1,10))
     K,(G11(1),LG(1,1,11)),(GT11(1),LGT(1,1,11)),(GTT11(1),LGTT(1,1,11))
     L,(G12(1),LG(1,1,12)),(GT12(1),LGT(1,1,12)),(GTT12(1),LGTT(1,1,12))
     M,(G13(1),LG(1,1,13)),(GT13(1),LGT(1,1,13)),(GTT13(1),LGTT(1,1,13))
     N,(G14(1),LG(1,1,14)),(GT14(1),LGT(1,1,14)),(GTT14(1),LGTT(1,1,14))

C  ***** THE FOLLOWING DATA CARDS CONTAIN THE FIELD COEFFICIENTS        ALMGL043
C  ***** FOR THE FOLLOWING SEVEN MODELS                                 ALMGL044
C  *****  G1,GT1  HENDRICKS + CAIN      99-TERM  GSFC  9/65  EPOCH 1960.ALMGL045
C  *****  G2,GT2,GTT2  CAIN ET. AL.    120-TERM  GSFC 12/66  EPOCH 1960.ALMGL046
C  *****  G3,GT3  CAIN + LANGEL        143-TERM  POGO 10/68  EPOCH 1960.ALMGL047
C  *****  G4,GT4  CAIN + SWEENEY       120-TERM  POGO  8/69  EPOCH 1960.ALMGL048
C  *****  G5,GT5  IGRF 1965.0           80-TERM       10/68  EPOCH 1965.ALMGL049
C  *****  G6,GT6  LEATON MALIN + EVANS 1965        80-TERM   EPOCH 1965.ALMGL050
C  *****          FOR MODEL 6 (LME 1965) SET RKM = 6371.2 + ALTITUDE    ALMGL051
C  *****  G7,GT7  HURWITZ US COAST + GEODETIC S.  168-TERM   EPOCH 1970.ALMGL052
C  *****   8      IGRF 1980                                  EPOC980. 
C  *****   9      IGRF 1975             80-TERM80-TERM    EPOCH 1975.
C  *****  10      BARRACLOUGH          168-TERM                   "
C  *****  11      AWC                     "                       "             
C  *****  12      IGRF 1985               "                     1985            
C  *****  13      IGRF 1990               "                     1990
C  *****  14      IGRF 1995               "                     1995
C                                                         
C     HENDRICKS / CAIN ET AL MODEL * 99-TERM GSFC 9/65 EPOCH 1960.0             
                                                                                
      DATA G1 / 10, -304249,-15361,13009,9576,-2277,498,709,48,99,3*0,  ALMGL053
     A 57748,-21616,30002,-19870,8028,3595,607,-572,67,29,3*0,-19498,   ALMGL054
     B 2043,15853,12904,5026,2313,45,56,-88,74,3*0,-4310,2308,-1300,8712ALMGL055
     C ,-3940,-312,-2417,75,-138,-156,3*0,1520,-2684,29,-2505,2714,     ALMGL056
     D -1573,-12,-244,-33,114,3*0,86,1212,-1160,-1104,799,-652,5,-15,71,ALMGL057
     E 111,3*0,-119,1028,609,-272,-124,-116,-1091,141,-56,10,3*0,-540,  ALMGL058
     F -244,-91,22,276,-211,-201,58,117,4*0,69,-122,58,-170,26,236,-25, ALMGL059
     G -160,64,16,3*0,-220,156,51,-35,-18,96,121,2,-25,15,42*0 /        ALMGL060
      DATA GT1 / 100, 2059,-2907,266,-86,255,-70,6*0,-394,602,121,-1003,ALMGL061
     H 194,-8,99,6*0,-1369,-1578,-70,163,-117,153,85,6*0,649,293,-924,  ALMGL062
     I -130,-54,-42,211,6*0,-177,-154,318,-548,-417,-72,157,6*0,304,288,ALMGL063
     J -186,125,80,164,-9,6*0,-139,12,153,-73,-6,45,6,84*0/             ALMGL064
      DATA GTT1 /1,168*0/                                               ALMGL065
                                                                                
C     CAIN ET AL 120-TERM GSFC 12/66 EPOCH 1960.0                               
                                                                                
      DATA G2 / 10, -304012,-15401,13071,9493,-2335,492,722,85,104,-29, ALMGL066
     A 2*0,57782,-21638,29979,-19889,8035,3557,575,-537,65,58,-9,       ALMGL067
     B 2*0,-19320,2029,15903,12768,5029,2284,-8,79,-93,75,-22,2*0,-4254,ALMGL068
     C 2278,-1338,8812,-3977,-288,-2383,156,-96,-151,8,2*0,1603,-2743,  ALMGL069
     D 23,-2466,2665,-1579,-15,-243,-61,121,-28,2*0,51,1178,-1148,-1089,ALMGL070
     E 824,-622,-20,-36,55,47,64,2*0,-121,1044,566,-234,-148,-133,-1089,ALMGL071
     F 155,-81,2,47,2*0,-537,-274,-81,70,243,-225,-214,36,130,16,-2,2*0,ALMGL072
     G 54,-117,42,-153,46,219,-7,-171,74,9,18,2*0,-224,138,63,-30,-19,  ALMGL073
     H 90,115,1,-15,2,20,2*0,-1,45,-10,26,-44,-13,-36,40,10,-20,11,28*0/ALMGL074
      DATA GT2 / 100, 1403,-2329,-93,145,161,-42,-57,35,-10,-1,2*0,-371,ALMGL075
     I 876,-9,-1062,90,60,82,-34,50,-13,-13,2*0,-1431,-1662,-456,231,   ALMGL076
     J -175,334,82,-144,170,-120,88,2*0,520,253,-698,-589,66,-4,235,-90,ALMGL077
     K -11,8,-18,2*0,-219,-14,188,-652,-301,-60,83,3,34,-8,17,2*0,224,  ALMGL078
     L 159,-261,50,-12,176,1,-60,-7,-39,-2,2*0,5,9,255,-119,33,84,23,-17ALMGL079
     M ,43,-36,5,2*0,-96,1,43,75,-33,49,90,-64,-15,47,17,2*0,-50,-21,3, ALMGL080
     N -79,5,10,-36,-43,-42,37,16,2*0,66,54,3,35,-3,-1,45,-5,75,-46,31, ALMGL081
     O 2*0,-61,-64,2,5,-63,-7,7,-3,-2,-45,-23,28*0/                     ALMGL082
      DATA GTT2 /1000,-62,-154,-123,1,45,-6,-14,6,-5,-3,2*0,-43,114,-18,ALMGL083
     P-27,-44,1,15,-6,8,-1,-3,2*0,54,-16,-253,28,17,75,10,-34,39,-27,20,ALMGL084
     Q 2*0,95,-7,79,-183,7,8,50,-4,-8,5,-8,2*0,4,56,-35,-47,-97,15,-11, ALMGL085
     R -6,15,-7,7,2*0,-46,7,-7,1,-24,56,26,-27,-2,-6,1,2*0,20,-11,15,   ALMGL086
     S -29,29,-10,23,-1,5,-9,1,2*0,-14,16,14,5,-8,16,11,-4,-8,6,1,2*0,  ALMGL087
     T -15,-12,5,-11,0,-3,-9,-3,-7,5,5,2*0,22,7,-2,9,6,-1,9,-4,19,-9,4, ALMGL088
     U 2*0,-12,-14,1,1,-11,-1,1,-1,1,-6,-2,28*0/                        ALMGL089
                                                                                
C     CAIN / LANGEL * 143-TERM POGO 10/68 EPOCH 1960.0                          
                                                                                
      DATA G3 / 10, -304650,-15414,13258,9591,-2343,491,759,74,110,-26, ALMGL090
     A 23,0,57910,-21633,29763,-19837,8196,3577,545,-524,60,66,-20,-18, ALMGL091
     B 0,-19772,1566,16075,13169,4864,2339,48,80,-81,18,10,-21,0,-4453, ALMGL092
     C 2334,-949,8420,-3724,-210,-2491,100,-92,-125,-55,55,0,1354,-2667,ALMGL093
     D 207,-2415,2562,-1471,17,-367,-8,158,-7,-15,0,169,1133,-1287,-1151ALMGL094
     E ,1303,-452,-37,-83,91,17,75,24,0,-96,1064,568,-272,-149,-43,-916,ALMGL095
     F 66,-114,26,78,-35,0,-579,-250,-8,63,95,-117,-376,-227,79,87,17,  ALMGL096
     G -13,0,101,-130,115,-164,55,223,-49,-262,351,51,-53,25,0,-204,144,ALMGL097
     H 6,-15,14,34,148,24,-9,-24,13,-12,0,11,9,-3,75,-23,14,-5,43,80,   ALMGL098
     I -137,-27,127,0,-8,44,-1,-39,-6,18,-32,8,-59,-17,105,50,14*0/     ALMGL099
      DATA GT3 / 100,2542,-2390,-559,-62,272,-61,-89,61,-24,-1,3,0,-466,ALMGL100
     J 988,350,-1152,-251,48,106,-21,-12,30,-9,11,0,-707,-1070,-214,-441ALMGL101
     K ,-122,317,62,-108,87,4,12,5,0,848,68,-1489,287,-296,-246,396,70, ALMGL102
     L -33,4,19,-30,0,345,-39,-87,-652,86,-89,-94,107,-14,-40,-20,1,0,5,ALMGL103
     M 300,32,311,-635,-315,149,96,-85,-28,-2,-34,0,-26,-48,258,-80,50, ALMGL104
     N 82,-167,101,99,-57,-43,48,0,-87,-46,-102,25,188,-243,232,523,81, ALMGL105
     O -132,-33,52,0,-15,-10,-122,-26,15,-37,29,91,-498,-14,103,-19,0,  ALMGL106
     P -38,16,67,-14,-83,130,-33,-38,99,50,22,-3,0,21,5,54,-26,-30,-3,  ALMGL107
     Q -39,-2,-104,79,46,-165,0,35,-26,-17,17,18,-50,23,-34,37,22,-155, ALMGL108
     R -40,14*0/                                                        ALMGL109
      DATA GTT3 /1,168*0/                                               ALMGL110
                                                                                
C     CAIN / SWEENEY * 120-TERM POGO 8/69 EPOCH 1960.0                          
                                                                                
      DATA G4 / 10,-304708,-15425,13334,9647,-2375,448,793,99,96,-17,   ALMGL111
     A 2*0,57571,-21702,29893,-19826,8108,3566,594,-516,32,93,-22,2*0,  ALMGL112
     B -19793,2661,15559,12922,5068,2498,-37,-3,-56,31,13,2*0,-4249,    ALMGL113
     C 2417,-1740,8336,-3978,-143,-2324,89,-165,-120,16,2*0,1344,-3037, ALMGL114
     D 194,-2764,2247,-1497,96,-335,-33,153,-22,2*0,51,1080,-1073,-1083,ALMGL115
     E 1171,-757,20,-33,50,7,94,2*0,-76,1181,583,-181,-270,1,-831,100,  ALMGL116
     F -120,8,87,2*0,-544,-212,-87,55,151,-236,-278,39,102,4,3,2*0,98,  ALMGL117
     G -162,99,-189,106,206,-2,-207,187,62,-24,2*0,-254,128,31,-25,-21, ALMGL118
     H 73,127,47,7,-38,-1,2*0,29,35,-7,66,-50,10,-28,21,42,-88,53,28*0/ ALMGL119
      DATA GT4 / 100,2682,-2366,-724,-157,359,12,-160,19,17,-3,2*0,225, ALMGL120
     I 1003,150,-1142,-118,58,38,-26,27,-8,-8,2*0,-684,-2832,792,84,    ALMGL121
     J -536,-27,235,72,33,-46,17,2*0,449,-96,177,327,102,-326,128,86,83,ALMGL122
     K -9,-87,2*0,369,564,-109,-205,834,-108,-277,84,42,-37,-12,2*0,234,ALMGL123
     L 401,-424,63,-503,504,8,-57,0,-3,-33,2*0,-65,-238,249,-170,234,   ALMGL124
     M -259,-130,101,49,-48,-33,2*0,-168,-114,58,123,94,40,60,-140,73,  ALMGL125
     N 54,-21,2*0,1,39,-106,-9,-49,56,-67,-8,-148,-13,27,2*0,48,42,17,  ALMGL126
     O -41,-22,21,1,-113,16,33,49,2*0,-14,-37,51,-2,4,-19,7,40,-53,31,  ALMGL127
     P -75,28*0/                                                        ALMGL128
      DATA GTT4 /1,168*0/                                               ALMGL129
                                                                                
C     IGRF 1965.0 * 80-TERM 10/68 EPOCH 1965.0                                  
                                                                                
      DATA G5 / 1, -30339,-1654,1297,958,-223,47,71,10,4*0,5758,-2123,  ALMGL130
     A 2994,-2036,805,357,60,-54,9,4*0,-2006,130,1567,1289,492,246,4,0, ALMGL131
     B -3,4*0,-403,242,-176,843,-392,-26,-229,12,-12,4*0,149,-280,8,-265ALMGL132
     C ,256,-161,3,-25,-4,4*0,16,125,-123,-107,77,-51,-4,-9,7,4*0,-14,  ALMGL133
     D 106,68,-32,-10,-13,-112,13,-5,4*0,-57,-27,-8,9,23,-19,-17,-2,12, ALMGL134
     E 4*0,3,-13,5,-17,4,22,-3,-16,6,56*0/                              ALMGL135
      DATA GT5 / 10, 153,-244,2,-7,19,-1,-5,1,4*0,-23,87,3,-108,2,11,-3,ALMGL136
     F -3,4,4*0,-118,-167,-16,7,-30,29,11,-7,6,4*0,42,7,-77,-38,-1,6,19,ALMGL137
     G -5,5*0,-1,16,29,-42,-21,0,-4,3,5*0,23,17,-24,8,-3,13,-4,0,-1,4*0,ALMGL138
     H-9,-4,20,-11,1,9,-2,-2,3,4*0,-11,3,4,2,4,2,3,-6,-3,4*0,1,-2,-3,-2,ALMGL139
     I -3,-4,-3,-3,-5,56*0/                                             ALMGL140
      DATA GTT5 /1,168*0/                                               ALMGL141
                                                                                
C     LEATON, MALIN + EVANS 1965 * 80-TERM EPOCH 1965.0                         
                                                                                
      DATA G6 / 1, -30375,-1648,1164,930,-179,42,77,11,4*0,5769,-2087,  ALMGL142
     A 2954,-2033,811,357,55,-56,23,4*0,-1995,116,1579,1299,490,248,12, ALMGL143
     B 8,-6,4*0,-389,230,-141,880,-402,-20,-239,5,-17,4*0,142,-276,5,   ALMGL144
     C -264,262,-171,16,-35,5,4*0,30,135,-123,-100,84,-64,8,-16,20,4*0, ALMGL145
     D -18,101,60,-32,-27,-12,-110,9,-1,4*0,-47,-35,-9,2,27,-17,-24,2,  ALMGL146
     E 12,4*0,5,-7,3,-20,8,26,10,-12,7,56*0/                            ALMGL147
      DATA GT6 / 10, 155,-266,0,6,8,7*0,6,83,-13,-95,10,4,-5,6*0,-114,  ALMGL148
     F -182,13,-19,-22,16,18,6*0,32,16,-85,-6,2,-3,14,6*0,30,-7,27,-27, ALMGL149
     G -30,-11,6,6*0,19,23,-18,14,5,17,2,6*0,-22,2,9,-21,-1,-2,-22,84*0/ALMGL150
      DATA GTT6 /1,168*0/                                               ALMGL151

C     HURWITZ (U S COAST / GEODETIC SURVEY) * 168-TERM EPOCH 1970.0             
                                                                                
      DATA G7/10,-302059,-17917,12899,9475,-2145,460,734,121,107,-39,16,ALMGL152
     A -4,57446,-20664,29971,-20708,8009,3595,651,-546,77,57,-26,-31,30,ALMGL153
     B -20582,430,16086,12760,4579,2490,95,46,-32,23,7,-36,5,-3699,2456,ALMGL154
     C -1880,8334,-3960,-290,-2188,175,-124,-110,-19,37,-3,1617,-2758,  ALMGL155
     D 185,-2788,2436,-1669,20,-210,-44,131,-15,-3,-13,157,1420,-1310,  ALMGL156
     E -911,808,-582,-22,-32,45,33,74,-6,4,-171,1146,625,-323,-78,38,   ALMGL157
     F -1125,143,34,2,46,-8,-14,-666,-265,-34,81,209,-240,-186,41,125,  ALMGL158
     G 15,6,1,-12,121,-160,22,-176,46,189,-46,-187,94,9,-8,2,-12,-174,  ALMGL159
     H 163,14,-27,-32,80,137,-4,-14,-4,22,-24,-1,27,19,0,35,-45,22,-31, ALMGL160
     I 56,-1,-63,14,4,10,-2,26,-26,-9,21,-1,18,-14,-28,-17,-14,6,-4,-3, ALMGL161
     J 4,9,-1,-10,26,-32,13,-6,-19,7,19,12/                             ALMGL162
      DATA GT7/10,231,-244,-19,-7,12,-7,0,3,4*0,-46,112,-1,-90,-6,7,6,  ALMGL163
     K -3,3,4*0,-104,-166,40,-20,-36,12,14,3,4,4*0,72,21,-52,-54,-11,0, ALMGL164
     L 17,6,1,4*0,22,-5,14,-24,-23,-15,6,3,-1,4*0,1,25,-14,9,1,11,-3,2, ALMGL165
     M -3,4*0,-5,11,2,-3,7,22,-5,1,9,4*0,-17,-3,7,1,-2,-3,-2,-1,-2,4*0, ALMGL166
     N 2,-6,-3,-4,1,-2,-2,-1,6,56*0/                                    ALMGL167
      DATA GTT7 /1,168*0/                                               ALMGL168
                                                                        
C     LANGEL FIELD COEFFICIENTS - 120 TERM POGO 8/71 EPOCH 1960              

C     THIS WILL BE MODEL 13 IF NEEDED
                                                                        
C       DATA G13/10,-304609,-15437,13085,9598,-2233,468,725,77,122,-25,   
C      A 2*0,58089,-21750,29974,-19844,8125,3588,583,-511,45,64,-22,2*0,  
C      B -19882,2124,15676,13110,5060,2408,-15,46,-59,30,33,2*0,-4408,    
C      C 2776,-1449,8684,-3826,-236,-2420,127,-119,-121,-26,2*0,1308,     
C      D -2806,32,-2678,2740,-1513,-18,-352,-18,151,-20,2*0,124,1156,     
C      E -1128,-1205,933,-488,-32,-74,70,15,90,2*0,-67,1085,681,-210,-250,
C      F -225,-719,147,-185,3,100,2*0,-547,-259,-74,133,212,-188,-320,17, 
C      G 182,10,-15,2*0,107,-135,60,-182,124,234,2,-231,209,68,-19,2*0,   
C      H -200,155,13,-52,5,94,152,-4,-78,28,40,2*0,30,19,11,61,-56,6,-50, 
C      I 50,6,-35,-11,28*0/
C                                               
C       DATA GT13 / 10,245,-234,-32,-9,12,0,-4,4,-1,0,2*0,-63,104,3,-111,  
C      J -12,1,4,-3,1,2,0,2*0,-50,-203,39,-16,-44,10,18,0,4,-4,-1,2*0,71, 
C      K 15,-17,-33,-12,-14,27,2,1,1,-2,2*0,38,21,9,-29,13,-9,-11,11,1,-4,
C      L 0,2*0,12,25,-36,19,-15,16,7,0,-2,-1,-2,2*0,-6,-9,12,-13,20,8,-27,
C      M 2,14,-3,-5,2*0,-16,-2,4,3,1,-1,13,-9,-2,3,0,2*0,-2,0,-5,-1,-7,0, 
C      N -6,3,-15,-1,1,2*0,-1,1,3,-1,-5,-1,-4,-3,12,-5,-1,2*0,-1,0,3,0,   
C      O 1,-1,3,0,0,-3,1,28*0/
C
C     IGRF 1980 FIELD COEFFICIENTS (MODEL = 8)

      DATA G8/10,-299880,-19970,12790,9380,-2190,490,700,200,60,-30,0,0,
     A 56060,-19570,30280,-21810,7830,3570,650,-590,70,110,-40,0,0,
     B -21290,-1990,16620,12510,3980,2610,420,20,10,20,20,0,0,
     C -3350,2710,-2520,8330,-4190,-740,-1920,200,-110,-120,-50,0,0,
     D 2120,-2570,530,-2980,1990,-1620,40,-130,-70,90,-20,0,0,
     E 460,1490,-1500,-780,920,-480,140,10,40,-30,50,0,0,
     F -150,930,710,-430,-20,170,-1080,110,30,-10,30,0,0,
     G -830,-280,-50,160,180,-230,-100,-20,70,70,10,0,0,
     H 70,-180,40,-220,90,160,-130,-150,-10,10,20,0,0,
     I -210,160,90,-50,-70,90,100,-60,20,-50,30,0,0,
     J 10,10,20,50,-40,-10,-20,40,-10,-60,29*0/

      DATA GT8 /10,224,-183,0,-14,15,4,-10,8,0,0,2*0,
     A -159,113,32,-65,-14,4,0,-8,-2,0,0,2*0,
     B -127,-252,70,-7,-82,-8,34,4,-3,0,0,2*0,
     C 2,27,-79,10,-18,-33,8,5,3,0,0,2*0,
     D 46,16,29,4,-50,2,8,16,-8,0,0,2*0,
     E 18,-4,0,13,21,14,3,1,-2,0,0,2*0,
     F -5,-14,0,-16,5,0,-1,1,7,0,0,2*0,
     G -4,4,2,14,-5,-1,11,0,-3,0,0,2*0,
     H -1,-7,0,-8,2,2,-11,8,12,0,0,2*0,
     I 52*0/
                                            
      DATA GTT8 / 1,168*0/                                              
                                                                        
C     IGRF 1975 FIELD COEFFICIENTS (MODEL = 9)                               
                                                                        
      DATA G9 /1,-30186,-1898,1299,951,-204,46,66,11,4*0,5735,-2036,    
     A2997,-2144,807,368,57,-57,13,4*0,-2124,-37,1551,1296,462,275,15,  
     B-7,3,4*0,-361,249,-253,805,-393,-20,-210,7,-12,4*0,148,-264,37,   
     C-307,235,-161,-1,-22,-4,4*0,39,142,-147,-99,74,-38,-8,-9,6,4*0,   
     D-23,102,88,-43,-9,-4,-114,11,-2,4*0,-68,-24,-4,11,27,-17,-14,-8,  
     E9,4*0,4,-15,2,-19,1,18,-6,-19,1,56*0/
                             
      DATA GT9 /10,256,-249,-38,-2,3,2,0,2,4*0,-102,100,7,-104,-20,-7,  
     F5,0,3,4*0,-30,-189,43,-41,-30,11,20,6*0,69,25,-50,-42,-21,-16,28, 
     G6,2,4*0,50,8,17,-10,-31,-5,0,9,-4,4*0,12,23,-20,13,11,10,9,3,-3,  
     H4*0,-5,-1,-2,-13,7,17,-1,3,6,4*0,-14,-1,3,3,-7,1,8,-5,-3,4*0,-2,  
     I-4,-2,-3,4,-3,-6,3,-1,56*0/
                                       
      DATA GTT9 /1,168*0/                                               
                                                                        
C     BARRACLOUGH FIELD COEFFICIENTS (MODEL = 10)                        
                                                                        
      DATA G10/10,-301036,-19067,12782,9469,-2206,441,715,110,93,-50,28,
     A-5,56826,-20165,30099,-21420,7925,3514,699,-533,51,100,-33,-19,8, 
     B-20647,-581,16330,12547,4438,2623,277,23,-26,16,24,-45,-13,-3298, 
     C2659,-2270,8310,-4039,-638,-1943,134,-126,-114,-60,29,1,1934,-2658
     D,530,-2852,2125,-1575,-9,-64,-138,106,-14,-12,-6,245,1484,-1613,  
     E-834,923,-402,38,32,-1,6,66,13,0,-112,1004,776,-403,-79,156,-1087,
     F170,-24,-2,46,-8,-18,-766,-247,-45,70,245,-218,-129,-59,123,6,12, 
     G19,-16,49,-139,50,-180,57,145,-111,-167,49,5,-18,34,-7,-196,157,49
     H,-31,-42,97,122,-2,3,5,34,-16,-9,13,20,26,28,-36,-2,3,32,30,-34,  
     I-10,17,4,5,7,-9,-15,3,6,-21,9,-25,-7,3,25,0,1,-5,4,0,8,-1,-2,-5,3,
     J-20,14,-2,11/
                                                     
      DATA GT10/10,268,-250,-38,-9,2,6,-4,4,4*0,-101,100,3,-105,-22,-10,
     K9,-2,3,4*0,-28,-189,55,-47,-40,13,23,-5,0,4*0,72,28,-64,-47,-21,  
     L-21,35,3,4,4*0,54,7,26,-7,-46,-6,0,8,-2,4*0,9,26,-27,13,11,13,8,6,
     M-4,4*0,-3,-2,2,-16,4,20,-4,5,6,4*0,-12,-2,0,3,-6,0,12,-8,-3,4*0,  
     N-2,-3,-3,-3,5,-5,-6,5,0,56*0/
                                     
      DATA GTT10/100,70,-20,-28,0,-13,7,6*0,-49,0,0,0,-17,-14,4,6*0,68, 
     O0,16,-32,0,-14,12,6*0,13,0,0,0,-14,-10,11,6*0,30,0,0,16,-17,0,0,  
     P6*0,-14,10,0,0,10,0,9,6*0,0,0,-10,88*0/                           
                                                                        
C     AWC FIELD COEFFICIENTS (MODEL = 11)                                
                                                                        
      DATA G11/10,-300557,-19320,12671,9538,-2142,417,743,124,125,-55,  
     A44,-7,56705,-20170,30013,-21272,7861,357,642,-497,84,52,-20,-31,  
     B0,-20444,-692,16197,12594,4378,2561,183,54,-37,15,20,-27,16,      
     C-3435,2632,-2085,8180,-4128,-428,-1994,248,-117,-74,-27,50,4,1967,
     D-2570,201,-2875,2323,-1667,32,-118,-85,119,-40,-3,-15,314,1507,   
     E-1374,-819,862,-589,106,-37,59,23,91,-13,-12,-199,1056,608,-392,  
     F18,146,-1085,150,29,-4,42,-10,-4,-730,-279,-52,82,139,-210,-93,   
     G3,71,29,-22,10,-36,69,-154,48,-170,110,155,-96,-191,31,16,-8,9,6, 
     H-172,172,42,-1,-49,91,122,-29,8,-8,50,-21,8,2,6,-5,43,-40,27,-20, 
     I43,1,-45,21,8,-3,10,24,-21,24,32,17,-7,-30,-21,2,-1,-5,13,13,-18, 
     J9,1,-18,14,-11,1,19,-33,6,15,7/
                                   
      DATA GT11 /10,244,-249,-37,5,3,-3,4,5*0,-103,99,12,-104,-18,-3,1, 
     K2,4,4*0,-31,-190,31,-34,-37,10,17,6,1,4*0,67,21,-35,-37,-21,-12,  
     L21,9,-1,4*0,47,10,9,-14,-16,-4,-1,9,-5,4*0,15,20,-13,13,10,6,9,0, 
     M-3,4*0,-7,0,-6,-10,11,15,2,1,5,4*0,-15,0,5,2,-8,2,4,-2,-4,4*0,-2, 
     N-4,0,-3,4,-2,-7,0,-2,56*0/
                                        
      DATA GTT11 /1,168*0/

C     IGRF 1985 FIELD COEFFICIENTS (MODEL = 12)

      DATA G12 / 10,-298770,-20730,13000,9370,-2150,520,750,210,50,-40,
     $2*0,54970,-19030,30450,-22080,7800,3560,650,-610,60,100,-40,2*0,
     $-21910,-3090,16910,12440,3630,2530,500,2,0,10,20,2*0,
     $-3120,2840,-2960,8350,-4260,-940,-1860,240,-110,-120,-50,2*0,
     $2330,-2500,680,-2980,1690,-1610,40,-60,-90,90,-20,2*0,
     $470,1480,-1550,-750,950,-480,170,40,20,-30,50,2*0,
     $-160,900,690,-500,-40,200,-1020,90,40,-10,30,2*0,
     $-820,-260,-10,230,170,-210,-60,0,40,70,10,2*0,
     $70,-210,50,-250,110,120,-160,-100,-60,20,20,2*0,
     $-210,160,90,-50,-60,90,100,-50,20,-50,30,2*0,
     $10,0,30,60,-40,0,-10,40,0,-60,0,28*0/

      DATA GT12 / 10,232,-137,51,1,13,14,2,7,4*0,
     $-245,100,34,-46,-6,1,-3,-6,0,4*0,
     $-115,-202,70,-6,-78,-15,17,-5,3,4*0,
     $53,23,-108,1,-14,-32,6,8,4,4*0,
     $38,22,25,9,-68,1,0,10,-3,4*0,
     $1,-2,-1,6,0,-1,9,4,-3,4*0,
     $-4,-11,-8,-23,-5,-1,12,-5,1,4*0,
     $2,10,11,19,3,2,9,-1,-5,4*0,
     $1,-10,1,-8,2,-8,-1,13,-8,4*0,52*0/

      DATA GTT12 / 1,168*0/
                                               
C     IGRF 1990 COEFFICIENTS (MODEL = 13)

      DATA G13/10,-297754,-21358,13146,9389,-2110,607,766,224,44,-36,
     *0,0,
     A54109,-18510,30582,-22402,7823,3525,639,-642,51,99,-39,0,0,
     B-22777,-3800,16932,12456,3239,2438,604,37,-9,8,24,0,0,
     C-2865,2933,-3485,8065,-4227,-1108,-1775,275,-108,-120,-53,0,0,
     D2481,-2395,870,-2994,1417,-1656,20,9,-124,93,-24,0,0,
     E472,1535,-1544,-692,977,-370,167,57,38,-39,44,0,0,
     F-158,827,683,-525,18,269,-963,98,38,-14,30,0,0,
     G-811,-273,6,204,164,-226,-50,-5,26,73,12,0,0,
     H97,-199,71,-221,119,110,-160,-107,-60,15,22,0,0,
     I-208,154,95,-57,-64,86,91,-66,19,-55,29,0,0,
     J13,4,31,56,-42,-5,-15,38,-5,-62,29*0/

      DATA GT13/10,180,-129,33,5,6,13,6,2,4*0,
     A -161,106,24,-67,6,-1,-2,-5,-7,4*0,
     B -158,-138,0,0,-70,-16,18,-3,-2,4*0,
     C 44,16,-106,-59,5,-31,13,6,1,4*0,
     D 26,18,31,-14,-55,0,-2,16,-11,4*0,
     E -1,5,4,17,4,23,1,2,0,4*0,
     F 2,-13,0,-9,5,12,12,2,0,4*0,
     G 6,2,8,-5,-2,0,0,3,-5,4*0,
     H 5,-2,3,3,4,-5,-3,6,-6,4*0,52*0/

      DATA GTT13 /1,168*0/

C     IGRF 1995 COEFFICIENTS (MODEL = 14)

      DATA G14/1,-29682,-2197,1329,941,-210,66,78,24,4,-3,0,0,
     A5318,-1789,3074,-2268,782,352,64,-67,4,9,-4,0,0,
     B-2356,-425,1685,1249,291,237,65,1,-1,1,2,0,0,
     C-263,302,-406,769,-421,-122,-172,29,-9,-12,-5,0,0,
     D262,-232,98,-301,116,-167,2,4,-14,9,-2,0,0,
     E44,157,-152,-64,99,-26,17,8,4,-4,4,0,0,
     F-16,77,67,-57,4,28,-94,10,5,-2,3,0,0,
     G-77,-25,3,22,16,-23,-3,-2,0,7,1,0,0,
     H12,-20,7,-21,12,10,-17,-10,-7,0,3,0,0,
     I-19,15,11,-7,-7,9,7,-8,1,-6,3,0,0,
     J2,1,3,6,-4,0,-2,3,-1,-6,29*0/

      DATA GT14/10,176,-132,15,8,8,5,-2,3,4*0,
     A -183,130,37,-64,9,1,-4,-8,-2,4*0,
     B -150,-88,-8,-2,-69,-15,6,-6,1,4*0,
     C 41,22,-121,-81,5,-20,19,6,4,4*0,
     D 18,12,27,-10,-46,-1,-2,12,-11,4*0,
     E 2,12,3,18,9,23,-2,1,3,4*0,
     F 3,-16,-2,-9,10,22,0,2,2,4*0,
     G 8,2,6,-4,0,-3,0,-6,-9,4*0,
     H 4,-2,2,7,0,-12,-7,-6,-3,4*0,52*0/

      DATA GTT14 /1,168*0/
C                                               
      DATA SHMIT(1,1)/0.0/,TMOLD/0.0/,MODOLD /0/,RAD/57.29578/          ALMGL169
C  *****  NON-SUBSCRIPTED, FIXED-INDEX VERSION BEGINS HERE (NO DO-LOOPS)ALMGL170
C  *****  BEGIN PROGRAM                                                 ALMGL171
      IF(SHMIT(1,1).EQ.-1.)   GO TO 8                                   ALMGL172
C  *****  INITIALIZE * ONCE ONLY, FIRST TIME SUBROUTINE IS CALLED       ALMGL173
      SHMIT(1,1)=-1.                                                    ALMGL174
      DO 2 N=2,13                                                       ALMGL175
      SHMIT(N,1) = (2*N-3) * SHMIT(N-1,1) / (N-1)                       ALMGL176
      JJ=2                                                              ALMGL177
      DO 2 M=2,N                                                        ALMGL178
      SHMIT(N,M) = SHMIT(N,M-1) * DSQRT(1.0D0*FLOAT((N-M+1)*JJ)/(N+M-2))ALMGL179
      SHMIT(M-1,N)=SHMIT(N,M)                                           ALMGL180
    2 JJ = 1                                                            ALMGL181
      DO 7 K=1,14                                                       ALMGL182
      F1=LG(1,1,K)                                                      ALMGL183
      F2=LGT(1,1,K)                                                     ALMGL184
      F3=LGTT(1,1,K)                                                    ALMGL185
      NMAX=NMX(K)                                                       ALMGL186
      L = 0                                                             ALMGL187
      DO 3 I=1,3                                                        ALMGL188
    3 ISUM(K,I) = 0                                                     ALMGL189
      DO 4 N=1,NMAX                                                     ALMGL190
      DO 4 M=1,NMAX                                                     ALMGL191
      L = L+1                                                           ALMGL192
      ISUM(K,1)=ISUM(K,1)+L*LG(N,M,K)                                   ALMGL193
      ISUM(K,2)=ISUM(K,2)+L*LGT(N,M,K)                                  ALMGL194
    4 ISUM(K,3)=ISUM(K,3)+L*LGTT(N,M,K)                                 ALMGL195
      DO 6 I=1,3                                                        ALMGL196
      IF(ISUM(K,I).EQ.LSUM(K,I))  GO TO 6                               ALMGL197
C  *****  ERROR IN DATA CARDS - NOTE WRITE AND STOP STATEMENTS          ALMGL198
      PRINT 5,   K,I,LSUM(K,I),ISUM(K,I)                                ALMGL199
    5 FORMAT(///29H DATA WRONG IN ALLMAG--MODEL ,I2,3X,2HI=,I1,3X,      ALMGL200
     A17HPRECALCULATED SUM,I10,3X,17HTHIS MACHINE GETS,I10)             ALMGL201
      STOP                                                              ALMGL202
    6 CONTINUE                                                          ALMGL203
      DO 7 N=1,NMAX                                                     ALMGL204
      DO 7 M=1,NMAX                                                     ALMGL205
      GG(N,M,K)=LG(N,M,K)*SHMIT(N,M)/F1                                 ALMGL206
      GGT(N,M,K)=LGT(N,M,K)*SHMIT(N,M)/F2                               ALMGL207
    7 GGTT(N,M,K)=LGTT(N,M,K)*SHMIT(N,M)/F3                             ALMGL208
    8 IF((MODEL.EQ.MODOLD).AND.(TM.EQ.TMOLD))  GO TO 11                 ALMGL209
C  *****  NOTE WRITE STATEMENT - NEW MODEL OR NEW TIME                  ALMGL210
C      PRINT 9,   MODEL,TM                                               ALMGL211
C    9 FORMAT('0 MODEL USED IS NUMBER ',I2,2X,'  FOR TM =',F9.3/)        ALMGL212
      IF(MODEL.LT.1.OR.MODEL.GT.14) STOP                                ALMGL213
      MODOLD=MODEL                                                      ALMGL214
      TMOLD=TM                                                          ALMGL215
      NMAX=NMX(MODEL)                                                   ALMGL216
      T=TM-T0(MODEL)                                                    ALMGL217
      DO 10 N=1,NMAX                                                    ALMGL218
      DO 10 M=1,NMAX                                                    ALMGL219
   10 G(N,M)=GG(N,M,MODEL)+T*(GGT(N,M,MODEL)+GGTT(N,M,MODEL)*T)         ALMGL220
C  *****  CALCULATION USUALLY BEGINS
      WLONG = -RAD * DATAN (G(1,2) / G(2,2))
      COLAT =  RAD * DATAN (SQRT(G(1,2)**2 + G(2,2)**2) / G(2,1))
      EM = DSQRT (G(1,2)**2 + G(2,2)**2 + G(2,1)**2)
      CONSTEM = EM / 100000.0
C      PRINT 19, WLONG, COLAT, EM
C   19 FORMAT(5X,'GEOGRAPHIC COORDINATES OF BOREAL MAGNETIC DIPOLE POLE'/
C     $10X,'WEST LONGITUDE =',F9.3/10X,'GEOC. COLATITUDE =',F9.3/
C     $10X,'EARTH''S MAGNETIC MOMENT =',F8.0,' GAMMA'/)
   11 P21=CT                                                            ALMGL222
      P22=ST                                                            ALMGL223
      AR=6371.2/RKM                                                     ALMGL224
      SP2=SPH                                                           ALMGL225
      CP2=CPH                                                           ALMGL226
      DP21=-P22                                                         ALMGL227
      DP22=P21                                                          ALMGL228
      AOR=AR*AR*AR                                                      ALMGL229
      C2=G(2,2)*CP2+G(1,2)*SP2                                          ALMGL230
      BR=-(AOR+AOR)*(G(2,1)*P21+C2*P22)                                 ALMGL231
      BT=AOR*(G(2,1)*DP21+C2*DP22)                                      ALMGL232
      BP=AOR*(G(1,2)*CP2-G(2,2)*SP2)*P22                                ALMGL233
      IF (NMAX.LE. 2) GO TO 1                                           ALMGL234
C                                                           N= 3        ALMGL235
      SP3=(SP2+SP2)*CP2                                                 ALMGL236
      CP3=(CP2+SP2)*(CP2-SP2)                                           ALMGL237
      P31=P21*P21-0.333333333                                           ALMGL238
      P32=P21*P22                                                       ALMGL239
      P33=P22*P22                                                       ALMGL240
      DP31=-P32-P32                                                     ALMGL241
      DP32=P21*P21-P33                                                  ALMGL242
      DP33=-DP31                                                        ALMGL243
      AOR=AOR*AR                                                        ALMGL244
      C2=G(3,2)*CP2+G(1,3)*SP2                                          ALMGL245
      C3=G(3,3)*CP3+G(2,3)*SP3                                          ALMGL246
      BR=BR-3.0*AOR*(G(3,1)*P31+C2*P32+C3*P33)                          ALMGL247
      BT=BT+AOR*(G(3,1)*DP31+C2*DP32+C3*DP33)                           ALMGL248
      BP=BP-AOR*((G(3,2)*SP2-G(1,3)*CP2)*P32+2.0*(G(3,3)*SP3-G(2,3)*CP3)ALMGL249
     1*P33)                                                             ALMGL250
      IF (NMAX.LE. 3) GO TO 1                                           ALMGL251
C                                                           N= 4        ALMGL252
      SP4=SP2*CP3+CP2*SP3                                               ALMGL253
      CP4=CP2*CP3-SP2*SP3                                               ALMGL254
      P41=P21*P31-0.26666666*P21                                        ALMGL255
      DP41=P21*DP31+DP21*P31-0.26666666*DP21                            ALMGL256
      P42=P21*P32-0.20000000*P22                                        ALMGL257
      DP42=P21*DP32+DP21*P32-0.20000000*DP22                            ALMGL258
      P43=P21*P33                                                       ALMGL259
      DP43=P21*DP33+DP21*P33                                            ALMGL260
      P44=P22*P33                                                       ALMGL261
      DP44=3.0*P43                                                      ALMGL262
      AOR=AOR*AR                                                        ALMGL263
      C2=G(4,2)*CP2+G(1,4)*SP2                                          ALMGL264
      C3=G(4,3)*CP3+G(2,4)*SP3                                          ALMGL265
      C4=G(4,4)*CP4+G(3,4)*SP4                                          ALMGL266
      BR=BR-4.0*AOR*(G(4,1)*P41+C2*P42+C3*P43+C4*P44)                   ALMGL267
      BT=BT+AOR*(G(4,1)*DP41+C2*DP42+C3*DP43+C4*DP44)                   ALMGL268
      BP=BP-AOR*((G(4,2)*SP2-G(1,4)*CP2)*P42+2.0*(G(4,3)*SP3-G(2,4)*CP3)ALMGL269
     1*P43+3.0*(G(4,4)*SP4-G(3,4)*CP4)*P44)                             ALMGL270
      IF (NMAX.LE. 4) GO TO 1                                           ALMGL271
C                                                           N= 5        ALMGL272
      SP5=(SP3+SP3)*CP3                                                 ALMGL273
      CP5=(CP3+SP3)*(CP3-SP3)                                           ALMGL274
      P51=P21*P41-0.25714285*P31                                        ALMGL275
      DP51=P21*DP41+DP21*P41-0.25714285*DP31                            ALMGL276
      P52=P21*P42-0.22857142*P32                                        ALMGL277
      DP52=P21*DP42+DP21*P42-0.22857142*DP32                            ALMGL278
      P53=P21*P43-0.14285714*P33                                        ALMGL279
      DP53=P21*DP43+DP21*P43-0.14285714*DP33                            ALMGL280
      P54=P21*P44                                                       ALMGL281
      DP54=P21*DP44+DP21*P44                                            ALMGL282
      P55=P22*P44                                                       ALMGL283
      DP55=4.0*P54                                                      ALMGL284
      AOR=AOR*AR                                                        ALMGL285
      C2=G(5,2)*CP2+G(1,5)*SP2                                          ALMGL286
      C3=G(5,3)*CP3+G(2,5)*SP3                                          ALMGL287
      C4=G(5,4)*CP4+G(3,5)*SP4                                          ALMGL288
      C5=G(5,5)*CP5+G(4,5)*SP5                                          ALMGL289
      BR=BR-5.0*AOR*(G(5,1)*P51+C2*P52+C3*P53+C4*P54+C5*P55)            ALMGL290
      BT=BT+AOR*(G(5,1)*DP51+C2*DP52+C3*DP53+C4*DP54+C5*DP55)           ALMGL291
      BP=BP-AOR*((G(5,2)*SP2-G(1,5)*CP2)*P52+2.0*(G(5,3)*SP3-G(2,5)*CP3)ALMGL292
     1*P53+3.0*(G(5,4)*SP4-G(3,5)*CP4)*P54+4.0*(G(5,5)*SP5-G(4,5)*CP5)*PALMGL293
     255)                                                               ALMGL294
      IF (NMAX.LE. 5) GO TO 1                                           ALMGL295
C                                                           N= 6        ALMGL296
      SP6=SP2*CP5+CP2*SP5                                               ALMGL297
      CP6=CP2*CP5-SP2*SP5                                               ALMGL298
      P61=P21*P51-0.25396825*P41                                        ALMGL299
      DP61=P21*DP51+DP21*P51-0.25396825*DP41                            ALMGL300
      P62=P21*P52-0.23809523*P42                                        ALMGL301
      DP62=P21*DP52+DP21*P52-0.23809523*DP42                            ALMGL302
      P63=P21*P53-0.19047619*P43                                        ALMGL303
      DP63=P21*DP53+DP21*P53-0.19047619*DP43                            ALMGL304
      P64=P21*P54-0.11111111*P44                                        ALMGL305
      DP64=P21*DP54+DP21*P54-0.11111111*DP44                            ALMGL306
      P65=P21*P55                                                       ALMGL307
      DP65=P21*DP55+DP21*P55                                            ALMGL308
      P66=P22*P55                                                       ALMGL309
      DP66=5.0*P65                                                      ALMGL310
      AOR=AOR*AR                                                        ALMGL311
      C2=G(6,2)*CP2+G(1,6)*SP2                                          ALMGL312
      C3=G(6,3)*CP3+G(2,6)*SP3                                          ALMGL313
      C4=G(6,4)*CP4+G(3,6)*SP4                                          ALMGL314
      C5=G(6,5)*CP5+G(4,6)*SP5                                          ALMGL315
      C6=G(6,6)*CP6+G(5,6)*SP6                                          ALMGL316
      BR=BR-6.0*AOR*(G(6,1)*P61+C2*P62+C3*P63+C4*P64+C5*P65+C6*P66)     ALMGL317
      BT=BT+AOR*(G(6,1)*DP61+C2*DP62+C3*DP63+C4*DP64+C5*DP65+C6*DP66)   ALMGL318
      BP=BP-AOR*((G(6,2)*SP2-G(1,6)*CP2)*P62+2.0*(G(6,3)*SP3-G(2,6)*CP3)ALMGL319
     1*P63+3.0*(G(6,4)*SP4-G(3,6)*CP4)*P64+4.0*(G(6,5)*SP5-G(4,6)*CP5)*PALMGL320
     265+5.0*(G(6,6)*SP6-G(5,6)*CP6)*P66)                               ALMGL321
      IF (NMAX.LE. 6) GO TO 1                                           ALMGL322
C                                                           N= 7        ALMGL323
      SP7=(SP4+SP4)*CP4                                                 ALMGL324
      CP7=(CP4+SP4)*(CP4-SP4)                                           ALMGL325
      P71=P21*P61-0.25252525*P51                                        ALMGL326
      DP71=P21*DP61+DP21*P61-0.25252525*DP51                            ALMGL327
      P72=P21*P62-0.24242424*P52                                        ALMGL328
      DP72=P21*DP62+DP21*P62-0.24242424*DP52                            ALMGL329
      P73=P21*P63-0.21212121*P53                                        ALMGL330
      DP73=P21*DP63+DP21*P63-0.21212121*DP53                            ALMGL331
      P74=P21*P64-0.16161616*P54                                        ALMGL332
      DP74=P21*DP64+DP21*P64-0.16161616*DP54                            ALMGL333
      P75=P21*P65-0.09090909*P55                                        ALMGL334
      DP75=P21*DP65+DP21*P65-0.09090909*DP55                            ALMGL335
      P76=P21*P66                                                       ALMGL336
      DP76=P21*DP66+DP21*P66                                            ALMGL337
      P77=P22*P66                                                       ALMGL338
      DP77=6.0*P76                                                      ALMGL339
      AOR=AOR*AR                                                        ALMGL340
      C2=G(7,2)*CP2+G(1,7)*SP2                                          ALMGL341
      C3=G(7,3)*CP3+G(2,7)*SP3                                          ALMGL342
      C4=G(7,4)*CP4+G(3,7)*SP4                                          ALMGL343
      C5=G(7,5)*CP5+G(4,7)*SP5                                          ALMGL344
      C6=G(7,6)*CP6+G(5,7)*SP6                                          ALMGL345
      C7=G(7,7)*CP7+G(6,7)*SP7                                          ALMGL346
      BR=BR-7.0*AOR*(G(7,1)*P71+C2*P72+C3*P73+C4*P74+C5*P75+C6*P76+C7*P7ALMGL347
     17)                                                                ALMGL348
      BT=BT+AOR*(G(7,1)*DP71+C2*DP72+C3*DP73+C4*DP74+C5*DP75+C6*DP76+C7*ALMGL349
     1DP77)                                                             ALMGL350
      BP=BP-AOR*((G(7,2)*SP2-G(1,7)*CP2)*P72+2.0*(G(7,3)*SP3-G(2,7)*CP3)ALMGL351

     1*P73+3.0*(G(7,4)*SP4-G(3,7)*CP4)*P74+4.0*(G(7,5)*SP5-G(4,7)*CP5)*PALMGL352
     275+5.0*(G(7,6)*SP6-G(5,7)*CP6)*P76+6.0*(G(7,7)*SP7-G(6,7)*CP7)*P77ALMGL353
     3)                                                                 ALMGL354
      IF (NMAX.LE. 7) GO TO 1                                           ALMGL355
C                                                           N= 8        ALMGL356
      SP8=SP2*CP7+CP2*SP7                                               ALMGL357
      CP8=CP2*CP7-SP2*SP7                                               ALMGL358
      P81=P21*P71-0.25174825*P61                                        ALMGL359
      DP81=P21*DP71+DP21*P71-0.25174825*DP61                            ALMGL360
      P82=P21*P72-0.24475524*P62                                        ALMGL361
      DP82=P21*DP72+DP21*P72-0.24475524*DP62                            ALMGL362
      P83=P21*P73-0.22377622*P63                                        ALMGL363
      DP83=P21*DP73+DP21*P73-0.22377622*DP63                            ALMGL364
      P84=P21*P74-0.18881118*P64                                        ALMGL365
      DP84=P21*DP74+DP21*P74-0.18881118*DP64                            ALMGL366
      P85=P21*P75-0.13986013*P65                                        ALMGL367
      DP85=P21*DP75+DP21*P75-0.13986013*DP65                            ALMGL368
      P86=P21*P76-0.07692307*P66                                        ALMGL369
      DP86=P21*DP76+DP21*P76-0.07692307*DP66                            ALMGL370
      P87=P21*P77                                                       ALMGL371
      DP87=P21*DP77+DP21*P77                                            ALMGL372
      P88=P22*P77                                                       ALMGL373
      DP88=7.0*P87                                                      ALMGL374
      AOR=AOR*AR                                                        ALMGL375
      C2=G(8,2)*CP2+G(1,8)*SP2                                          ALMGL376
      C3=G(8,3)*CP3+G(2,8)*SP3                                          ALMGL377
      C4=G(8,4)*CP4+G(3,8)*SP4                                          ALMGL378
      C5=G(8,5)*CP5+G(4,8)*SP5                                          ALMGL379
      C6=G(8,6)*CP6+G(5,8)*SP6                                          ALMGL380
      C7=G(8,7)*CP7+G(6,8)*SP7                                          ALMGL381
      C8=G(8,8)*CP8+G(7,8)*SP8                                          ALMGL382
      BR=BR-8.0*AOR*(G(8,1)*P81+C2*P82+C3*P83+C4*P84+C5*P85+C6*P86+C7*P8ALMGL383
     17+C8*P88)                                                         ALMGL384
      BT=BT+AOR*(G(8,1)*DP81+C2*DP82+C3*DP83+C4*DP84+C5*DP85+C6*DP86+C7*ALMGL385
     1DP87+C8*DP88)                                                     ALMGL386
      BP=BP-AOR*((G(8,2)*SP2-G(1,8)*CP2)*P82+2.0*(G(8,3)*SP3-G(2,8)*CP3)ALMGL387
     1*P83+3.0*(G(8,4)*SP4-G(3,8)*CP4)*P84+4.0*(G(8,5)*SP5-G(4,8)*CP5)*PALMGL388
     285+5.0*(G(8,6)*SP6-G(5,8)*CP6)*P86+6.0*(G(8,7)*SP7-G(6,8)*CP7)*P87ALMGL389
     3+7.0*(G(8,8)*SP8-G(7,8)*CP8)*P88)                                 ALMGL390
      IF (NMAX.LE. 8) GO TO 1                                           ALMGL391
C                                                           N= 9        ALMGL392
      SP9=(SP5+SP5)*CP5                                                 ALMGL393
      CP9=(CP5+SP5)*(CP5-SP5)                                           ALMGL394
      P91=P21*P81-0.25128205*P71                                        ALMGL395
      DP91=P21*DP81+DP21*P81-0.25128205*DP71                            ALMGL396
      P92=P21*P82-0.24615384*P72                                        ALMGL397
      DP92=P21*DP82+DP21*P82-0.24615384*DP72                            ALMGL398
      P93=P21*P83-0.23076923*P73                                        ALMGL399
      DP93=P21*DP83+DP21*P83-0.23076923*DP73                            ALMGL400
      P94=P21*P84-0.20512820*P74                                        ALMGL401
      DP94=P21*DP84+DP21*P84-0.20512820*DP74                            ALMGL402
      P95=P21*P85-0.16923076*P75                                        ALMGL403
      DP95=P21*DP85+DP21*P85-0.16923076*DP75                            ALMGL404
      P96=P21*P86-0.12307692*P76                                        ALMGL405
      DP96=P21*DP86+DP21*P86-0.12307692*DP76                            ALMGL406
      P97=P21*P87-0.06666666*P77                                        ALMGL407
      DP97=P21*DP87+DP21*P87-0.06666666*DP77                            ALMGL408
      P98=P21*P88                                                       ALMGL409
      DP98=P21*DP88+DP21*P88                                            ALMGL410
      P99=P22*P88                                                       ALMGL411
      DP99=8.0*P98                                                      ALMGL412
      AOR=AOR*AR                                                        ALMGL413
      C2=G(9,2)*CP2+G(1,9)*SP2                                          ALMGL414
      C3=G(9,3)*CP3+G(2,9)*SP3                                          ALMGL415
      C4=G(9,4)*CP4+G(3,9)*SP4                                          ALMGL416
      C5=G(9,5)*CP5+G(4,9)*SP5                                          ALMGL417
      C6=G(9,6)*CP6+G(5,9)*SP6                                          ALMGL418
      C7=G(9,7)*CP7+G(6,9)*SP7                                          ALMGL419
      C8=G(9,8)*CP8+G(7,9)*SP8                                          ALMGL420
      C9=G(9,9)*CP9+G(8,9)*SP9                                          ALMGL421
      BR=BR-9.0*AOR*(G(9,1)*P91+C2*P92+C3*P93+C4*P94+C5*P95+C6*P96+C7*P9ALMGL422
     17+C8*P98+C9*P99)                                                  ALMGL423
      BT=BT+AOR*(G(9,1)*DP91+C2*DP92+C3*DP93+C4*DP94+C5*DP95+C6*DP96+C7*ALMGL424
     1DP97+C8*DP98+C9*DP99)                                             ALMGL425
      BP=BP-AOR*((G(9,2)*SP2-G(1,9)*CP2)*P92+2.0*(G(9,3)*SP3-G(2,9)*CP3)ALMGL426
     1*P93+3.0*(G(9,4)*SP4-G(3,9)*CP4)*P94+4.0*(G(9,5)*SP5-G(4,9)*CP5)*PALMGL427
     295+5.0*(G(9,6)*SP6-G(5,9)*CP6)*P96+6.0*(G(9,7)*SP7-G(6,9)*CP7)*P97ALMGL428
     3+7.0*(G(9,8)*SP8-G(7,9)*CP8)*P98+8.0*(G(9,9)*SP9-G(8,9)*CP9)*P99) ALMGL429
      IF (NMAX.LE. 9) GO TO 1                                           ALMGL430
C                                                           N=10        ALMGL431
      SP10=SP2*CP9+CP2*SP9                                              ALMGL432
      CP10=CP2*CP9-SP2*SP9                                              ALMGL433
      P101=P21*P91-0.25098039*P81                                       ALMGL434
      DP101=P21*DP91+DP21*P91-0.25098039*DP81                           ALMGL435
      P102=P21*P92-0.24705882*P82                                       ALMGL436
      DP102=P21*DP92+DP21*P92-0.24705882*DP82                           ALMGL437
      P103=P21*P93-0.23529411*P83                                       ALMGL438
      DP103=P21*DP93+DP21*P93-0.23529411*DP83                           ALMGL439
      P104=P21*P94-0.21568627*P84                                       ALMGL440
      DP104=P21*DP94+DP21*P94-0.21568627*DP84                           ALMGL441
      P105=P21*P95-0.18823529*P85                                       ALMGL442
      DP105=P21*DP95+DP21*P95-0.18823529*DP85                           ALMGL443
      P106=P21*P96-0.15294117*P86                                       ALMGL444
      DP106=P21*DP96+DP21*P96-0.15294117*DP86                           ALMGL445
      P107=P21*P97-0.10980392*P87                                       ALMGL446
      DP107=P21*DP97+DP21*P97-0.10980392*DP87                           ALMGL447
      P108=P21*P98-0.05882352*P88                                       ALMGL448
      DP108=P21*DP98+DP21*P98-0.05882352*DP88                           ALMGL449
      P109=P21*P99                                                      ALMGL450
      DP109=P21*DP99+DP21*P99                                           ALMGL451
      P1010=P22*P99                                                     ALMGL452
      DP1010=9.0*P109                                                   ALMGL453
      AOR=AOR*AR                                                        ALMGL454
      C2=G(10,2)*CP2+G(1,10)*SP2                                        ALMGL455
      C3=G(10,3)*CP3+G(2,10)*SP3                                        ALMGL456
      C4=G(10,4)*CP4+G(3,10)*SP4                                        ALMGL457
      C5=G(10,5)*CP5+G(4,10)*SP5                                        ALMGL458
      C6=G(10,6)*CP6+G(5,10)*SP6                                        ALMGL459
      C7=G(10,7)*CP7+G(6,10)*SP7                                        ALMGL460
      C8=G(10,8)*CP8+G(7,10)*SP8                                        ALMGL461
      C9=G(10,9)*CP9+G(8,10)*SP9                                        ALMGL462
      C10=G(10,10)*CP10+G(9,10)*SP10                                    ALMGL463
      BR=BR-10.0*AOR*(G(10,1)*P101+C2*P102+C3*P103+C4*P104+C5*P105+C6*P1ALMGL464
     106+C7*P107+C8*P108+C9*P109+C10*P1010)                             ALMGL465
      BT=BT+AOR*(G(10,1)*DP101+C2*DP102+C3*DP103+C4*DP104+C5*DP105+C6*DPALMGL466
     1106+C7*DP107+C8*DP108+C9*DP109+C10*DP1010)                        ALMGL467
      BP=BP-AOR*((G(10,2)*SP2-G(1,10)*CP2)*P102+2.0*(G(10,3)*SP3-G(2,10)ALMGL468
     1*CP3)*P103+3.0*(G(10,4)*SP4-G(3,10)*CP4)*P104+4.0*(G(10,5)*SP5-G(4ALMGL469
     2,10)*CP5)*P105+5.0*(G(10,6)*SP6-G(5,10)*CP6)*P106+6.0*(G(10,7)*SP7ALMGL470
     3-G(6,10)*CP7)*P107+7.0*(G(10,8)*SP8-G(7,10)*CP8)*P108+8.0*(G(10,9)ALMGL471
     4*SP9-G(8,10)*CP9)*P109+9.0*(G(10,10)*SP10-G(9,10)*CP10)*P1010)    ALMGL472
      IF (NMAX.LE.10) GO TO 1                                           ALMGL473
C                                                           N=11        ALMGL474
      SP11=(SP6+SP6)*CP6                                                ALMGL475
      CP11=(CP6+SP6)*(CP6-SP6)                                          ALMGL476
      P111=P21*P101-0.25077399*P91                                      ALMGL477
      DP111=P21*DP101+DP21*P101-0.25077399*DP91                         ALMGL478
      P112=P21*P102-0.24767801*P92                                      ALMGL479
      DP112=P21*DP102+DP21*P102-0.24767801*DP92                         ALMGL480
      P113=P21*P103-0.23839009*P93                                      ALMGL481
      DP113=P21*DP103+DP21*P103-0.23839009*DP93                         ALMGL482
      P114=P21*P104-0.22291021*P94                                      ALMGL483
      DP114=P21*DP104+DP21*P104-0.22291021*DP94                         ALMGL484
      P115=P21*P105-0.20123839*P95                                      ALMGL485
      DP115=P21*DP105+DP21*P105-0.20123839*DP95                         ALMGL486
      P116=P21*P106-0.17337461*P96                                      ALMGL487
      DP116=P21*DP106+DP21*P106-0.17337461*DP96                         ALMGL488
      P117=P21*P107-0.13931888*P97                                      ALMGL489
      DP117=P21*DP107+DP21*P107-0.13931888*DP97                         ALMGL490
      P118=P21*P108-0.09907120*P98                                      ALMGL491
      DP118=P21*DP108+DP21*P108-0.09907120*DP98                         ALMGL492
      P119=P21*P109-0.05263157*P99                                      ALMGL493
      DP119=P21*DP109+DP21*P109-0.05263157*DP99                         ALMGL494
      P1110=P21*P1010                                                   ALMGL495
      DP1110=P21*DP1010+DP21*P1010                                      ALMGL496
      P1111=P22*P1010                                                   ALMGL497
      DP1111=10.0*P1110                                                 ALMGL498
      AOR=AOR*AR                                                        ALMGL499
      C2=G(11,2)*CP2+G(1,11)*SP2                                        ALMGL500
      C3=G(11,3)*CP3+G(2,11)*SP3                                        ALMGL501
      C4=G(11,4)*CP4+G(3,11)*SP4                                        ALMGL502
      C5=G(11,5)*CP5+G(4,11)*SP5                                        ALMGL503
      C6=G(11,6)*CP6+G(5,11)*SP6                                        ALMGL504
      C7=G(11,7)*CP7+G(6,11)*SP7                                        ALMGL505
      C8=G(11,8)*CP8+G(7,11)*SP8                                        ALMGL506
      C9=G(11,9)*CP9+G(8,11)*SP9                                        ALMGL507
      C10=G(11,10)*CP10+G(9,11)*SP10                                    ALMGL508
      C11=G(11,11)*CP11+G(10,11)*SP11                                   ALMGL509
      BR=BR-11.0*AOR*(G(11,1)*P111+C2*P112+C3*P113+C4*P114+C5*P115+C6*P1ALMGL510
     116+C7*P117+C8*P118+C9*P119+C10*P1110+C11*P1111)                   ALMGL511
      BT=BT+AOR*(G(11,1)*DP111+C2*DP112+C3*DP113+C4*DP114+C5*DP115+C6*DPALMGL512
     1116+C7*DP117+C8*DP118+C9*DP119+C10*DP1110+C11*DP1111)             ALMGL513
      BP=BP-AOR*((G(11,2)*SP2-G(1,11)*CP2)*P112+2.0*(G(11,3)*SP3-G(2,11)ALMGL514
     1*CP3)*P113+3.0*(G(11,4)*SP4-G(3,11)*CP4)*P114+4.0*(G(11,5)*SP5-G(4ALMGL515
     2,11)*CP5)*P115+5.0*(G(11,6)*SP6-G(5,11)*CP6)*P116+6.0*(G(11,7)*SP7ALMGL516
     3-G(6,11)*CP7)*P117+7.0*(G(11,8)*SP8-G(7,11)*CP8)*P118+8.0*(G(11,9)ALMGL517
     4*SP9-G(8,11)*CP9)*P119+9.0*(G(11,10)*SP10-G(9,11)*CP10)*P1110+10.0ALMGL518
     5*(G(11,11)*SP11-G(10,11)*CP11)*P1111)                             ALMGL519
      IF (NMAX.LE.11) GO TO 1                                           ALMGL520
C                                                           N=12        ALMGL521
      SP12=SP2*CP11+CP2*SP11                                            ALMGL522
      CP12=CP2*CP11-SP2*SP11                                            ALMGL523
      P121=P21*P111-0.25062656*P101                                     ALMGL524
      DP121=P21*DP111+DP21*P111-0.25062656*DP101                        ALMGL525
      P122=P21*P112-0.24812030*P102                                     ALMGL526
      DP122=P21*DP112+DP21*P112-0.24812030*DP102                        ALMGL527
      P123=P21*P113-0.24060150*P103                                     ALMGL528
      DP123=P21*DP113+DP21*P113-0.24060150*DP103                        ALMGL529
      P124=P21*P114-0.22807017*P104                                     ALMGL530
      DP124=P21*DP114+DP21*P114-0.22807017*DP104                        ALMGL531
      P125=P21*P115-0.21052631*P105                                     ALMGL532
      DP125=P21*DP115+DP21*P115-0.21052631*DP105                        ALMGL533
      P126=P21*P116-0.18796992*P106                                     ALMGL534
      DP126=P21*DP116+DP21*P116-0.18796992*DP106                        ALMGL535
      P127=P21*P117-0.16040100*P107                                     ALMGL536
      DP127=P21*DP117+DP21*P117-0.16040100*DP107                        ALMGL537
      P128=P21*P118-0.12781954*P108                                     ALMGL538
      DP128=P21*DP118+DP21*P118-0.12781954*DP108                        ALMGL539
      P129=P21*P119-0.09022556*P109                                     ALMGL540
      DP129=P21*DP119+DP21*P119-0.09022556*DP109                        ALMGL541
      P1210=P21*P1110-0.04761904*P1010                                  ALMGL542
      DP1210=P21*DP1110+DP21*P1110-0.04761904*DP1010                    ALMGL543
      P1211=P21*P1111                                                   ALMGL544
      DP1211=P21*DP1111+DP21*P1111                                      ALMGL545
      P1212=P22*P1111                                                   ALMGL546
      DP1212=11.0*P1211                                                 ALMGL547
      AOR=AOR*AR                                                        ALMGL548
      C2=G(12,2)*CP2+G(1,12)*SP2                                        ALMGL549
      C3=G(12,3)*CP3+G(2,12)*SP3                                        ALMGL550
      C4=G(12,4)*CP4+G(3,12)*SP4                                        ALMGL551
      C5=G(12,5)*CP5+G(4,12)*SP5                                        ALMGL552
      C6=G(12,6)*CP6+G(5,12)*SP6                                        ALMGL553
      C7=G(12,7)*CP7+G(6,12)*SP7                                        ALMGL554
      C8=G(12,8)*CP8+G(7,12)*SP8                                        ALMGL555
      C9=G(12,9)*CP9+G(8,12)*SP9                                        ALMGL556
      C10=G(12,10)*CP10+G(9,12)*SP10                                    ALMGL557
      C11=G(12,11)*CP11+G(10,12)*SP11                                   ALMGL558
      C12=G(12,12)*CP12+G(11,12)*SP12                                   ALMGL559
      BR=BR-12.0*AOR*(G(12,1)*P121+C2*P122+C3*P123+C4*P124+C5*P125+C6*P1ALMGL560
     126+C7*P127+C8*P128+C9*P129+C10*P1210+C11*P1211+C12*P1212)         ALMGL561
      BT=BT+AOR*(G(12,1)*DP121+C2*DP122+C3*DP123+C4*DP124+C5*DP125+C6*DPALMGL562
     1126+C7*DP127+C8*DP128+C9*DP129+C10*DP1210+C11*DP1211+C12*DP1212)  ALMGL563
      BP=BP-AOR*((G(12,2)*SP2-G(1,12)*CP2)*P122+2.0*(G(12,3)*SP3-G(2,12)ALMGL564
     1*CP3)*P123+3.0*(G(12,4)*SP4-G(3,12)*CP4)*P124+4.0*(G(12,5)*SP5-G(4ALMGL565
     2,12)*CP5)*P125+5.0*(G(12,6)*SP6-G(5,12)*CP6)*P126+6.0*(G(12,7)*SP7ALMGL566
     3-G(6,12)*CP7)*P127+7.0*(G(12,8)*SP8-G(7,12)*CP8)*P128+8.0*(G(12,9)ALMGL567
     4*SP9-G(8,12)*CP9)*P129+9.0*(G(12,10)*SP10-G(9,12)*CP10)*P1210+10.0ALMGL568
     5*(G(12,11)*SP11-G(10,12)*CP11)*P1211+11.0*(G(12,12)*SP12-G(11,12)*ALMGL569
     6CP12)*P1212)                                                      ALMGL570
      IF (NMAX.LE.12) GO TO 1                                           ALMGL571
C                                                           N=13        ALMGL572
      SP13=(SP7+SP7)*CP7                                                ALMGL573
      CP13=(CP7+SP7)*(CP7-SP7)                                          ALMGL574
      P131=P21*P121-0.25051759*P111                                     ALMGL575
      DP131=P21*DP121+DP21*P121-0.25051759*DP111                        ALMGL576
      P132=P21*P122-0.24844720*P112                                     ALMGL577
      DP132=P21*DP122+DP21*P122-0.24844720*DP112                        ALMGL578
      P133=P21*P123-0.24223602*P113                                     ALMGL579
      DP133=P21*DP123+DP21*P123-0.24223602*DP113                        ALMGL580
      P134=P21*P124-0.23188405*P114                                     ALMGL581
      DP134=P21*DP124+DP21*P124-0.23188405*DP114                        ALMGL582
      P135=P21*P125-0.21739130*P115                                     ALMGL583
      DP135=P21*DP125+DP21*P125-0.21739130*DP115                        ALMGL584
      P136=P21*P126-0.19875776*P116                                     ALMGL585
      DP136=P21*DP126+DP21*P126-0.19875776*DP116                        ALMGL586
      P137=P21*P127-0.17598343*P117                                     ALMGL587
      DP137=P21*DP127+DP21*P127-0.17598343*DP117                        ALMGL588
      P138=P21*P128-0.14906832*P118                                     ALMGL589
      DP138=P21*DP128+DP21*P128-0.14906832*DP118                        ALMGL590
      P139=P21*P129-0.11801242*P119                                     ALMGL591
      DP139=P21*DP129+DP21*P129-0.11801242*DP119                        ALMGL592
      P1310=P21*P1210-0.08281573*P1110                                  ALMGL593
      DP1310=P21*DP1210+DP21*P1210-0.08281573*DP1110                    ALMGL594
      P1311=P21*P1211-0.04347826*P1111                                  ALMGL595
      DP1311=P21*DP1211+DP21*P1211-0.04347826*DP1111                    ALMGL596
      P1312=P21*P1212                                                   ALMGL597
      DP1312=P21*DP1212+DP21*P1212                                      ALMGL598
      P1313=P22*P1212                                                   ALMGL599
      DP1313=12.0*P1312                                                 ALMGL600
      AOR=AOR*AR                                                        ALMGL601
      C2=G(13,2)*CP2+G(1,13)*SP2                                        ALMGL602
      C3=G(13,3)*CP3+G(2,13)*SP3                                        ALMGL603
      C4=G(13,4)*CP4+G(3,13)*SP4                                        ALMGL604
      C5=G(13,5)*CP5+G(4,13)*SP5                                        ALMGL605
      C6=G(13,6)*CP6+G(5,13)*SP6                                        ALMGL606
      C7=G(13,7)*CP7+G(6,13)*SP7                                        ALMGL607
      C8=G(13,8)*CP8+G(7,13)*SP8                                        ALMGL608
      C9=G(13,9)*CP9+G(8,13)*SP9                                        ALMGL609
      C10=G(13,10)*CP10+G(9,13)*SP10                                    ALMGL610
      C11=G(13,11)*CP11+G(10,13)*SP11                                   ALMGL611
      C12=G(13,12)*CP12+G(11,13)*SP12                                   ALMGL612
      C13=G(13,13)*CP13+G(12,13)*SP13                                   ALMGL613
      BR=BR-13.0*AOR*(G(13,1)*P131+C2*P132+C3*P133+C4*P134+C5*P135+C6*P1ALMGL614
     136+C7*P137+C8*P138+C9*P139+C10*P1310+C11*P1311+C12*P1312+C13*P1313ALMGL615
     2)                                                                 ALMGL616
      BT=BT+AOR*(G(13,1)*DP131+C2*DP132+C3*DP133+C4*DP134+C5*DP135+C6*DPALMGL617
     1136+C7*DP137+C8*DP138+C9*DP139+C10*DP1310+C11*DP1311+C12*DP1312+C1ALMGL618
     23*DP1313)                                                         ALMGL619
      BP=BP-AOR*((G(13,2)*SP2-G(1,13)*CP2)*P132+2.0*(G(13,3)*SP3-G(2,13)ALMGL620
     1*CP3)*P133+3.0*(G(13,4)*SP4-G(3,13)*CP4)*P134+4.0*(G(13,5)*SP5-G(4ALMGL621
     2,13)*CP5)*P135+5.0*(G(13,6)*SP6-G(5,13)*CP6)*P136+6.0*(G(13,7)*SP7ALMGL622
     3-G(6,13)*CP7)*P137+7.0*(G(13,8)*SP8-G(7,13)*CP8)*P138+8.0*(G(13,9)ALMGL623
     4*SP9-G(8,13)*CP9)*P139+9.0*(G(13,10)*SP10-G(9,13)*CP10)*P1310+10.0ALMGL624
     5*(G(13,11)*SP11-G(10,13)*CP11)*P1311+11.0*(G(13,12)*SP12-G(11,13)*ALMGL625
     6CP12)*P1312+12.0*(G(13,13)*SP13-G(12,13)*CP13)*P1313)             ALMGL626
    1 BR = BR / 100000.                                                 ALMGL627
      BT = BT / 100000.                                                 ALMGL628
      BP = BP / ST / 100000.                                            ALMGL629
      B =  DSQRT(BR*BR+BT*BT+BP*BP)                                     ALMGL632
      RETURN                                                            ALMGL634
      END                                                               ALMGL635
