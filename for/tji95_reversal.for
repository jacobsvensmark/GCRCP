      PROGRAM tji95_reversal
C
C........+.........+.........+.........+.........+.........+.........+..
C     Multi-platform COSMIC-RAY TRAJECTORY PROGRAM
C     FORTRAN 77     transportable version
C     Read in control card; LAT, LON, RIG, ZENITH, AZIMUTH, DELPC, INDO
C          Then calculate        INDO    trajectories
C               Starting at      PC
C               Incrementing at  DELPC   intervals
C     Includes conversion from Geodetic to Geocentric coordinates
C     Includes re-entrant albedo calculations
C     Uses subroutine SINGLTJE to do trajectory calculations
C     Magnetic field - IGRF 1995 (order 10)                          ###
C........+.........+.........+.........+.........+.........+.........+..
C     Restrictions: Cannot run over N or S pole; will get BETA blowup
C........+.........+.........+.........+.........+.........+.........+..
C     Mod History
CLast Mod 21 Dec 00  Make all intrinsic function double precision for PC
C     Mod 20 Dec 00  Insert 8 character format 1000 with AZ & ZE
C     Mod 17 Feb 99  set limit to 600000
C     Mod 17 Feb 99  if (ymax.lt.6.6) IFATE = 3
C     Mod    Aug 97  Adjust step size to minimize beta problems
C     Mod    Jan 97  High latitude step size adjust, introduce AHLT
C     Mod    Jun 96  EDIF limit set to 1.0e-5
C     Mod    Jun 96  IERRPT formats, Boundary and look ahead 
C     Mod    Feb 96  Standard reference TJ1V line check
C     Mod    Dec 94  Print out start and end times of PC run
C     **************************************************************
C          Timing estimates base on COMPAQ Digital FORTRAN
C     Will run on PIII PC at 850 MHZ        55000 steps/sec (Real*8)
C     Will run on PIII PC at 700 MHZ        39000 steps/sec (Real*8)
C     Will run on PIII PC at 550 MHZ        32000 steps/sec (Real*8)
C     Will run on PII  PC at 400 MHZ        23000 steps/sec (Real*8)
C     **************************************************************
C     *  TAPE*       Monitor program operation
C     *  TAPE1       Trajectory control cards
C     *  TAPE7       80  character line (card image)  output
C     *  TAPE8       132 character line printer output
C     *  TAPE16      Diagnostic output for trouble shooting
C     *              Normally turned off (open statement commented out)
C     **************************************************************
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
      REAL TIME,TIME_TARGET
      INTEGER JDATA,MGNMAX
      REAL GMSUM
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ The following used for timing PC runs
C        Can use on PC (IBM and COMPACQ FORTRAN) and on SUN's
C        Cannot use on IBM SP2 or DEC VAX
C     ...+.........+.........+.........+.........+.........+.........+..
C
      INTEGER*2 ISYEAR,ISMONTH,ISDAY, ISHOUR,ISMIN,ISSEC,ISHSEC
      INTEGER*2 IEYEAR,IEMONTH,IEDAY, IEHOUR,IEMIN,IESEC,IEHSEC
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
      COMMON /GHSV/ G(4,4),GMSUM,JDATA,MGNMAX
C........+.........+.........+.........+.........+.........+.........+..
C
      CHARACTER*1 CFF
      CHARACTER outfile*120 ! Added by JSV, May 2021
      CHARACTER infile*120 ! Added by JSV, May 2021
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Use CFF / '1'/  on 'Main Fames"; 
C        use CFF /Z'0C'/ for files to be printed on word processors
C........+.........+.........+.........+.........+.........+.........+..
C
      DATA CFF / '1'/
C     DATA CFF /Z'0C'/
C
C........+.........+.........+.........+.........+.........+.........+..
      CALL getarg(1,outfile) ! Added by JSV, May 2021
      CALL getarg(2,infile) ! Added by JSV, May 2021

      !OPEN (1, FILE='TAPE1', STATUS='OLD') ! Commented by JSV, May 2021
      OPEN (1, FILE=infile, STATUS='OLD') ! Added by JSV, May 2021
      OPEN (7, FILE=outfile, STATUS='UNKNOWN') ! Added by JSV, May 2021
      OPEN (8, FILE='TAPE8', STATUS='UNKNOWN')
      OPEN (16,FILE='TAPE16',STATUS='UNKNOWN')

      OPEN (UNIT=191,FILE='time_target.txt')
          READ(UNIT=191,FMT='(F6.2)') TIME_TARGET
      CLOSE(191) 
      TIME = 0.0
      OPEN (UNIT=192,FILE='reversal_gh3order.txt',STATUS='UNKNOWN')
          READ (192,*)
          DO WHILE (ABS(TIME-TIME_TARGET) .GT. 0.01)
              READ (192,FMT='(F12.6,F13.2,15F12.2)') TIME,
     c            G(1,1),G(2,1),G(3,1),G(4,1),
     c            G(1,2),G(2,2),G(3,2),G(4,2),
     c            G(1,3),G(2,3),G(3,3),G(4,3),
     c            G(1,4),G(2,4),G(3,4),G(4,4)
         END DO
      CLOSE(192)
      MGNMAX = 4 ! JSV: Set to order of G
      JDATA = 0 ! Set to 0
      GMSUM = 1.0e-3 ! s

C 1000 FORMAT (A1,' RUN START DATE ',  I4, '/',I2, '/',I2,'@',I2,':',I2,
C     *           ':',I2)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ User defined program control 
C........+.........+.........+.........+.........+.........+.........+..

      FSTEP = 4.0E08
      LIMIT = 600000
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
      PI   = ACOS(-1.0)
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
C                Geocentric coordinates used for output
C        All calculation are done in Geocentric coordinates!
C     \/ Conversion from Geodetic to Geocentric coordinates
C.......+.........+.........+.........+.........+.........+.........+..
C
      CALL GDGC (TCD, TSD)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Remember positron of initial point on trajectory
C                 in Geocentric coordinates
C        Y(1) is distance in earth radii from geocenter
C             Start with height above geoid and convert to earth radii
C                   The initial values of Y(1), Y(2), and Y(3) are
C                   calculated in subroutine GDGC
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
      TSGDZE = SIN(GDZE)
      TCGDZE = COS(GDZE)
      TSGDAZ = SIN(GDAZ)
      TCGDAZ = COS(GDAZ)
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
C     Trajectories are calculated in subroutine SINGLTJ
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
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Get date and time of run end  (PC routine)
C........+.........+.........+.........+.........+.........+.........+..
CUDKOMENTERET HSV
C      CALL GETDAT (IEYEAR,IEMONTH,IEDAY)
C      CALL GETTIM (IEHOUR,IEMIN,IESEC,IEHSEC)
C
      WRITE (8, 1100)  IEYEAR,IEMONTH,  IEDAY,IEHOUR,IEMIN,IESEC
      WRITE (16,1100)  IEYEAR,IEMONTH,  IEDAY,IEHOUR,IEMIN,IESEC
 1100 FORMAT (//'  RUN END   DATE ',  I4, '/',I2, '/',I2,'@',I2,':',I2,
     *          ':',I2)
      WRITE (8, 1110)  ISYEAR,ISMONTH,  ISDAY,ISHOUR,ISMIN,ISSEC
 1110 FORMAT ('  RUN START DATE ',  I4, '/',I2, '/',I2,'@',I2,':',I2,
     *        ':',I2)
C
      WRITE (*, 1120) TSTEP,NTRAJC
      WRITE (8, 1120) TSTEP,NTRAJC
      WRITE (16,1120) TSTEP,NTRAJC
 1120 FORMAT (//' TOTAL NUMBER OF STEPS        ',F15.0///
     *          ' TOTAL NUMBER OF TRAJECTORIES',I15///)
      Write (*,1130)
 1130 format (' End program TJI95T')
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
      ERADPL = SQRT(ERPLSQ)
      ERECSQ = EREQSQ/ERPLSQ - 1.0
C
      GDCLT = PIO2-GDLATD/RAD
      TSGDCLT = SIN(GDCLT)
      TCGDCLT = COS(GDCLT)
      ONE = EREQSQ*TSGDCLT*TSGDCLT
      TWO = ERPLSQ*TCGDCLT*TCGDCLT
      THREE = ONE+TWO
      RHO = SQRT(THREE)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Get geocentric distance from geocenter in kilometers
C........+.........+.........+.........+.........+.........+.........+..
C
      DISTKM = SQRT(SALT*(SALT+2.0*RHO)+(EREQSQ*ONE+ERPLSQ*TWO)/THREE)
C
C........+.........+.........+.........+.........+.........+.........+..
C     TCD and TSD are sine and cosine of the angle the Geodetic vertical
C         must be rotated to form the Geocentric vertical
C........+.........+.........+.........+.........+.........+.........+..
C
      TCD = (SALT+RHO)/DISTKM
      TSD = (EREQSQ-ERPLSQ)/RHO*TCGDCLT*TSGDCLT/DISTKM
      TCY2 = TCGDCLT*TCD-TSGDCLT*TSD
      TSY2 = TSGDCLT*TCD+TCGDCLT*TSD
C
      Y(2) = ACOS(TCY2)
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
      DATA CNAME / '  I95 '/                                            
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
      RC1O6 = 1.0/6.0
      SR2 = SQRT(2.0)
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
C........+.........+.........+.........+.........+.........+.........+..
C
      Y(1) = RY1
      Y(2) = RY2
      Y(3) = RY3
      GRNDKM = (ERADPL/SQRT(1.0-ERECSQ*TSY2SQ))
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
      TENG = SQRT((PC*ZCHARGE)**2+(ANUC*EMCSQ)**2)
      EOMC = -8987.566297*ZCHARGE/TENG
      GMA = SQRT(((PC*ZCHARGE)/(EMCSQ*ANUC))**2+1.0)
      BETA = SQRT(1.0-1.0/(GMA*GMA))
      PVEL = VEL*BETA
      HMAX = 1.0/PVEL
      disck = disout - 1.1*hmax*pvel
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Set max step length ("HMAX") to 1 earth radii
C        PVEL  is particle velocity in earth radii per second
C        DISCK is check for approaching termination boundary
C                 (within 1.1 steps)
C........+.........+.........+.........+.........+.........+.........+..
C
      EDIF = BETA*1.0E-4
      if (edif.lt.1.0-5)  edif = 1.0e-5
      if (beta.lt.0.1)    edif = 1.0e-4
C
      Y(4) = BETA*Y1GC
      Y(5) = BETA*Y2GC
      Y(6) = BETA*Y3GC
C
       azd = gdazd
       zed = gdzed
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
      CALL FGRAD
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
      CALL FGRAD
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
      CALL FGRAD
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
      CALL FGRAD
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
         TSY2SQ = SIN(Y(2))**2
         GRNDKM = (ERADPL/SQRT(1.0-ERECSQ*TSY2SQ))
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
      RCKBETA = SQRT(Y(4)*Y(4)+Y(5)*Y(5)+Y(6)*Y(6))
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
      ACCER = SQRT(F(4)*F(4)+F(5)*F(5)+F(6)*F(6))
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
            DO 220 I = 1, 6
               Y(I) = YOLD(I)
               F(I) = FOLD(I)
  220       CONTINUE
         endif
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
         IF (IERRPT.GT.3)  WRITE (16, 2140) y(1),disck,pvel,H,nstep
C
         if (h.lt.1.0e-5 .or. hck.lt.1.0e-5 .or. hcng.lt.1.0e-5)  then
            IRT = 1
            go to 260
         else
            hck  = hck/2.0
            hcng = hcng/2.0
            TAU = TAU - H
            DO 240 I = 1, 6
               Y(I) = YOLD(I)
               F(I) = FOLD(I)
  240       CONTINUE
            GO TO 130
         ENDIF
      ENDIF
C
 2140 FORMAT (' 2140 ',2x,'y(1),disck,pvel,H', 
     *        4x,1pe12.6,4x,e12.6,4x,e12.6,4x,e9.2,27x,i6)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Store values of  Y  and  F  as  FOLD & YOLD
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
C     \/ Write out results
C        IRT     +1     ALLOWED          (FATE = 0)
C        IRT      0     FAILED           (FATE = 2)
C        IRT     -1     RE-ENTRANT       (FATE = 1)
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (IRT.GT.0) THEN
         TCY2 = COS(Y(2))
         TSY2 = SIN(Y(2))
         YDA5 = Y(5)*TCY2+Y(4)*TSY2
         ATRG1 = Y(4)*TCY2-Y(5)*TSY2
         ATRG2 = SQRT(Y(6)*Y(6)+YDA5*YDA5)
         FASLAT = 0.0
         IF (ATRG1.NE.0.0.AND.ATRG2.NE.0.0) FASLAT = 
     *                                      ATAN2(ATRG1,ATRG2)*RAD
         FASLON = Y(3)*RAD
         IF (Y(6).NE.0.0.AND.YDA5.NE.0.0) FASLON = (Y(3)+ATAN2(Y(6),
     *                                    YDA5))*RAD
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
 2160 FORMAT (F7.2,F8.2,F9.3,2F6.1,I7,F7.2,3X,F8.2,3X,I7,3X,I3,3X,A6) ! Added by HSV and JSV, May 2021
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
 2180 FORMAT (F7.2,F8.2,F9.3,2F6.1,I7,4X,'R',10X,'R',3X,I9,3X,I3,3X,A6) ! Added by HSV and JSV, May 2021
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
 2200 FORMAT (F7.2,F8.2,F9.3,2F6.1,I7,4X,'F',10X,'F',3X,I9,3X,I3,3X,A6) ! Added by HSV and JSV, May 2021
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
      SUBROUTINE FGRAD
C
C........+.........+.........+.........+.........+.........+.........+..
C     Mod    Feb 96  standard reference TJ1V (line check 17 Feb)
C     Mod 27 Jan 1999  Change MAGNEW to NEWMAG95                        ###
C........+.........+.........+.........+.........+.........+.........+..
C     Programmer   -  Don F. Smart; FORTRAN77
C     Note -  The programming adheres to the conventional FORTRAN
C             default standard that variables beginning with
C            'i','j','k','l','m',or 'n' are integer variables
C             Variables beginning with "c" are character variables
C             All other variables are real
C........+.........+.........+.........+.........+.........+.........+..
C        Do not mix different type variables in same common block
C           Some computers do not allow this
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
C
C........+.........+.........+.........+.........+.........+.........+..
C
      F(1) = VEL*Y(4)
      F(2) = VEL*Y(5)/Y(1)
      TSY2 = SIN(Y(2))
      TCY2 = COS(Y(2))
      F(3) = VEL*Y(6)/(Y(1)*TSY2)
      SQY6 = Y(6)*Y(6)/Y(1)
      Y5OY1 = Y(5)/Y(1)
      TAY2 = TSY2/TCY2
      CALL MAGNEW95
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
      SUBROUTINE MAGNEW95                                               
C
C........+.........+.........+.........+.........+.........+.........+..
C     Compute Magnetic field 
C     Derived from NASA (NSSDC) routine NEWMAG version of December 1965
C             modified for 10 order field
C     Coefficients for IGRF 1995 loaded into this subroutine            ###
C     Coefficients obtained from program CNGMAGN
C........+.........+.........+.........+.........+.........+.........+..
CLast Mod 27 Jan 1999  IGRF 95 coefficients                             ###
C     Mod    Feb 1996  standard reference TJ1V (line check 17 Feb)
C     Mod    Nov 1980  for arguments in labeled common
C........+.........+.........+.........+.........+.........+.........+..
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
      INTEGER MGNMAX,JDATA
      REAL GMSUM
C
C........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /WRKVLU/ F(6),Y(6),ERAD,EOMC,VEL,BR,BT,BP,B
      COMMON /WRKTSC/ TSY2,TCY2,TSY3,TCY3
C
C........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /GHSV/ G(4,4),GMSUM,JDATA,MGNMAX
C
C........+.........+.........+.........+.........+.........+.........+..
C
C     \/ Load in data constants if this is the first time called
C        otherwise, skip to evaluation of magnetic field
C          designed to drop high order terms if contribution
C                   would be less then "BERR"
C           also designed so the maximum order of expansion
C               can be specified
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (JDATA.EQ.77) GO TO 120
C
C........+.........+.........+.........+.........+.........+.........+..
C     Gauss normalized Schmidt coefficients ordered for fast computation
C
Cards for FORTRAN
C     1995.00     Coef in CNGMAG format                             IGRF95  
C     1995.00     Coef in CNGMAG format                             IGRF95  
C
C   ******************************************************************
C   *    The array G contains Gauss normalized Schmidt coefficients
c   *    the array G contains both the G and H coefficients
C   *     G(1,1) = 0.0
C   *Schmidt G(N,M) corresponds to -G(NN+1,MM+1) Gauss normalized coef
C   *Schmidt H(N,M) corresponds to -G( MI ,NN+1) Gauss normalized coef
C   *                            where MI = M
C   ******************************************************************
C
C      WRITE(*,*) GMSUM,GSUM 
!      IF (GMSUM.EQ.0) GO TO 110
!      print *,"GMSUM must have been for this to happen"
!      P22 = 0.
!      BERR = 0.0001
!      AR = 0.
!      DO 100 L = 1, MGNMAX
!         DO 100 M = 1, MGNMAX
!            AR = AR+1.
!            P22 = P22+AR*G(M,L)
!  100 CONTINUE
!C            write(*,*)GMSUM,P22
!      GMSUM = (GSUM-P22)/GSUM
!C
!C........+.........+.........+.........+.........+.........+.........+..
!C**** \/ Note following print and stop statements
!C........+.........+.........+.........+.........+.........+.........+..
!C
!      IF (ABS(GMSUM).GT.1.E2) THEN
!         WRITE (*, 2200) GMSUM
!         WRITE (7, 2200) GMSUM
!         WRITE (8, 2200) GMSUM
!         STOP
!      ENDIF
! 2200 FORMAT (' DATA WRONG IN MAGNEW',E15.6)
C
  110 CONTINUE
C
      GMSUM = 0.
      JDATA = 77
C
  120 CONTINUE
      P21 = TCY2
      P22 = TSY2
      AR = 1.0/Y(1)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ N= 2
C........+.........+.........+.........+.........+.........+.........+..
C
      DP22 = P21
      TSY3 = SIN(Y(3))
      TCY3 = COS(Y(3))
      TSP2 = TSY3
      TCP2 = TCY3
      DP21 = -P22
      AOR = AR*AR*AR
      RC2 = G(2,2)*TCP2+G(1,2)*TSP2
      BR = -(AOR+AOR)*(G(2,1)*P21+RC2*P22)
      BT = AOR*(G(2,1)*DP21+RC2*DP22)
      BP = AOR*(G(1,2)*TCP2-G(2,2)*TSP2)*P22
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ N = 3
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (MGNMAX.LT.3) GO TO 130
      AOR = AOR*AR
      ERR = BERR*SQRT((BP/P22)**2+BR**2+BT**2)
      TSP3 = (TSP2+TSP2)*TCP2
      TCP3 = (TCP2+TSP2)*(TCP2-TSP2)
      P31 = P21*P21-0.333333333
      P32 = P21*P22
      P33 = P22*P22
      DP31 = -P32-P32
      DP32 = P21*P21-P33
      DP33 = -DP31
      RC2 = G(3,2)*TCP2+G(1,3)*TSP2
      RC3 = G(3,3)*TCP3+G(2,3)*TSP3
      BR = BR-3.0*AOR*(G(3,1)*P31+RC2*P32+RC3*P33)
      BT = BT+AOR*(G(3,1)*DP31+RC2*DP32+RC3*DP33)
      BP = BP-AOR*((G(3,2)*TSP2-G(1,3)*TCP2)*P32+
     *     2.0*(G(3,3)*TSP3-G(2,3)*TCP3)*P33)
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ N = 4
C........+.........+.........+.........+.........+.........+.........+..
C
      IF (MGNMAX.LT.4) GO TO 130
      AOR = AOR*AR
      TSP4 = TSP2*TCP3+TCP2*TSP3
      TCP4 = TCP2*TCP3-TSP2*TSP3
      P41  = P21*P31-0.26666666*P21
      DP41 = P21*DP31+DP21*P31-0.26666666*DP21
      P42  = P21*P32-0.20000000*P22
      DP42 = P21*DP32+DP21*P32-0.20000000*DP22
      P43  = P21*P33
      DP43 = P21*DP33+DP21*P33
      P44  = P22*P33
      DP44 = 3.0*P43
      RC2 = G(4,2)*TCP2+G(1,4)*TSP2
      RC3 = G(4,3)*TCP3+G(2,4)*TSP3
      RC4 = G(4,4)*TCP4+G(3,4)*TSP4
      BR = BR-4.0*AOR*(G(4,1)*P41+RC2*P42+RC3*P43+RC4*P44)
      BT = BT+AOR*(G(4,1)*DP41+RC2*DP42+RC3*DP43+RC4*DP44)
      BP = BP-AOR*((G(4,2)*TSP2-G(1,4)*TCP2)*P42+
     *     2.0*(G(4,3)*TSP3-G(2,4)*TCP3)*P43+
     *     3.0*(G(4,4)*TSP4-G(3,4)*TCP4)*P44)
C
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Convert to units of Gauss
C........+.........+.........+.........+.........+.........+.........+..
C
  130 BP = BP/P22*1.E-5
      BT = BT*1.E-5
      BR = BR*1.E-5
      B = SQRT(BR*BR+BT*BT+BP*BP)
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
