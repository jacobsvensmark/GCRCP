=========================================================================
National Space Science Data Center                         March 13, 2001
=========================================================================

NAME:                   Geomagnetic Cutoff Rigidity Computer Program

TECHNICAL CONTACTS:	Don F. Smart and Margaret (Peggy) A. Shea
                        Center for Space Plasmas and Aeronomic Research
                        The University of Alabama in Huntsville
                        Huntsville, Alabama 35889
                        Phone: (256)824-6660
                        E-mail: sssrc@msn.com

NSSDC CONTACT:		Dr. John F. Cooper
                        Raytheon ITSS / SSDOO Project
                        Code 632
                        NASA Goddard Space Flight Center
                        Greenbelt, Maryland 20771
                        Phone: (301) 286-1193
                        Fax: (301) 286-1771

                        E-mail: jfcooper@pop600.gsfc.nasa.gov

FILES:                  *.doc  

                          Final report to NASA in Microsoft Word format on
documentation 
                          for the computer codes submitted for archiving at
NSSDC.

                        Tji95.for

                          Fortran-77 source code using IGRF95 geomagnetic field
model.

                        Tji95t.for

                          Same source code as Tji95.for but also includes timing
routines
                          for use with Compaq Digital Fortran on PC desktop
computers.
                          May also work with FORTRAN compilers used on SUN and
SGI systems.

                        Tjallmag.for

                          Fortran-77 source code using the NASA ALLMAG
subroutine for a
                          variety of geomagnetic field models including IGRF95.

                        Erg-rig.for

                          Sample Fortran-77 source subroutine for conversion
from kinetic
                          energy to magnetiic rigidity.

                        Rig_erg.for

                          Sample Fortran-77 source subroutine for conversion
from magnetic
                          rigidity to kinetic energy.


BRIEF DESCRIPTION:

Prediction of magnetic rigidity (momentum/charge) cutoffs or thresholds 
for penetration of energetic cosmic ray particles from interplanetary 
space through Earth's geomagnetic field to the lower atmosphere was 
pioneered by Carl Stoermer in 1930 as the true nature of these particles 
was just becoming known. In the 1950's and thereafter the development of

increasingly fast computers allowed detailed predictions of charged 
particle trajectories through complex multipolar models of the 
geomagnetic field. From the 1970's to the present the standard software 
for such computations has been that of Margaret (Peggy) A. Shea and Don 
F. Smart, both originally at the old Air Force Geophysics Laboratory in

Massachusetts  and now affiliated with the University of Alabama at 
Huntsville. Updated FORTRAN 77 versions of this software, test files for
input and output, and accompanying documentation in Microsoft Word format, 
have recently been archived at the National Space Science Data Center 
through funding support from NASA Grant NAG5-8009. This software 
determines allowed trajectories reaching specified locations near or 
above the Earth's surface for a given range of magnetic rigidities and

incidence directions. The lower limit of the rigidity range for allowed

trajectories approximately defines the local cosmic ray  cutoff. Cutoff

averaging is needed at some locations and directions for penumbral 
effects of partial shadowing by the body of the Earth. Choices of models
for the geomagnetic field include IGRF-95 and other older models in the 
NASA ALLMAG subroutines.

----------------------------- REFERENCES --------------------------------

Included in software documentation.

=========================================================================
National Space Science Data Center

==========================================================================



