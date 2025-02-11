#include <fintrf.h>

C     Gateway function for subroutine pvawsme in pcaws.f
C     for FORTRAN subroutine pvawsme in pcaws.f
C     
C     [bi, theta, smoothedmpmData] = pvawslast(y,       1
C                                              yd,      2
C                                              pos,     3
C                                              nv,      4
C                                              nvd,     5
C                                              nd,      6
C                                              n1,      7
C                                              n2,      8
C                                              n3,      9
C                                              hakt,   10 
C                                              lambda, 11
C                                              theta,  12
C                                              bi,     13
C                                              bin,    14 create space in memory for PLHS
C                                              thnew,  15 create space in memory for PLHS
C                                              ydnew,  16 create space in memory for PLHS
C                                              invcov, 17
C                                              ncores, 18
C                                              spmin,  19
C                                              lwght,  20 create space in memory
C                                              wght,   21
C                                              swjy,   22 create space in memory
C                                              swjd,   23 create space in memory
C                                              np1,    24
C                                              np2,    25
C                                              np3);   26

      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)

      IMPLICIT NONE

      MWPOINTER PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS

C     local integer variables

      MWSIZE NV_IN
      MWSIZE NVD_IN
      MWSIZE ND_IN
      MWSIZE N1_IN
      MWSIZE N2_IN
      MWSIZE N3_IN 
      MWSIZE NCORES_IN
      MWSIZE SPMIN_IN 
      MWSIZE NP1_IN
      MWSIZE NP2_IN
      MWSIZE NP3_IN

C     local real variables (non-portable syntax)

      REAL*8 HAKT_IN
      REAL*8 LAMBDA_IN

C     local input pointers

      MWPOINTER Y_IN_PR
      MWPOINTER YD_IN_PR
      MWPOINTER POS_IN_PR
      MWPOINTER THETA_IN_PR
      MWPOINTER BI_IN_PR
      MWPOINTER BIN_IN_PR
      MWPOINTER THNEW_IN_PR
      MWPOINTER YDNEW_OUT_PR
      MWPOINTER INVCOV_IN_PR
      MWPOINTER LWGHT_IN_PR
      MWPOINTER WGHT_IN_PR
      MWPOINTER SWJY_IN_PR
      MWPOINTER SWJD_IN_PR

C     Check for proper number of arguments

      IF (NRHS .NE. 26) THEN
         CALL MEXERRMSGTXT('pvawslast requires 26 input arguments')
      ELSEIF (NLHS .NE. 4) THEN
         CALL MEXERRMSGTXT('pvawslast requires 4 output argument')
      ENDIF

C     get all the PRHS

      Y_IN_PR = MXGETPR(PRHS(1))
      YD_IN_PR = MXGETPR(PRHS(2))
      POS_IN_PR = MXGETPR(PRHS(3))
      NV_IN = INT(MXGETSCALAR(PRHS(4)))
      NVD_IN = INT(MXGETSCALAR(PRHS(5)))
      ND_IN = INT(MXGETSCALAR(PRHS(6)))
      N1_IN = INT(MXGETSCALAR(PRHS(7)))
      N2_IN = INT(MXGETSCALAR(PRHS(8)))
      N3_IN = INT(MXGETSCALAR(PRHS(9)))
      HAKT_IN = MXGETSCALAR(PRHS(10))
      LAMBDA_IN = MXGETSCALAR(PRHS(11))
      THETA_IN_PR = MXGETPR(PRHS(12))
      BI_IN_PR = MXGETPR(PRHS(13))
      BIN_IN_PR = MXGETPR(PRHS(14))          DO WE NEED THIS AT ALL??
      THNEW_IN_PR = MXGETPR(PRHS(15))        DO WE NEED THIS AT ALL??
      YDNEW_IN_PR = MXGETPR(PRHS(16))        DO WE NEED THIS AT ALL??
      INVCOV_IN_PR = MXGETPR(PRHS(17))
      NCORES_IN = INT(MXGETSCALAR(PRHS(18)))
      SPMIN_IN = INT(MXGETSCALAR(PRHS(19)))
      LWGHT_IN_PR = MXGETPR(PRHS(20))        DO WE NEED THIS AT ALL??
      WGHT_IN_PR = MXGETPR(PRHS(21))
      SWJY_IN_PR = MXGETPR(PRHS(22))         DO WE NEED THIS AT ALL??
      SWJD_IN_PR = MXGETPR(PRHS(23))         DO WE NEED THIS AT ALL??
      NP1_IN = INT(MXGETSCALAR(PRHS(24)))
      NP2_IN = INT(MXGETSCALAR(PRHS(25)))
      NP3_IN = INT(MXGETSCALAR(PRHS(26)))

C     create all the PLHS
      
      PLHS(1) = MXCREATEDOUBLEMATRIX(1, NV_IN, 0)           CORRECT?
      PLHS(2) = MXCREATEDOUBLEMATRIX(NV_IN, NV_IN, 0)       CORRECT?

      BI_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))
      YDNEW_OUT_PR = MXGETPR(PLHS(3))

      LWGHT_LOCAL = MXCREATEDOUBLEMATRIX(ONE, DWL, 0)       CORRECT? DWL UNDEFINED
      SWJY_LOCAL = MXCREATEDOUBLEMATRIX(NV, NCORES_IN, 0)   CORRECT?
      SWJD_LOCAL = MXCREATEDOUBLEMATRIX(NV, NCORES_IN, 0)   CORRECT?

      LWGHT_LOCAL_PR = MXGETPR(LWGHT_LOCAL)
      SWJY_LOCAL_PR = MXGETPR(SWJY_LOCAL)
      SWJD_LOCAL_PR = MXGETPR(SWJD_LOCAL)

C     calculate NVD here? 
C     MWSIZE NVD
C     NVD = NV_IN * (NV_IN + 1) / 2
C
C     others?

C     do we need the BIN_IN_PR, THNEW_IN_PR, LWGTH_IN_PR, SWJY_IN_PR? in the interface?

C     do we need hakt as return value?

C     Call the actual FORTRAN subroutine

      CALL pvawsme(%val(Y_IN_PR), 
     +             %val(YD_IN_PR), 
     +             %val(POS_IN_PR), 
     +             NV_IN, 
     +             NVD_IN, 
     +             ND_IN, 
     +             N1_IN, 
     +             N2_IN, 
     +             N3_IN,
     +             HAKT_IN, 
     +             LAMBDA_IN, 
     +             %val(THETA_IN_PR), 
     +             %val(BI_IN_PR), 
     +             %val(BI_OUT_PR), 
     +             %val(THNEW_OUT_PR),
     +             %val(YDNEW_OUT_PR),
     +             %val(INVCOV_IN_PR), 
     +             NCORES_IN, 
     +             SPMIN_IN, 
     +             %val(LWGHT_LOCAL_PR), 
     +             %val(WGHT_IN_PR), 
     +             %val(SWJY_LOCAL_PR), 
     +             %val(SWJD_LOCAL_PR)
     +             NP1_IN, 
     +             NP2_IN, 
     +             NP3_IN)

      RETURN
      END 

      