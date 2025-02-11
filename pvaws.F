#include <fintrf.h>

C     Gateway function pvaws 
C     for FORTRAN subroutine pvaws2 in pcaws.f
C     
C     [bi, theta] = pvaws(y,         1
C                         pos,       2
C                         nv,        3
C                         nvd,       4
C                         n1,        5
C                         n2,        6
C                         n3,        7
C                         hakt,      8 
C                         lambda,    9
C                         theta,    10
C                         bi,       11
C                         bin,      12 create space in memory for PLHS
C                         thnew,    13 create space in memory for PLHS
C                         invcov,   14
C                         ncores,   15
C                         spmin,    16
C                         lwght,    17 create space in memory
C                         wght,     18
C                         swjy,     19 create space in memory
C                         np1,      20
C                         np2,      21
C                         np3);     22
C


      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)

      IMPLICIT NONE

C     variables for the interface

      MWPOINTER PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS

C     local integer variables

      MWSIZE NV_IN
      MWSIZE NVD_IN
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
      MWPOINTER POS_IN_PR
      MWPOINTER THETA_IN_PR
      MWPOINTER BI_IN_PR
      MWPOINTER BIN_IN_PR
      MWPOINTER THNEW_IN_PR
      MWPOINTER INVCOV_IN_PR
      MWPOINTER LWGHT_IN_PR
      MWPOINTER WGHT_IN_PR
      MWPOINTER SWJY_IN_PR

C     Check for proper number of arguments

      IF (NRHS .NE. 26) THEN
         CALL MEXERRMSGTXT('pvaws requires 22 input arguments')
      ELSEIF (NLHS .NE. 4) THEN
         CALL MEXERRMSGTXT('pvaws requires 3 output argument')
      ENDIF

C     get all the PRHS

      Y_IN_PR = MXGETPR(PRHS(1))
      POS_IN_PR = MXGETPR(PRHS(2))
      NV_IN = INT(MXGETSCALAR(PRHS(3)))
      NVD_IN = INT(MXGETSCALAR(PRHS(4)))
      N1_IN = INT(MXGETSCALAR(PRHS(5)))
      N2_IN = INT(MXGETSCALAR(PRHS(6)))
      N3_IN = INT(MXGETSCALAR(PRHS(7)))
      HAKT_IN = MXGETSCALAR(PRHS(8))
      LAMBDA_IN = MXGETSCALAR(PRHS(9))
      THETA_IN_PR = MXGETPR(PRHS(10))
      BI_IN_PR = MXGETPR(PRHS(11))
      BIN_IN_PR = MXGETPR(PRHS(12))          DO WE NEED THIS AT ALL??
      THNEW_IN_PR = MXGETPR(PRHS(13))        DO WE NEED THIS AT ALL??
      INVCOV_IN_PR = MXGETPR(PRHS(14))
      NCORES_IN = INT(MXGETSCALAR(PRHS(15)))
      SPMIN_IN = INT(MXGETSCALAR(PRHS(16)))
      LWGHT_IN_PR = MXGETPR(PRHS(17))        DO WE NEED THIS AT ALL??
      WGHT_IN_PR = MXGETPR(PRHS(18))
      SWJY_IN_PR = MXGETPR(PRHS(19))         DO WE NEED THIS AT ALL??
      NP1_IN = INT(MXGETSCALAR(PRHS(20)))
      NP2_IN = INT(MXGETSCALAR(PRHS(21)))
      NP3_IN = INT(MXGETSCALAR(PRHS(22)))

C     create all the PLHS
      
      PLHS(1) = MXCREATEDOUBLEMATRIX(1, NV_IN, 0)           CORRECT?
      PLHS(2) = MXCREATEDOUBLEMATRIX(NV_IN, NV_IN, 0)       CORRECT?

      BI_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))

      LWGHT_LOCAL = MXCREATEDOUBLEMATRIX(ONE, DWL, 0)       CORRECT? DWL UNDEFINED
      SWJY_LOCAL = MXCREATEDOUBLEMATRIX(NV, NCORES_IN, 0)   CORRECT?

      LWGHT_LOCAL_PR = MXGETPR(LWGHT_LOCAL)
      SWJY_LOCAL_PR = MXGETPR(SWJY_LOCAL)

C     calculate NVD here? 
C     MWSIZE NVD
C     NVD = NV_IN * (NV_IN + 1) / 2
C
C     others?

C     do we need the BIN_IN_PR, THNEW_IN_PR, LWGTH_IN_PR, SWJY_IN_PR? in the interface?

C     do we need hakt as return value?

C     Call the actual FORTRAN subroutine

      CALL PVAWS2(%val(Y_IN_PR), 
     +            %val(POS_IN_PR), 
     +            NV_IN, 
     +            NVD_IN, 
     +            N1_IN, 
     +            N2_IN, 
     +            N3_IN,
     +            HAKT_IN, 
     +            LAMBDA_IN, 
     +            %val(THETA_IN_PR), 
     +            %val(BI_IN_PR), 
     +            %val(BI_OUT_PR), 
     +            %val(THNEW_OUT_PR),
     +            %val(INVCOV_IN_PR), 
     +            NCORES_IN, 
     +            SPMIN_IN, 
     +            %val(LWGHT_LOCAL_PR), 
     +            %val(WGHT_IN_PR), 
     +            %val(SWJY_LOCAL_PR),
     +            NP1_IN, 
     +            NP2_IN, 
     +            NP3_IN)

      RETURN
      END 
