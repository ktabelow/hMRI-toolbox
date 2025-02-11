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
C                         invcov,   12
C                         ncores,   13
C                         spmin,    14
C                         wght,     15
C                         np1,      16
C                         np2,      17
C                         np3);     18
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
      MWPOINTER BIN_OUT_PR
      MWPOINTER THNEW_OUT_PR
      MWPOINTER INVCOV_IN_PR
      MWPOINTER LWGHT_LOCAL_PR
      MWPOINTER WGHT_IN_PR
      MWPOINTER SWJY_LOCAL_PR

C     Check for proper number of arguments

      IF (NRHS .NE. 26) THEN
         CALL MEXERRMSGTXT('pvaws requires 18 input arguments')
      ELSEIF (NLHS .NE. 4) THEN
         CALL MEXERRMSGTXT('pvaws requires 2 output argument')
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
      INVCOV_IN_PR = MXGETPR(PRHS(12))
      NCORES_IN = INT(MXGETSCALAR(PRHS(13)))
      SPMIN_IN = INT(MXGETSCALAR(PRHS(14)))
      WGHT_IN_PR = MXGETPR(PRHS(15))
      NP1_IN = INT(MXGETSCALAR(PRHS(16)))
      NP2_IN = INT(MXGETSCALAR(PRHS(17)))
      NP3_IN = INT(MXGETSCALAR(PRHS(18)))

C     create all the PLHS
      
      PLHS(1) = MXCREATEDOUBLEMATRIX(1, NV_IN, 0)           CORRECT?
      PLHS(2) = MXCREATEDOUBLEMATRIX(NV_IN, NV_IN, 0)       CORRECT?

      BIN_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))

      dlw = 2 * floor(hakt ./ [1, wghts]) + 1;

      LWGHT_LOCAL = MXCREATEDOUBLEMATRIX(ONE, DWL, 0)       CORRECT? DWL UNDEFINED
      SWJY_LOCAL = MXCREATEDOUBLEMATRIX(NV, NCORES_IN, 0)   CORRECT?

      LWGHT_LOCAL_PR = MXGETPR(LWGHT_LOCAL)
      SWJY_LOCAL_PR = MXGETPR(SWJY_LOCAL)

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
     +            %val(BIN_OUT_PR), 
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
