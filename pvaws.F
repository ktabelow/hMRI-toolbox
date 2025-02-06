#include <fintrf.h>

C     Gateway function for subroutine pvaws2 in pcaws.f

      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)

      IMPLICIT NONE

      MWPOINTER PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS

C     Check for proper number of arguments

      IF (NRHS .NE. 26) THEN
         CALL MEXERRMSGTXT('pvaws requires 22 input arguments')
      ELSEIF (NLHS .NE. 4) THEN
         CALL MEXERRMSGTXT('pvaws requires 3 output argument')
      ENDIF

C     Call the actual FORTRAN subroutine

      CALL pvaws2(y, pos, nv, nvd, n1, n2, n3,
     1            hakt, lambda, theta, bi, bin, thnew,
     2            invcov, ncores, spmin, lwght, wght, swjy,
     3            np1, np2, np3)

      RETURN
      END 
