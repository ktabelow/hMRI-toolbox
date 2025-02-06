#include <fintrf.h>

C     Gateway function for subroutine pvawsme in pcaws.f

      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)

      IMPLICIT NONE

      MWPOINTER PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS

C     Check for proper number of arguments

      IF (NRHS .NE. 26) THEN
         CALL MEXERRMSGTXT('pvawslast requires 26 input arguments')
      ELSEIF (NLHS .NE. 4) THEN
         CALL MEXERRMSGTXT('pvawslast requires 4 output argument')
      ENDIF

C     Call the actual FORTRAN subroutine

      CALL pvawsme(y, yd, pos, nv, nvd, nd, n1, n2, n3,
     1             hakt, lambda, theta, bi, bin, thnew, ydnew,
     2             invcov, ncores, spmin, lwght, wght, swjy, swjd,
     3             np1, np2, np3)

      RETURN
      END 

      