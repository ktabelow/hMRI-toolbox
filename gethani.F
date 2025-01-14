#include <fintrf.h>

C
C Gateway function for gethani.F
C  to be called:
C  bw = gethani(x,y,kern,value,wght,eps)
C
C   in R
C
C   gethani <- function(x,y,lkern,value,wght,eps=1e-2){
C  .Fortran("gethani",
C           as.double(x),
C           as.double(y),
C           as.integer(lkern),
C           as.double(value),
C           as.double(wght),
C           as.double(eps),
C           bw=double(1),
C           PACKAGE="qMRI")$bw
C}

C 
C
C 
C 
C
        SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     
C
       IMPLICIT NONE
       MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

       INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     
C
       MWPOINTER MXGETPR
       MWPOINTER MXCREATEDOUBLEMATRIX, mxCreateDoubleScalar
       MWPOINTER MXCREATENUMERICMATRIX
       MWPOINTER MXCREATENUMERICARRAY
       
       mwPointer mxGetNumberOfElements
       mwPointer mxGetDimensions
       real*8  mxGetScalar
       integer*4 mxIsNumeric, mxClassIDFromClassName
       character*8 double

C-----------------------------------------------------------------------
C

       MWSIZE MXGETM, MXGETN
       mwSize mxGetNumberOfDimensions
C

C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     
C
C       MWPOINTER X_IN_PR, X_OUT_PR, Y_IN_PR, Y_OUT_PR
       MWPOINTER BW_IN_PR, BW_OUT_PR
       MWPOINTER WGHT_IN_PR

      

C      POINTERS TO M (rows) AND N (columns) of INPUTS

       MWPOINTER M_X, N_X, M_Y, N_Y, M_BW, N_BW,
     + M_WGHT, N_WGHT, M_KERN, N_KERN, M_VALUE, N_VALUE,
     + M_EPS, N_EPS
     
       
C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES

       INTEGER*4 KERN_IN       
       
        REAL*8 X_IN, Y_IN, BW_INOUT, VALUE_IN, EPS_IN, WGHT_IN(2)  


       
C      SIZE OF THE INPUT ARGUMENTS (MXN)      

       MWSIZE DIM_X, DIM_Y, DIM_BW, DIM_WGHT, DIM_KERN, 
     + TWO, ONE, DIM_EPS, DIM_VALUE
     
       

       
       MWSIZE PRODN 

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 6) THEN
         CALL MEXERRMSGTXT('GETHANI requires 6 input arguments')
       ELSEIF (NLHS .NE. 1) THEN
         CALL MEXERRMSGTXT('GETHANI requires 1 output argument')
       ENDIF
      
C      CHECK IF THE IMPUTS ARE NUMERIC

C       CALL MEXPRINTF(' BEFORE Checking numerics ')

       IF(mxIsNumeric(prhs(1)) .eq. 0) then 
           CALL MEXERRMSGTXT('X HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(2)) .eq. 0) then 
           CALL MEXERRMSGTXT('Y HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(3)) .eq. 0) then 
           CALL MEXERRMSGTXT('KERN HAS TO BE A NUMBER')
        ENDIF
       
        IF(mxIsNumeric(prhs(4)) .eq. 0) then 
           CALL MEXERRMSGTXT('VALUE HAS TO BE A NUMBER')
        ENDIF
  
        IF(mxIsNumeric(prhs(5)) .eq. 0) then 
           CALL MEXERRMSGTXT('WGHT HAS TO BE A VECTOR')
        ENDIF

         IF(mxIsNumeric(prhs(6)) .eq. 0) then 
           CALL MEXERRMSGTXT('EPS HAS TO BE A NUMBER')
        ENDIF

C        IF(mxIsNumeric(prhs(7)) .eq. 0) then 
C           CALL MEXERRMSGTXT('BW HAS TO BE A NUMBER')
C        ENDIF

       

        
       
C        CALL MEXPRINTF(' after checking numerics ')
C
C CHECK THE DIMENSIONS OF INPUT VALUES.  
C N1,N2,N3,... HAVE TO BE 1
C 
C 
C
       M_X = MXGETM(PRHS(1))
       N_X = MXGETN(PRHS(1))
       M_Y = MXGETM(PRHS(2))
       N_Y = MXGETN(PRHS(2))
C       M_BW = MXGETM(PRHS(7))
C       N_BW = MXGETN(PRHS(7))
       M_WGHT = MXGETM(PRHS(5))
       N_WGHT = MXGETN(PRHS(5))
      
    
       M_KERN = MXGETM(PRHS(3))
       N_KERN = MXGETN(PRHS(3))
       M_VALUE = MXGETM(PRHS(4))
       N_VALUE = MXGETN(PRHS(4))
       M_EPS = MXGETM(PRHS(6))
       N_EPS = MXGETN(PRHS(6))
      
       

      
C       NARRAY = MXGETDIMENSIONS(PRHS(1))
       
       
       
       DIM_X = M_X*N_X
       DIM_Y = M_Y*N_Y
C       DIM_BW = M_BW*N_BW       
       DIM_WGHT = M_WGHT*N_WGHT
       DIM_VALUE = M_VALUE*N_VALUE
       DIM_KERN = M_KERN*N_KERN
       DIM_EPS = M_EPS*N_EPS
     
C       CALL MEXPRINTF(' checking dimensions ')
C
      
      IF ((DIM_X .NE. 1).OR.(DIM_Y .NE. 1)) THEN
        CALL MEXERRMSGTXT('X, Y, BW have to be scalar')
      ENDIF

      IF (DIM_KERN .NE. 1) THEN
        CALL MEXERRMSGTXT('KERN has to be a SCALAR')
      ENDIF
     
      IF (DIM_VALUE .NE. 1) THEN
        CALL MEXERRMSGTXT('VALUE has to be a SCALAR')
      ENDIF

      IF (DIM_EPS .NE. 1) THEN
        CALL MEXERRMSGTXT('EPS has to be a SCALAR')
      ENDIF

      IF (DIM_WGHT .NE. 2) THEN
        CALL MEXERRMSGTXT('WGHT has to be a VECTOR OF LENGHT 2')
      ENDIF


C        CALL MEXPRINTF(' BEFORE getting SCALARS ')


        KERN_IN = INT(MXGETSCALAR(PRHS(3)))
        X_IN = MXGETSCALAR(PRHS(1))
        Y_IN = MXGETSCALAR(PRHS(2))
       
        
         VALUE_IN = MXGETSCALAR(PRHS(4))
         EPS_IN = MXGETSCALAR(PRHS(6))

C       CALL MEXPRINTF(' AFTER getting SCALARS ')

       
C       TWO= MXGETNUMBEROFELEMENTS(PRHS(5))
       ONE=1

C       getting all the other INPUT pointers
       

C       BW_IN_PR = MXGETPR(PRHS(7))
       WGHT_IN_PR = MXGETPR(PRHS(5))
       
       
     
C       CALL MEXPRINTF(' AFTER getting OTHER input ')

C       CALL MEXPRINTF(' BEFORE CREATING OUTPUT ')
C CREATE MATRICES FOR RETURN ARGUMENT
C

      
      PLHS(1) = mxCreateDoubleScalar(1)
      

      
C
C ASSIGN POINTERS TO THE OUTPUT PARAMETERS
C
      BW_OUT_PR = MXGETPR(PLHS(1))
      
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
     


C       CALL MXCOPYPTRTOREAL8(BW_IN_PR, BW_INOUT,ONE)
   
      CALL MXCOPYPTRTOREAL8(WGHT_IN_PR, WGHT_IN, DIM_WGHT)
     
C      CALL MEXPRINTF('BEFORE CALL')
C
C DO THE ACTUAL COMPUTATIONS 
C

        CALL cgethani( X_IN,Y_IN,KERN_IN,VALUE_IN,
     + WGHT_IN,EPS_IN,BW_INOUT)
C      CALL MEXPRINTF(' after call ')
      
C    
      CALL MXCOPYREAL8TOPTR(BW_INOUT,BW_OUT_PR,ONE)


      
     

C         CALL MEXPRINTF(' time to return ')

      RETURN
      END
