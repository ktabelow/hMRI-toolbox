#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//     Gateway function pvawslast
//     for C function pvawsme in pcaws.c
//     
//     [bi, theta, smoothedmpmData] = pvawslast(y,       1
//                                              yd,      2
//                                              pos,     3
//                                              nv,      4
//                                              nvd,     5
//                                              nd,      6
//                                              n1,      7
//                                              n2,      8
//                                              n3,      9
//                                              hakt,   10 
//                                              lambda, 11
//                                              theta,  12
//                                              bi,     13
//                                              invcov, 14
//                                              ncores, 15
//                                              spmin,  16
//                                              wght,   17
//                                              dlw,    18
//                                              np1,    19
//                                              np2,    20
//                                              np3);   21

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

     // local integer variables

      mwSize NV_IN;
      mwSize NVD_IN;
      mwSize ND_IN;
      mwSize N1_IN;
      mwSize N2_IN;
      mwSize N3_IN; 
      mwSize NCORES_IN;
      mwSize SPMIN_IN; 
      mwSize DLW_IN;
      mwSize NP1_IN;
      mwSize NP2_IN;
      mwSize NP3_IN;

      // local real variables (non-portable syntax)

      double HAKT_IN;
      double LAMBDA_IN;

      // local input pointers

      double *Y_IN_PR;
      double *YD_IN_PR;
      double *POS_IN_PR;
      double *THETA_IN_PR;
      double *BI_IN_PR;
      double *BI_OUT_PR;
      double *THNEW_OUT_PR;
      double *YDNEW_OUT_PR;
      double *INVCOV_IN_PR;
      double *LWGHT_LOCAL_PR;
      double *WGHT_IN_PR;
      double *SWJY_LOCAL_PR;
      double *SWJD_LOCAL_PR;

      // Check for proper number of arguments

      if (nrhs != 21)
         mexErrMsgTxt("pvawslast requires 21 input arguments");
      if (nlhs != 4)
         mexErrMsgTxt("pvawslast requires 4 output arguments");

      // Validate output arguments
      if (mxGetM(plhs[0]) != 1 || mxGetN(plhs[0]) != NV_IN)
         mexErrMsgTxt("Output argument 1 must be a 1xNV_IN matrix.");
      if (mxGetM(plhs[1]) != NV_IN || mxGetN(plhs[1]) != NV_IN)
         mexErrMsgTxt("Output argument 2 must be a NV_INxNV_IN matrix.");
      if (mxGetM(plhs[2]) != NV_IN || mxGetN(plhs[2]) != NCORES_IN)
         mexErrMsgTxt("Output argument 3 must be a NV_INxNCORES_IN matrix.");
      if (!mxIsDouble(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2)
         mexErrMsgTxt("Input y must be a 2D double array.");
      if (!mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2)
         mexErrMsgTxt("Input yd must be a 2D double array.");
      if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) != 2)
         mexErrMsgTxt("Input pos must be a 2D double array.");
      if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
         mexErrMsgTxt("Input nv must be a scalar.");
      if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1)
         mexErrMsgTxt("Input nvd must be a scalar.");
      if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1)
         mexErrMsgTxt("Input nd must be a scalar.");
      if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1)
         mexErrMsgTxt("Input n1 must be a scalar.");
      if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1)
         mexErrMsgTxt("Input n2 must be a scalar.");
      if (!mxIsDouble(prhs[8]) || mxGetNumberOfElements(prhs[8]) != 1)
         mexErrMsgTxt("Input n3 must be a scalar.");
      if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1)
         mexErrMsgTxt("Input hakt must be a scalar.");
      if (!mxIsDouble(prhs[10]) || mxGetNumberOfElements(prhs[10]) != 1)
         mexErrMsgTxt("Input lambda must be a scalar.");
      if (!mxIsDouble(prhs[11]) || mxGetNumberOfDimensions(prhs[11]) != 2)
         mexErrMsgTxt("Input theta must be a 2D double array.");
      if (!mxIsDouble(prhs[12]) || mxGetNumberOfDimensions(prhs[12]) != 2)
         mexErrMsgTxt("Input bi must be a 2D double array.");
      if (!mxIsDouble(prhs[13]) || mxGetNumberOfDimensions(prhs[13]) != 2)
         mexErrMsgTxt("Input invcov must be a 2D double array.");
      if (!mxIsDouble(prhs[14]) || mxGetNumberOfElements(prhs[14]) != 1)
         mexErrMsgTxt("Input ncores must be a scalar.");
      if (!mxIsDouble(prhs[15]) || mxGetNumberOfElements(prhs[15]) != 1)
         mexErrMsgTxt("Input spmin must be a scalar.");
      if (!mxIsDouble(prhs[16]) || mxGetNumberOfDimensions(prhs[16]) != 2)
         mexErrMsgTxt("Input wght must be a 2D double array.");
      if (!mxIsDouble(prhs[17]) || mxGetNumberOfElements(prhs[17]) != 1)
         mexErrMsgTxt("Input dlw must be a scalar.");
      if (!mxIsDouble(prhs[18]) || mxGetNumberOfElements(prhs[18]) != 1)
         mexErrMsgTxt("Input np1 must be a scalar.");
      if (!mxIsDouble(prhs[19]) || mxGetNumberOfElements(prhs[19]) != 1)
         mexErrMsgTxt("Input np2 must be a scalar.");
      if (!mxIsDouble(prhs[20]) || mxGetNumberOfElements(prhs[20]) != 1)
         mexErrMsgTxt("Input np3 must be a scalar.");
      if (nlhs != 4)
         mexErrMsgTxt("pvawslast requires 4 output arguments");


      // get all the prhs 

      Y_IN_PR = mxGetPr(prhs[0]);
      YD_IN_PR = mxGetPr(prhs[1]);
      POS_IN_PR = mxGetPr(prhs[2]);
      NV_IN = (mwSize)mxGetScalar(prhs[3]);
      NVD_IN = (mwSize)mxGetScalar(prhs[4]);
      ND_IN = (mwSize)mxGetScalar(prhs[5]);
      N1_IN = (mwSize)mxGetScalar(prhs[6]);
      N2_IN = (mwSize)mxGetScalar(prhs[7]);
      N3_IN = (mwSize)mxGetScalar(prhs[8]);
      HAKT_IN = mxGetScalar(prhs[9]);
      LAMBDA_IN = mxGetScalar(prhs[10]);
      THETA_IN_PR = mxGetPr(prhs[11]);
      BI_IN_PR = mxGetPr(prhs[12]);
      INVCOV_IN_PR = mxGetPr(prhs[13]);
      NCORES_IN = (mwSize)mxGetScalar(prhs[14]);
      SPMIN_IN = (mwSize)mxGetScalar(prhs[15]);
      WGHT_IN_PR = mxGetPr(prhs[16]);
      DLW_IN = (mwSize)mxGetScalar(prhs[17]);
      NP1_IN = (mwSize)mxGetScalar(prhs[18]);
      NP2_IN = (mwSize)mxGetScalar(prhs[19]);
      NP3_IN = (mwSize)mxGetScalar(prhs[20]);

      // Create all the plhs (output arrays)

      plhs[0] = mxCreateDoubleMatrix(1, NV_IN, mxREAL);           // Correct?
      plhs[1] = mxCreateDoubleMatrix(NV_IN, NV_IN, mxREAL);       // Correct?
      plhs[2] = mxCreateDoubleMatrix(NV_IN, NCORES_IN, mxREAL);   // Correct?
      
      BI_OUT_PR = mxGetPr(plhs[0]);
      THNEW_OUT_PR = mxGetPr(plhs[1]);
      YDNEW_OUT_PR = mxGetPr(plhs[2]);

      // create all arrays for internal calculations 

      mxArray *LWGHT_LOCAL = mxCreateDoubleMatrix(1, DLW_IN, mxREAL);       // Correct?
      mxArray *SWJY_LOCAL = mxCreateDoubleMatrix(NV_IN, NCORES_IN, mxREAL); // Correct?
      mxArray *SWJD_LOCAL = mxCreateDoubleMatrix(NV_IN, NCORES_IN, mxREAL); // Correct?

      LWGHT_LOCAL_PR = mxGetPr(LWGHT_LOCAL);
      SWJY_LOCAL_PR = mxGetPr(SWJY_LOCAL);
      SWJD_LOCAL_PR = mxGetPr(SWJD_LOCAL);

      // Call the actual C function
      // pvawsme is a computational function that performs smoothing and statistical analysis
      // on the input data. It takes various inputs such as data arrays, scalar parameters, 
      // and matrices, and produces smoothed outputs and updated statistical parameters.
      // Inputs:
      // - Y_IN_PR, YD_IN_PR, POS_IN_PR: Input data arrays
      // - NV_IN, NVD_IN, ND_IN, N1_IN, N2_IN, N3_IN: Scalar parameters for dimensions
      // - HAKT_IN, LAMBDA_IN: Scalar parameters for smoothing
      // - THETA_IN_PR, BI_IN_PR, INVCOV_IN_PR: Input matrices for statistical calculations
      // - NCORES_IN, SPMIN_IN: Scalar parameters for computational settings
      // - WGHT_IN_PR: Weight matrix
      // - NP1_IN, NP2_IN, NP3_IN: Scalar parameters for additional dimensions
      // Outputs:
      // - BI_OUT_PR: Smoothed output data
      // - THNEW_OUT_PR: Updated statistical parameters
      // - YDNEW_OUT_PR: Smoothed derivative data
      // - Internal arrays (LWGHT_LOCAL_PR, SWJY_LOCAL_PR, SWJD_LOCAL_PR) are used for intermediate calculations.

      pvawslast(Y_IN_PR, 
                YD_IN_PR, 
                POS_IN_PR, 
                NV_IN, 
                NVD_IN, 
                ND_IN, 
                N1_IN, 
                N2_IN, 
                N3_IN,
                HAKT_IN, 
                LAMBDA_IN, 
                THETA_IN_PR, 
                BI_IN_PR, 
                BI_OUT_PR, 
                THNEW_OUT_PR,
                YDNEW_OUT_PR,
                INVCOV_IN_PR, 
                NCORES_IN, 
                SPMIN_IN, 
                LWGHT_LOCAL_PR, 
                WGHT_IN_PR, 
                SWJY_LOCAL_PR, 
                SWJD_LOCAL_PR,
                NP1_IN, 
                NP2_IN, 
                NP3_IN);

// Free allocated memory for internal calculations
mxDestroyArray(LWGHT_LOCAL);
mxDestroyArray(SWJY_LOCAL);
mxDestroyArray(SWJD_LOCAL);
}

      