#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//    Gateway function pvaws
//    for C function pvaws in pcaws.c
//
//     [bi, theta] = pvaws(y,         1
//                         pos,       2
//                         nv,        3
//                         nvd,       4
//                         n1,        5
//                         n2,        6
//                         n3,        7
//                         hakt,      8
//                         lambda,    9
//                         theta,    10
//                         bi,       11
//                         invcov,   12
//                         ncores,   13
//                         spmin,    14
//                         wght,     15
//                         dlw,      16
//                         np1,      17
//                         np2,      18
//                         np3);     19
//

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // local integer variables

    mwSize NV_IN;
    mwSize NVD_IN;
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
    double *POS_IN_PR;
    double *THETA_IN_PR;
    double *BI_IN_PR;
    double *BI_OUT_PR;
    double *THNEW_OUT_PR;
    double *INVCOV_IN_PR;
    double *LWGHT_LOCAL_PR;
    double *WGHT_IN_PR;
    double *SWJY_LOCAL_PR;

    // Check for proper number of arguments

    if (nrhs != 19)
        mexErrMsgTxt("pvaws requires 19 input arguments");
    if (nlhs != 2)
        mexErrMsgTxt("pvaws requires 2 output arguments");

    // get all the prhs

    Y_IN_PR = mxGetPr(prhs[0]);
    POS_IN_PR = mxGetPr(prhs[1]);
    NV_IN = (mwSize)mxGetScalar(prhs[2]);
    NVD_IN = (mwSize)mxGetScalar(prhs[3]);
    N1_IN = (mwSize)mxGetScalar(prhs[4]);
    N2_IN = (mwSize)mxGetScalar(prhs[5]);
    N3_IN = (mwSize)mxGetScalar(prhs[6]);
    HAKT_IN = mxGetScalar(prhs[7]);
    LAMBDA_IN = mxGetScalar(prhs[8]);
    THETA_IN_PR = mxGetPr(prhs[9]);
    BI_IN_PR = mxGetPr(prhs[10]);
    INVCOV_IN_PR = mxGetPr(prhs[11]);
    NCORES_IN = (mwSize)mxGetScalar(prhs[12]);
    SPMIN_IN = (mwSize)mxGetScalar(prhs[13]);
    WGHT_IN_PR = mxGetPr(prhs[14]);
    DLW_IN = (mwSize)mxGetScalar(prhs[15]);
    NP1_IN = (mwSize)mxGetScalar(prhs[16]);
    NP2_IN = (mwSize)mxGetScalar(prhs[17]);
    NP3_IN = (mwSize)mxGetScalar(prhs[18]);

    // create all the plhs

    plhs[0] = mxCreateNumericArray(4, NV_IN, mxSINGLE_CLASS, mxREAL);  // CORRECT ? 
    plhs[1] = mxCreateNumericArray(3, NV_IN, mxSINGLE_CLASS, mxREAL);  // CORRECT ?
 
    BI_OUT_PR = mxGetPr(plhs[0]);
    THNEW_OUT_PR = mxGetPr(plhs[1]);

    // create all arrays for internal calculations

    mxArray *LWGHT_LOCAL = mxCreateDoubleMatrix(1, DLW_IN, mxREAL);
    mxArray *SWJY_LOCAL = mxCreateDoubleMatrix(NV_IN, NCORES_IN, mxREAL);

    LWGHT_LOCAL_PR = mxGetPr(LWGHT_LOCAL);
    SWJY_LOCAL_PR = mxGetPr(SWJY_LOCAL);

    // Call the actual C function 'pvaws2' which performs the main computation for the pvaws algorithm.
    // This function processes the input data and computes the outputs 'BI_OUT_PR' and 'THNEW_OUT_PR'.

    pvaws(
        Y_IN_PR,
        POS_IN_PR,
        NV_IN,
        NVD_IN,
        N1_IN,
        N2_IN,
        N3_IN,
        HAKT_IN,
        LAMBDA_IN,
        THETA_IN_PR,
        BI_IN_PR,
        BI_OUT_PR,
        THNEW_OUT_PR,
        INVCOV_IN_PR,
        NCORES_IN,
        SPMIN_IN,
        LWGHT_LOCAL_PR,
        WGHT_IN_PR,
        SWJY_LOCAL_PR,
        NP1_IN,
        NP2_IN,
        NP3_IN
    );

    // Free allocated memory for temporary arrays used in internal calculations.
    // These arrays were created to store intermediate results and are no longer needed.
    mxDestroyArray(LWGHT_LOCAL);
    mxDestroyArray(SWJY_LOCAL);
}