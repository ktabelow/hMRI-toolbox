#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// Forward declaration of the pvaws function from aws.c
extern void pvaws(
    double *y,     // 1
    int *pos,      // 2
    int nv,        // 3
    int nvd,       // 4
    int n1,        // 5
    int n2,        // 6
    int n3,        // 7
    double hakt,   // 8
    double lambda, // 9
    double *theta, // 10
    double *bi,    // 11
    double *bin,   // 12
    double *thnew, // 13
    double *invcov,// 14
    double spmin,  // 16
    double *lwght, // 17
    double *wght,  // 18
    double *swjy,  // 19
    int np1,       // 20
    int np2,       // 21
    int np3        // 22
);


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
    int *POS_IN_PR;
    double *THETA_IN_PR;
    double *BI_IN_PR;
    double *BI_OUT_PR;
    double *THNEW_OUT_PR;
    double *INVCOV_IN_PR;
    double *LWGHT_LOCAL_PR;
    double *WGHT_IN_PR;
    double *SWJY_LOCAL_PR;

    // Check for proper number of arguments

    if (nrhs != 18)
        mexErrMsgTxt("pvaws requires 18 input arguments");
    if (nlhs != 2)
        mexErrMsgTxt("pvaws requires 2 output arguments");

    // get all the prhs

    // Get input data pointers with dimension validation
    Y_IN_PR = mxGetPr(prhs[0]);
    if (mxGetNumberOfElements(prhs[0]) < N1_IN * N2_IN * N3_IN) {
        mexErrMsgTxt("Y input array size does not match expected dimensions");
    }

    POS_IN_PR = (int*)mxGetData(prhs[1]);
    if (mxGetNumberOfElements(prhs[1]) < N1_IN * N2_IN * N3_IN) {
        mexErrMsgTxt("POS input array size does not match expected dimensions (should be at least N1*N2*N3)");
    }
    // Get scalar parameters with type checking
    if (!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgTxt("NV_IN must be a scalar numeric value");
    }
    NV_IN = (mwSize)mxGetScalar(prhs[2]);

    if (!mxIsNumeric(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1) {
        mexErrMsgTxt("NVD_IN must be a scalar numeric value");
    }
    NVD_IN = (mwSize)mxGetScalar(prhs[3]);

    if (!mxIsNumeric(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgTxt("N1_IN must be a scalar numeric value");
    }
    N1_IN = (mwSize)mxGetScalar(prhs[4]);

    if (!mxIsNumeric(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgTxt("N2_IN must be a scalar numeric value");
    }
    N2_IN = (mwSize)mxGetScalar(prhs[5]);

    if (!mxIsNumeric(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgTxt("N3_IN must be a scalar numeric value");
    }
    N3_IN = (mwSize)mxGetScalar(prhs[6]);
    HAKT_IN = mxGetScalar(prhs[7]);
    LAMBDA_IN = mxGetScalar(prhs[8]);
    THETA_IN_PR = mxGetPr(prhs[9]);
    BI_IN_PR = mxGetPr(prhs[10]);
    INVCOV_IN_PR = mxGetPr(prhs[11]);

    if (!mxIsNumeric(prhs[12]) || mxGetNumberOfElements(prhs[12]) != 1) {
        mexErrMsgTxt("SPMIN_IN must be a scalar numeric value");
    }
    SPMIN_IN = (mwSize)mxGetScalar(prhs[12]);

    WGHT_IN_PR = mxGetPr(prhs[13]);

    if (!mxIsNumeric(prhs[14]) || mxGetNumberOfElements(prhs[14]) != 1) {
        mexErrMsgTxt("DLW_IN must be a scalar numeric value");
    }
    DLW_IN = (mwSize)mxGetScalar(prhs[14]);

    if (!mxIsNumeric(prhs[15]) || mxGetNumberOfElements(prhs[15]) != 1) {
        mexErrMsgTxt("NP1_IN must be a scalar numeric value");
    }
    NP1_IN = (mwSize)mxGetScalar(prhs[15]);

    if (!mxIsNumeric(prhs[16]) || mxGetNumberOfElements(prhs[16]) != 1) {
        mexErrMsgTxt("NP2_IN must be a scalar numeric value");
    }
    NP2_IN = (mwSize)mxGetScalar(prhs[16]);

    if (!mxIsNumeric(prhs[17]) || mxGetNumberOfElements(prhs[17]) != 1) {
        mexErrMsgTxt("NP3_IN must be a scalar numeric value");
    }
    NP3_IN = (mwSize)mxGetScalar(prhs[17]);

    // create all the plhs

    mwSize dims0[4] = {1, 1, 1, NV_IN};
    plhs[0] = mxCreateNumericArray(4, dims0, mxDOUBLE_CLASS, mxREAL);  // Changed to mxDOUBLE_CLASS for consistency
    mwSize dims1[3] = {1, 1, NV_IN};
    plhs[1] = mxCreateNumericArray(3, dims1, mxDOUBLE_CLASS, mxREAL);  // Changed to mxDOUBLE_CLASS for consistency

    if (plhs[0] == NULL || plhs[1] == NULL) {
        mexErrMsgTxt("Memory allocation failed for output arrays");
    }

    BI_OUT_PR = mxGetPr(plhs[0]);
    THNEW_OUT_PR = mxGetPr(plhs[1]);

    // create all arrays for internal calculations

    mxArray *LWGHT_LOCAL = mxCreateDoubleMatrix(1, DLW_IN, mxREAL);
    mxArray *SWJY_LOCAL = mxCreateDoubleMatrix(NV_IN, NCORES_IN, mxREAL);

    if (LWGHT_LOCAL == NULL || SWJY_LOCAL == NULL) {
        mxDestroyArray(LWGHT_LOCAL);
        mxDestroyArray(SWJY_LOCAL);
        mexErrMsgTxt("Memory allocation failed for temporary arrays");
    }

    LWGHT_LOCAL_PR = mxGetPr(LWGHT_LOCAL);
    SWJY_LOCAL_PR = mxGetPr(SWJY_LOCAL);

     // Call the actual C function 'pvaws2' which performs the main computation for the pvaws algorithm.
     // This function processes the input data and computes the outputs 'BI_OUT_PR' and 'THNEW_OUT_PR'.

     pvaws(
         Y_IN_PR,           // 1
         POS_IN_PR,         // 2
         (int)NV_IN,        // 3
         (int)NVD_IN,       // 4
         (int)N1_IN,        // 5
         (int)N2_IN,        // 6
         (int)N3_IN,        // 7
         HAKT_IN,           // 8
         LAMBDA_IN,         // 9
         THETA_IN_PR,       // 10
         BI_IN_PR,          // 11
         BI_OUT_PR,         // 12
         THNEW_OUT_PR,      // 13
         INVCOV_IN_PR,      // 14
         (double)SPMIN_IN,  // 15
         LWGHT_LOCAL_PR,    // 16
         WGHT_IN_PR,        // 17
         SWJY_LOCAL_PR,     // 18
         (int)NP1_IN,       // 19
         (int)NP2_IN,       // 20
         (int)NP3_IN        // 21
     );

    // Free allocated memory for temporary arrays used in internal calculations.
    // These arrays were created to store intermediate results and are no longer needed.
    mxDestroyArray(LWGHT_LOCAL);
    mxDestroyArray(SWJY_LOCAL);
}