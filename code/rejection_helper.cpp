#include "mex.h"
#include <random>
#include <cmath>

extern "C" {
    void dscal(const mwSize* n, const double* alpha, double* x, const mwSize* incx);
    void dsyr(const char* uplo, const mwSize* n, const double* alpha, const double* x, const mwSize* incx, double* A, const mwSize* lda);
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("rejection_helper:inputError", "Two inputs required: matrix A and vector w.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("rejection_helper:outputError", "Two outputs required.");
    }

    /* Ensure the first input is a double matrix */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("rejection_helper:inputError", "First input must be a non-complex double matrix.");
    }

    /* Ensure the second input is a double vector */
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("rejection_helper:inputError", "Second input must be a non-complex double vector.");
    }

    /* Get the dimensions of matrix A */
    mwSize n = mxGetN(prhs[0]);  // assuming A is n by n
    if (n != mxGetM(prhs[0])) {
        mexErrMsgIdAndTxt("rejection_helper:inputError", "Matrix A must be square (n x n).");
    }

    /* Get the length of vector w */
    mwSize w_len = mxGetM(prhs[1]);
    if (w_len != n) {
        mexErrMsgIdAndTxt("rejection_helper:inputError", "Vector w must be of length n.");
    }

    /* Get pointers to the data */
    double *A = mxGetPr(prhs[0]);  // pointer to matrix A
    double *w = mxGetPr(prhs[1]);  // pointer to vector w

    /* Create the output vector */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);  // output vector of length n
    double *out = mxGetPr(plhs[0]);  // pointer to output data
    plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);  // output matrix of size n x n
    double *A_copy = mxGetPr(plhs[1]);  // pointer to output matrix (modified A)

    /* Create a copy of A to avoid modifying the original in place */
    std::copy(A, A + n * n, A_copy);

    int idx = 0;

    /* Random number generation */
    std::random_device rd;   // Seed
    std::mt19937 gen(rd());  // Mersenne Twister engine
    std::uniform_real_distribution<> dis(0.0, 1.0);  // Range [0, 1)

    /* Modify A in place and compute the output */
    for (mwSize i = 0; i < n; i++) {
        // mexPrintf("%f\n", A_copy[i + i*n]);  // Debugging line (can be removed)
        if (w[i] * dis(gen) >= A_copy[i + i*n]) continue; /* Rejection sampling */
        out[idx] = i+1; idx++;
        double scale = 1.0 / std::sqrt(A_copy[i + i*n]);
        mwSize inc = 1;  // Increment

        /* Scale A_copy[i:end,i] by 1/sqrtAii */
        mwSize len = n - i;  // Number of elements in the column from i to end
        dscal(&len, &scale, &A_copy[i + i*n], &inc);

        /* Modify A_copy[i+1:end,i+1:end] -= A_copy[i+1:end,i]*A_copy[i+1:end,i].T */
        char uplo = 'L';  // Only update the lower triangular part
        double alpha = -1.0;
        len = n - i - 1;
        dsyr(&uplo, &len, &alpha, &A_copy[i+1+i*n], &inc, &A_copy[(i + 1) + (i + 1)*n], &n);
    }
}