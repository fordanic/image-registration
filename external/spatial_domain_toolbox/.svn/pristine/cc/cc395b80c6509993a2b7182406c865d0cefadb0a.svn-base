#include "mex.h"

/*
 * P = MIN_QUADFORM(Qtot)
 * 
 * Minimize quadratic forms according to equations 6.20 -- 6.24 in
 * Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation and
 * Motion Estimation".
 *
 * Qtot   - A collection of quadratic forms, having the size
 *          HEIGTH x WIDTH x N x N.
 * 
 * P      - A collection of optimal parameters, having the size
 *          HEIGHT x WIDTH x (N-1).
 * 
 * Author: Gunnar Farnebäck
 *         Computer Vision Laboratory
 *         Linköping University, Sweden
 *         gf@isy.liu.se
 */

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, k;
    const int *Qtotdim;
    double *Qtot;
    mxArray *paramsarray;
    double *params;
    int paramsdim[3];
    mxArray *Qmatrix;
    double *Q;
    mxArray *qvector;
    double *q;
    double *p;
    mxArray *input[2];
    mxArray *output[1];
    int num_in, num_out;
    int N, M;
    
    /* Check the number of input and output arguments. */
    if (nrhs < 1)
	mexErrMsgTxt("Too few input arguments.");
    if (nrhs > 1)
	mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1)
	mexErrMsgTxt("Too many output arguments.");

    /* Check the formats of the input arguments. */
    for (i=0; i<nrhs; i++)
    {
	if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i])
	    || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
	{
	    mexErrMsgTxt("All input arguments must be real and full numeric arrays, stored as doubles.");
	}
    }

    if (mxGetNumberOfDimensions(prhs[0]) != 4)
	mexErrMsgTxt("Qtot must be 4D.");
    Qtotdim = mxGetDimensions(prhs[0]);

    if (Qtotdim[2] != Qtotdim[3])
	mexErrMsgTxt("Qtot must be HEIGHT x WIDTH x n x n.");
    M = Qtotdim[2] - 1;

    Qtot = mxGetPr(prhs[0]);
    
    /* Create the output array. */
    paramsdim[0] = Qtotdim[0];
    paramsdim[1] = Qtotdim[1];
    paramsdim[2] = M;
    paramsarray = mxCreateNumericArray(3, paramsdim, mxDOUBLE_CLASS, mxREAL);
    params = mxGetPr(paramsarray);

    Qmatrix = mxCreateDoubleMatrix(M, M, mxREAL);
    Q = mxGetPr(Qmatrix);
    input[0] = Qmatrix;
    qvector = mxCreateDoubleMatrix(M, 1, mxREAL);
    q = mxGetPr(qvector);
    input[1] = qvector;
    num_in = 2;
    num_out = 1;
    
    N = Qtotdim[0] * Qtotdim[1];
    for (k=0; k<N; k++)
    {
	for (i=0; i<M; i++)
	{
	    for (j=0; j<M; j++)
	    {
		Q[i+j*M] = Qtot[k+(i+j*(M+1))*N];
	    }
	    q[i] = -Qtot[k+(i+M*(M+1))*N];
	}
	mexCallMATLAB(num_out, output, num_in, input, "\\");
	p = mxGetPr(output[0]);
	for (i=0; i<M; i++)
	    params[k+i*N] = p[i];
	mxDestroyArray(output[0]);
    }
    
    /* Output the computed result. */
    plhs[0] = paramsarray;
}
