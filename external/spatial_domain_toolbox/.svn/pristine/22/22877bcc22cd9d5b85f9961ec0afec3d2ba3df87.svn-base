#include "mex.h"
#include "math.h"

/* See prepare_displacement_matrices.m for documentation.
 *
 *  Author: Gunnar Farnebäck
 *          Computer Vision Laboratory
 *          Linköping University, Sweden
 *          gf@isy.liu.se
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int ndims;
    const int *dims;
    int height;
    int width;
    int result_dims[4];
    
    double *A1;
    double *b1;
    double *A2;
    double *b2;
    double *displacement;
    mxArray *A_array;
    mxArray *b_array;
    double *A;
    double *b;

    int i, j;
    int di, dj;
    int index, index2;
    int pixels;
    double aa1, aa2, aa3, aa4;
    
    /* Check the number of input and output arguments. */
    if (nrhs < 4)
	mexErrMsgTxt("Too few input arguments.");
    if (nrhs > 5)
	mexErrMsgTxt("Too many input arguments.");
    if (nlhs < 2)
	mexErrMsgTxt("Too few output arguments.");
    if (nlhs > 2)
	mexErrMsgTxt("Too many output arguments.");

    /* Check the formats of the input arguments. */
    for (i = 0; i < nrhs; i++)
    {
	if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i])
	    || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
	{
	    mexErrMsgTxt("The input arguments must be real and full numeric arrays, stored as doubles.");
	}
    }

    /* Get the dimensionalities of A1. Check size consistency for the
     * rest of the arguments.
     */
    ndims = mxGetNumberOfDimensions(prhs[0]);
    if (ndims != 4)
	mexErrMsgTxt("A1 must be MxNx2x2");
    dims = mxGetDimensions(prhs[0]);
    height = dims[0];
    width = dims[1];
    if (dims[2] != 2 || dims[3] != 2)
	mexErrMsgTxt("A1 must be MxNx2x2.");
	
    ndims = mxGetNumberOfDimensions(prhs[1]);
    if (ndims != 3)
	mexErrMsgTxt("b1 must be MxNx2");
    dims = mxGetDimensions(prhs[1]);
    if (dims[0] != height || dims[1] != width)
	mexErrMsgTxt("b1 and A1 have inconsistent sizes.");
    if (dims[2] != 2)
	mexErrMsgTxt("b1 must be MxNx2.");
	
    ndims = mxGetNumberOfDimensions(prhs[2]);
    if (ndims != 4)
	mexErrMsgTxt("A2 must be MxNx2x2");
    dims = mxGetDimensions(prhs[2]);
    if (dims[0] != height || dims[1] != width)
	mexErrMsgTxt("A2 and A1 have inconsistent sizes.");
    if (dims[2] != 2 || dims[3] != 2)
	mexErrMsgTxt("A2 must be MxNx2x2.");
	
    ndims = mxGetNumberOfDimensions(prhs[3]);
    if (ndims != 3)
	mexErrMsgTxt("b2 must be MxNx2");
    dims = mxGetDimensions(prhs[3]);
    if (dims[0] != height || dims[1] != width)
	mexErrMsgTxt("b2 and A1 have inconsistent sizes.");
    if (dims[2] != 2)
	mexErrMsgTxt("b2 must be MxNx2.");
	
    if (nrhs == 5)
    {
	ndims = mxGetNumberOfDimensions(prhs[4]);
	if (ndims != 3)
	    mexErrMsgTxt("displacement must be MxNx2");
	dims = mxGetDimensions(prhs[4]);
	if (dims[0] != height || dims[1] != width)
	    mexErrMsgTxt("displacement and A1 have inconsistent sizes.");
	if (dims[2] != 2)
	    mexErrMsgTxt("displacement must be MxNx2.");
    }

    /* Extract the double arrays from the arguments. */
    A1 = mxGetPr(prhs[0]);
    b1 = mxGetPr(prhs[1]);
    A2 = mxGetPr(prhs[2]);
    b2 = mxGetPr(prhs[3]);

    if (nrhs == 5)
	displacement = mxGetPr(prhs[4]);
    else
	displacement = NULL;

    /* Create the output arrays. */
    result_dims[0] = height;
    result_dims[1] = width;
    result_dims[2] = 2;
    result_dims[3] = 2;
    A_array = mxCreateNumericArray(4, result_dims, mxDOUBLE_CLASS, mxREAL);
    b_array = mxCreateNumericArray(3, result_dims, mxDOUBLE_CLASS, mxREAL);
    A = mxGetPr(A_array);
    b = mxGetPr(b_array);

    /* Number of pixels in the images. Used for indexing in third and
     * fourth dimensions.
     */
    pixels = height * width;
    
    /* Do the computations. */

    for (j = 0; j < width; j++)
	for (i = 0; i < height; i++)
	{
	    index = i + j * height;
	    /* Get displacement vector (di, dj). */
	    if (displacement)
	    {
		di = floor(0.5 + displacement[index]);
		if (i + di < 0)
		    di = -i;
		else if (i + di >= height)
		    di = height - i - 1;
		
		dj = floor(0.5 + displacement[index + pixels]);
		if (j + dj < 0)
		    dj = -j;
		else if (j + dj >= width)
		    dj = width - j - 1;

		index2 = index + di + dj * height;
	    }
	    else
	    {
		di = 0;
		dj = 0;
		index2 = index;
	    }


	    /* A(i,j,:,:) = (A1(i,j,:,:) + A2(i+di,j+dj,:,:)) / 2; */
    	    /* AA = squeeze(A(i,j,:,:)); */
	    aa1 = A[index] = (A1[index] + A2[index2])/2.0;
	    aa2 = A[index + pixels] = (A1[index + pixels]
				      + A2[index2 + pixels])/2.0;
	    aa3 = A[index + 2*pixels] = (A1[index + 2*pixels]
					+ A2[index2 + 2*pixels])/2.0;
	    aa4 = A[index + 3*pixels] = (A1[index + 3*pixels]
					+ A2[index2 + 3*pixels])/2.0;

    	    /* bb2 = squeeze(b2(i+di,j+dj,:)) - 2 * AA * [di;dj]; */
    	    /* b(i,j,:) = -(shiftdim(bb2,-2) - b1(i,j,:)) / 2; */

	    b[index] = -0.5 * (b2[index2] - b1[index]) + (aa1 * di + aa3 * dj);
	    b[index + pixels] = (-0.5 * (b2[index2 + pixels]
					 - b1[index + pixels])
				 + (aa2 * di + aa4 * dj));
	    
	}
    
    
    /* Output the computed result. */
    plhs[0] = A_array;
    plhs[1] = b_array;
}
