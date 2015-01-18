#include "mex.h"

#define DEBUG 0
#define DIMENSIONALITY 3

/*
 *  conv3_readable.c
 *
 *  Perform convolution in 3 dimensions. Only real images and kernels
 *  are supported. An optimized version of this code is available as
 *  conv3.c. 
 *
 *  Call:
 *  result=conv3(signal,kernel)
 *  or
 *  result=conv3(signal,kernel,region_of_interest)
 *
 *  signal             - signal values
 *  kernel             - filter kernel
 *  region_of_interest - where to compute the convolution
 *  result             - filter output
 *
 *  Formats:
 *  signal is a 3-dimensional array.
 *  kernel is a 3-dimensional array.
 *  region of interest is a Nx2 matrix, where each row gives start
 *          and end indices for each signal dimension. N must be the
 *          same as the number of non-trailing singleton dimensions
 *          for signal.
 *  result is a 3-dimensional array. The size is the same as for
 *          signal, unless region_of_interest is specified. The
 *          signal is assumed to be surrounded by zeros.
 *
 *  Note 1: Only double, nonsparse arrays are currently supported.
 *  Note 2: If signal or kernel has fewer than 3 dimensions,
 *          trailing singleton dimensions are added.
 *
 *
 *  Author: Gunnar Farnebäck
 *          Computer Vision Laboratory
 *          Linköping University, Sweden
 *          gf@isy.liu.se
 */

static void
conv3(const double *signal,
      const double *kernel,
      double *result,
      const int *signal_dimensions,
      const int *kernel_dimensions,
      const int *result_dimensions,
      const int *start_indices,
      const int *stop_indices)
{
    int result_indices[DIMENSIONALITY];
    int kernel_indices[DIMENSIONALITY];
    int signal_indices[DIMENSIONALITY];
    int displacements[DIMENSIONALITY];
    
    int signalindex;
    int kernelindex;
    int resultindex;
    
    double ip;

    int i;

    /* The center of the kernel is supposed to be the local origin.
       Compute the corresponding offsets. */ 
    for (i = 0; i < DIMENSIONALITY; i++)
    {
	/* Note that this is integer division. */
	displacements[i] = (kernel_dimensions[i] - 1) / 2;
    }

    /* Loop over the signal dimensions */
    for (result_indices[2] = start_indices[2];
	 result_indices[2] <= stop_indices[2];
	 result_indices[2]++)
    {
	for (result_indices[1] = start_indices[1];
	     result_indices[1] <= stop_indices[1];
	     result_indices[1]++)
	{
	    for (result_indices[0] = start_indices[0];
		 result_indices[0] <= stop_indices[0];
		 result_indices[0]++)
	    {
		/* Reset the accumulated inner product. */
		ip = 0.0;
		
		/* Loop over the dimensions for the kernel. */
		for (kernel_indices[2] = 0;
		    kernel_indices[2] < kernel_dimensions[2];
		    kernel_indices[2]++)
		{
		    /* Compute the signal index corresponding
		       to the current result index and kernel
		       index. */
		    signal_indices[2] = (result_indices[2]
					 + kernel_indices[2]
					 - displacements[2]);
		    /* Check if we are outside the signal
		       boundary. It is implied that the
		       signal is zero then. */
		    if (signal_indices[2] < 0
			|| signal_indices[2] >= signal_dimensions[2])
		    {
			continue;
		    }
		    for (kernel_indices[1] = 0;
			 kernel_indices[1] < kernel_dimensions[1];
			 kernel_indices[1]++)
		    {
			signal_indices[1] = (result_indices[1]
					     + kernel_indices[1]
					     - displacements[1]);
			if (signal_indices[1] < 0
			    || signal_indices[1] >= signal_dimensions[1])
			{
			    continue;
			}
			for (kernel_indices[0] = 0;
			     kernel_indices[0] < kernel_dimensions[0];
			     kernel_indices[0]++)
			{
			    signal_indices[0] = (result_indices[0]
						 + kernel_indices[0]
						 - displacements[0]);
			    if (signal_indices[0] < 0
				|| signal_indices[0] >= signal_dimensions[0])
			    {
				continue;
			    }

			    signalindex = ((signal_indices[2]
					    * signal_dimensions[1]
					    + signal_indices[1])
					   * signal_dimensions[0]
					   + signal_indices[0]);
			    kernelindex = ((kernel_indices[2]
					    * kernel_dimensions[1]
					    + kernel_indices[1])
					   * kernel_dimensions[0]
					   + kernel_indices[0]);
			    ip += signal[signalindex] * kernel[kernelindex];
			}
		    }
		}
		resultindex = (((result_indices[2]-start_indices[2])
				* result_dimensions[1]
				+ (result_indices[1]-start_indices[1]))
			       * result_dimensions[0]
			       + (result_indices[0] - start_indices[0]));
		result[resultindex] = ip;
	    }
	}
    }
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i;
    int sn, kn;
    const int *sdim;
    const int *kdim;
    int signal_dimensions[DIMENSIONALITY];
    int kernel_dimensions[DIMENSIONALITY];
    int result_dimensions[DIMENSIONALITY];
    mxArray *resultarray;
    double *region_of_interest;
    int startindex, stopindex;
    int start_indices[DIMENSIONALITY];
    int stop_indices[DIMENSIONALITY];
    int signal_is_complex = 0;
    
    /* Check the number of input and output arguments. */
    if (nrhs < 2)
	mexErrMsgTxt("Too few input arguments.");
    if (nrhs > 3)
	mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1)
	mexErrMsgTxt("Too many output arguments.");

    /* Check the formats of the input arguments. */
    for (i = 0; i < nrhs; i++)
    {
	if (!mxIsNumeric(prhs[i]) || mxIsSparse(prhs[i])
	    || !mxIsDouble(prhs[i]))
	{
	    mexErrMsgTxt("All input arguments must be full numeric arrays, stored as doubles.");
	}
	if (mxIsComplex(prhs[i]))
	{
	    if (i == 0)
		signal_is_complex = 1;
	    else
		mexErrMsgTxt("Only the signal may be complex.");
	}
    }

    /* We won't deal with empty arrays. */
    for (i = 0; i < nrhs; i++)
	if (mxIsEmpty(prhs[i]))
	    mexErrMsgTxt("Empty arrays not allowed.");

    /* Get the dimensionalities. Trailing singleton dimensions are ignored. */
    /* signal */
    sdim = mxGetDimensions(prhs[0]);
    for (sn = mxGetNumberOfDimensions(prhs[0]); sn > 0; sn--)
	if (sdim[sn-1] > 1)
	    break;
    /* kernel */
    kdim = mxGetDimensions(prhs[1]);
    for (kn = mxGetNumberOfDimensions(prhs[1]); kn > 0; kn--)
	if (kdim[kn-1] > 1)
	    break;

    /* Verify that there are no more than 3 dimensions. */
    if (sn > DIMENSIONALITY)
	mexErrMsgTxt("More than 3-dimensional signal.");
    if (kn > DIMENSIONALITY)
	mexErrMsgTxt("More than 3-dimensional kernel.");

    /* Build dimension vectors for signal and kernel. */
    for (i = 0; i < sn; i++)
	signal_dimensions[i] = sdim[i];
    for (; i < DIMENSIONALITY; i++)
	signal_dimensions[i] = 1;

    for (i = 0; i < kn; i++)
	kernel_dimensions[i] = kdim[i];
    for (; i < DIMENSIONALITY; i++)
	kernel_dimensions[i] = 1;

    if (nrhs == 3)
    {
	if (mxGetNumberOfDimensions(prhs[2]) != 2
	    || mxGetM(prhs[2]) != sn || mxGetN(prhs[2]) != 2)
	{
	    mexErrMsgTxt("Region of interest must be an N by 2 matrix, where N is the dimensionality of the signal.");
	}
	region_of_interest = mxGetPr(prhs[2]);
	for (i = 0; i < sn; i++)
	{
	    startindex = region_of_interest[i];
	    stopindex = region_of_interest[i+sn];
	    if (startindex < 1
		|| startindex > stopindex
		|| stopindex > sdim[i])
	    {
		mexErrMsgTxt("Invalid region of interest.");
	    }
	}
    }

    /* Create the start and stop indices. */
    if (nrhs == 2)
    {
	for (i = 0; i < sn; i++)
	{
	    start_indices[i] = 0;
	    stop_indices[i] = sdim[i] - 1;
	}
    }
    else
    {
	region_of_interest = mxGetPr(prhs[2]);
	for (i = 0; i < sn; i++)
	{
	    start_indices[i] = region_of_interest[i] - 1;
	    stop_indices[i] = region_of_interest[i+sn] - 1;
	}
    }	
    for (i = sn; i < DIMENSIONALITY; i++)
    {
	start_indices[i] = 0;
	stop_indices[i] = 0;
    }

    /* Compute the output dimensions. */
    for (i = 0; i < DIMENSIONALITY; i++)
    {
	result_dimensions[i] = stop_indices[i] - start_indices[i] + 1;
    }

    /* Create the output array. */
    resultarray = mxCreateNumericArray(DIMENSIONALITY, result_dimensions,
				       mxDOUBLE_CLASS,
				       signal_is_complex ? mxCOMPLEX : mxREAL);

    /* Call the computational routine. */
    conv3(mxGetPr(prhs[0]),
	  mxGetPr(prhs[1]),
	  mxGetPr(resultarray),
	  signal_dimensions,
	  kernel_dimensions,
	  result_dimensions,
	  start_indices,
	  stop_indices);

    if (signal_is_complex)
	conv3(mxGetPi(prhs[0]),
	      mxGetPr(prhs[1]),
	      mxGetPi(resultarray),
	      signal_dimensions,
	      kernel_dimensions,
	      result_dimensions,
	      start_indices,
	      stop_indices);

    /* Output the computed result. */
    plhs[0] = resultarray;
}
