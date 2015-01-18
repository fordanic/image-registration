#include "mex.h"
#include <string.h>

#define DEBUG 0

/*
 * NEIGHBORHOODLOOP
 * 
 * See the help text in neighborhoodloop.m for a description.
 * 
 * Author: Gunnar Farnebäck
 *         Medical Informatics
 *         Linköping University, Sweden
 *         gunnar@imt.liu.se
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int N;
    int *neighborhoodsize;
    int k, r;
    int x;
    int number_of_arrays;
    int number_of_points;
    int *number_of_elements;
    int *number_of_elements2;
    int max_number_of_dimensions;

    /* Set if the callback function is a function handle
     * or an inline object.
     */
    int function_is_handle;

    int function_name_length;
    char *function_name;
    int first_extra_parameter;
    int number_of_extra_parameters;
    mxArray **input_arrays;
    mxArray **output_arrays;
    int *dims;
    int num_dims;

    int *neighborhood_coords;
    int *neighborhood_delta;
    int number_of_neighborhood_points;
    int *xcoords;

    /* Check the input and output arguments. */
    
    /* First we expect a scalar. */
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0])
	|| mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0])
	|| mxGetNumberOfDimensions(prhs[0]) > 2
	|| mxGetDimensions(prhs[0])[0] != 1
	|| mxGetDimensions(prhs[0])[1] != 1)
    {
	mexErrMsgTxt("N is expected to be a scalar.");
    }

    N = (int) mxGetScalar(prhs[0]);
    if ((double) N != mxGetScalar(prhs[0]))
	mexErrMsgTxt("N is expected to be an integer.");

    if (N < 0)
	mexErrMsgTxt("N must not be negative.");

    /* Second we expect another scalar or a vector of length N. */
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1])
	|| mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1])
	|| mxGetNumberOfDimensions(prhs[1]) > 2
	|| ((mxGetDimensions(prhs[1])[0] != 1
	     || mxGetDimensions(prhs[1])[1] != 1)
	    && (mxGetDimensions(prhs[1])[0] 
		* mxGetDimensions(prhs[1])[1] != N)))
    {
	mexErrMsgTxt("size is expected to be a scalar or a vector of length N.");
    }

    neighborhoodsize = mxCalloc(N, sizeof(*neighborhoodsize));
    if (neighborhoodsize == NULL)
	mexErrMsgTxt("Failed to allocate an array.");
    number_of_neighborhood_points = 1;
    for (r = 0; r < N; r++)
    {
	double size;
	if (mxGetDimensions(prhs[1])[0] * mxGetDimensions(prhs[1])[1] == 1)
	    size = mxGetScalar(prhs[1]);
	else
	    size = mxGetPr(prhs[1])[r];
	neighborhoodsize[r] = (int) size;
	number_of_neighborhood_points *= neighborhoodsize[r];
	if ((double) neighborhoodsize[r] != size)
	    mexErrMsgTxt("size is expected to be integer.");

	if (neighborhoodsize[r] < 0 || neighborhoodsize[r] % 2 != 1)
	mexErrMsgTxt("size must be positive and odd.");
    }

    /* Then a number of arrays. */
    number_of_arrays = 0;
    number_of_points = 1;
    max_number_of_dimensions = 0;

    for (k = 2; k < nrhs; k++)
    {
	int n;

	if (mxIsChar(prhs[k]))
	{
	    function_is_handle = 0;
	    break;
	}

	if (mxGetNumberOfElements(prhs[k]) == 1)
	{
	    function_is_handle = 1;
	    break;
	}
	
	if (!mxIsNumeric(prhs[k]) || mxIsSparse(prhs[k])
	    || !mxIsDouble(prhs[k]))
	{
	    mexPrintf("Expected array for argument %d.", k + 1);
	    mexErrMsgTxt("");
	}

	n = mxGetNumberOfDimensions(prhs[k]);
	if (n < N)
	{
	    mexPrintf("Argument %d has fewer dimensions than N.", k + 1);
	    mexErrMsgTxt("");
	}
	
	if (max_number_of_dimensions < n)
	    max_number_of_dimensions = n;
	
	for (r = 0; r < N; r++)
	{
	    int dim1, dim2;
	    
	    if (mxGetNumberOfDimensions(prhs[2]) <= r)
		dim1 = 1;
	    else
		dim1 = mxGetDimensions(prhs[2])[r];

	    if (k == 2)
		number_of_points *= dim1;
	    else
	    {
		if (mxGetNumberOfDimensions(prhs[k]) <= r)
		    dim2 = 1;
		else
		    dim2 = mxGetDimensions(prhs[k])[r];
		
		if (dim1 != dim2)
		{
		    mexPrintf("Argument %d and argument 3 have incompatible sizes.",
			      k + 2);
		    mexErrMsgTxt("");
		}
	    }
	}
	number_of_arrays++;
    }

    if (number_of_arrays == 0)
	mexErrMsgTxt("Expected array for argument 3.");

    if (k == nrhs)
	mexErrMsgTxt("No function name provided.");

    if (!function_is_handle)
    {
	function_name_length = mxGetM(prhs[k]) * mxGetN(prhs[k]);
	function_name = mxCalloc(function_name_length + 1, 1);
	if (function_name == NULL)
	    mexErrMsgTxt("Failed to allocate space for function name.");
	if (mxGetString(prhs[k], function_name, function_name_length + 1) != 0)
	    mexErrMsgTxt("Failed to convert function name to string.");
    }
    else
	function_name = "feval";

    first_extra_parameter = k + 1;
    number_of_extra_parameters = nrhs - first_extra_parameter;

    if (nlhs == 0)
	nlhs = 1;
    
    /* Create some auxiliary arrays. */
    number_of_elements = mxCalloc(number_of_arrays,
				  sizeof(*number_of_elements));
    if (number_of_elements == NULL)
	mexErrMsgTxt("Failed to allocate an array.");

    number_of_elements2 = mxCalloc(nlhs, sizeof(*number_of_elements2));
    if (number_of_elements2 == NULL)
	mexErrMsgTxt("Failed to allocate an array.");
    
    
    /* Create array of input arrays for the callback. */
    input_arrays = mxCalloc(number_of_arrays + function_is_handle
			    + number_of_extra_parameters,
			    sizeof(*input_arrays));
    if (input_arrays == NULL)
	mexErrMsgTxt("Failed to allocate an array.");

    /* If the function is a handle, we need to pass the handle as the
     * first argument to the feval call.
     */
    if (function_is_handle)
	input_arrays[0] = (mxArray *) prhs[k];
    
    /* Populate it with neighborhood arrays. */
    dims = mxCalloc(max_number_of_dimensions, sizeof(*dims));
    if (dims == NULL)
	mexErrMsgTxt("Failed to allocate an array.");

    for (k = 0; k < number_of_arrays; k++)
    {
	num_dims = 0;
	number_of_elements[k] = 1;
	for (r = 0; r < mxGetNumberOfDimensions(prhs[k + 2]); r++)
	{
	    int d = mxGetDimensions(prhs[k + 2])[r];

	    if (r >= N) {
		dims[r] = d;
		number_of_elements[k] *= d;
	    }
	    else
		dims[r] = neighborhoodsize[r];
	    
	    num_dims++;
	}
	for (; num_dims < 2; num_dims++)
	    dims[num_dims] = 1;
	
	input_arrays[k + function_is_handle] =
	    mxCreateNumericArray(num_dims, dims, mxDOUBLE_CLASS,
				 mxIsComplex(prhs[k + 2]) ?
				 mxCOMPLEX : mxREAL);
	if (input_arrays[k + function_is_handle] == NULL)
	    mexErrMsgTxt("Failed to create an array.");
    }
    
    /* Continue population with the fixed parameters. */
    for (k = 0; k < number_of_extra_parameters; k++)
	input_arrays[number_of_arrays + function_is_handle + k] =
	    (mxArray *) prhs[first_extra_parameter + k];

    /* Create array of output arrays for the callback. */
    output_arrays = mxCalloc(nlhs, sizeof(*output_arrays));
    if (output_arrays == NULL)
	mexErrMsgTxt("Failed to allocate an array.");

    /* Set up the first N dimension sizes. */
    for (r = 0; r < N; r++)
    {
	if (r < mxGetNumberOfDimensions(prhs[2]))
	    dims[r] = mxGetDimensions(prhs[2])[r];
	else
	    dims[r] = 1;
    }

    /* Create array to store x coordinates in. */
    xcoords = mxCalloc(N, sizeof(*xcoords));
    if (xcoords == NULL)
	mexErrMsgTxt("Failed to allocate an array.");

    /* Set up arrays to keep track of the neighborhood points. */
    neighborhood_coords = mxCalloc(number_of_neighborhood_points * N,
				   sizeof(*neighborhood_coords));
    if (neighborhood_coords == NULL)
	mexErrMsgTxt("Failed to allocate an array.");
    
    neighborhood_delta = mxCalloc(number_of_neighborhood_points,
				   sizeof(*neighborhood_delta));
    if (neighborhood_delta == NULL)
	mexErrMsgTxt("Failed to allocate an array.");
    
    for (k = 0; k < number_of_neighborhood_points; k++)
    {
	int kk = k;
	int delta = 0;
	for (r = 0; r < N; r++)
	{
	    int coord = kk % neighborhoodsize[r] - neighborhoodsize[r] / 2;
	    kk /= neighborhoodsize[r];
	    neighborhood_coords[k * N + r] = coord;
	}
	for (r = N - 1; r >= 0; r--)
	{
	    int coord = neighborhood_coords[k * N + r];
	    delta *= dims[r];
	    delta += coord;
	}
	neighborhood_delta[k] = delta;
    }
    
    /* Time to start looping. */
    for (x = 0; x < number_of_points; x++)
    {
	int xx = x;
	for (k = 0; k < N; k++)
	{
	    xcoords[k] = xx % dims[k];
	    xx /= dims[k];
	}
	
	/* Copy the content of neighborhoods at this point of each array. */
	for (k = 0; k < number_of_arrays; k++)
	{
	    int components;
	    double *source;
	    double *target;
	    /* Loop over the real and complex components. */
	    for (components = 0; components < 2; components++)
	    {
		if (components == 0)
		{
		    source = mxGetPr(prhs[k + 2]);
		    target = mxGetPr(input_arrays[k + function_is_handle]);
		}
		else
		{
		    if (!mxIsComplex(input_arrays[k + function_is_handle]))
			break;
		    source = mxGetPi(prhs[k + 2]);
		    target = mxGetPi(input_arrays[k + function_is_handle]);
		}
		/* Copy the neighborhood for this component. */
		for (r = 0; r < number_of_elements[k]; r++)
		{
		    int s;
		    for (s = 0; s < number_of_neighborhood_points; s++)
		    {
			double value = 1.0;
			int t;
			for (t = 0; t < N; t++)
			{
			    int coord = xcoords[t] + neighborhood_coords[t + s * N];
			    if (coord < 0 || coord >= dims[t])
			    {
				value = 0.0;
				break;
			    }
			}
			
			if (value != 0.0)
			    value = source[r * number_of_points + x
					   + neighborhood_delta[s]];
			
			target[r * number_of_neighborhood_points + s] = value;
		    }
		}
	    }
	}

	/* Call back to matlab. */
	mexCallMATLAB(nlhs, output_arrays,
		      number_of_arrays + number_of_extra_parameters
		      + function_is_handle,
		      input_arrays, function_name);

	/* If this is the first turn of the loop, allocate the output
         * arrays from this function. We had to wait with this until
         * we could know their sizes.
	 */
	if (x == 0) {
	    /* Check the maximum number of dimensions. */
	    max_number_of_dimensions = 0;
	    for (k = 0; k < nlhs; k++)
	    {
		int M = N + mxGetNumberOfDimensions(output_arrays[k]);
		if (max_number_of_dimensions < M)
		    max_number_of_dimensions = M;
	    }

	    /* Reallocate dims. */
	    mxFree(dims);
	    dims = mxCalloc(max_number_of_dimensions, sizeof(*dims));
	    if (dims == NULL)
		mexErrMsgTxt("Failed to allocate an array.");
	    
	    /* Set up the first N dimension sizes. */
	    for (r = 0; r < N; r++)
	    {
		if (r < mxGetNumberOfDimensions(prhs[2]))
		    dims[r] = mxGetDimensions(prhs[2])[r];
		else
		    dims[r] = 1;
	    }
	    
	    /* Loop over each array. */
	    for (k = 0; k < nlhs; k++)
	    {
		/* Add the remaining dimensions. */
		num_dims = N;
		number_of_elements2[k] = 1;
		for (r = 0;
		     r < mxGetNumberOfDimensions(output_arrays[k]);
		     r++)
		{
		    int d = mxGetDimensions(output_arrays[k])[r];
		    dims[N + r] = d;
		    num_dims++;
		    number_of_elements2[k] *= d;
		}
		plhs[k] = mxCreateNumericArray(num_dims, dims, mxDOUBLE_CLASS,
					       mxIsComplex(output_arrays[k]) ?
					       mxCOMPLEX : mxREAL);
		if (plhs[k] == NULL)
		    mexErrMsgTxt("Failed to create an array.");
	    }
	}

	/* Copy output. */
	for (k = 0; k < nlhs; k++)
	{
	    double *source = mxGetPr(output_arrays[k]);
	    double *target = mxGetPr(plhs[k]);
	    for (r = 0; r < number_of_elements2[k]; r++)
		target[r * number_of_points + x] = source[r];

	    if (mxIsComplex(output_arrays[k]))
	    {
		if (!mxIsComplex(plhs[k]))
		{
		    /* Oh, looks like we need to be complex after all.
		     * Allocate a new array and copy the real data.
		     */
		    mxArray *new_array =
			mxCreateNumericArray(mxGetNumberOfDimensions(plhs[k]),
					     mxGetDimensions(plhs[k]),
					     mxDOUBLE_CLASS, mxCOMPLEX);
		    if (new_array == NULL)
			mexErrMsgTxt("Failed to create an array.");
		    
		    for (r = 0;
			 r < number_of_points * number_of_elements2[k];
			 r++)
		    {
			mxGetPr(new_array)[r] = mxGetPr(plhs[k])[r];
		    }

		    mxDestroyArray(plhs[k]);
		    plhs[k] = new_array;
		}
		
		source = mxGetPi(output_arrays[k]);
		target = mxGetPi(plhs[k]);
		for (r = 0; r < number_of_elements2[k]; r++)
		    target[r * number_of_points + x] = source[r];
	    }

	    /* Free array space. */
	    mxDestroyArray(output_arrays[k]);
	}
    }

    mxFree(neighborhoodsize);
    if (!function_is_handle)
	mxFree(function_name);
    mxFree(number_of_elements);
    mxFree(number_of_elements2);
    mxFree(input_arrays);
    mxFree(output_arrays);
    mxFree(dims);
    mxFree(xcoords);
    mxFree(neighborhood_coords);
    mxFree(neighborhood_delta);
}
