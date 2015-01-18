#include "mex.h"
#include <string.h>

#define DEBUG 0

/*
 * ARRAYLOOP
 * 
 * See the help text in arrayloop.m for a description.
 * 
 * Author: Gunnar Farnebäck
 *         Medical Informatics
 *         Linköping University, Sweden
 *         gunnar@imt.liu.se
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int N;
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

    /* Then a number of arrays. */
    number_of_arrays = 0;
    number_of_points = 1;
    max_number_of_dimensions = 0;

    for (k = 1; k < nrhs; k++)
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
	    
	    if (mxGetNumberOfDimensions(prhs[1]) <= r)
		dim1 = 1;
	    else
		dim1 = mxGetDimensions(prhs[1])[r];

	    if (k == 1)
		number_of_points *= dim1;
	    else
	    {
		if (mxGetNumberOfDimensions(prhs[k]) <= r)
		    dim2 = 1;
		else
		    dim2 = mxGetDimensions(prhs[k])[r];
		
		if (dim1 != dim2)
		{
		    mexPrintf("Argument %d and argument 2 have incompatible sizes.",
			      k + 1);
		    mexErrMsgTxt("");
		}
	    }
	}
	number_of_arrays++;
    }

    if (number_of_arrays == 0)
	mexErrMsgTxt("Expected array for argument 2.");

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
    
    /* Populate it with squeezed arrays. */
    max_number_of_dimensions -= N;
    if (max_number_of_dimensions < 2)
	max_number_of_dimensions = 2;
    
    dims = mxCalloc(max_number_of_dimensions, sizeof(*dims));
    if (dims == NULL)
	mexErrMsgTxt("Failed to allocate an array.");

    for (k = 0; k < number_of_arrays; k++)
    {
	num_dims = 0;
	number_of_elements[k] = 1;
	for (r = 0; r < mxGetNumberOfDimensions(prhs[k + 1]) - N; r++)
	{
	    int d = mxGetDimensions(prhs[k + 1])[N + r];
	    dims[r] = d;
	    num_dims++;
	    number_of_elements[k] *= d;
	}
	for (; num_dims < 2; num_dims++)
	    dims[num_dims] = 1;
	
	input_arrays[k + function_is_handle] =
	    mxCreateNumericArray(num_dims, dims, mxDOUBLE_CLASS,
				 mxIsComplex(prhs[k + 1]) ?
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

    /* Time to start looping. */
    for (x = 0; x < number_of_points; x++)
    {
	/* Copy the content of each array at this point. */
	for (k = 0; k < number_of_arrays; k++)
	{
	    double *source = mxGetPr(prhs[k + 1]);
	    double *target = mxGetPr(input_arrays[k + function_is_handle]);
	    for (r = 0; r < number_of_elements[k]; r++)
		target[r] = source[r * number_of_points + x];

	    if (mxIsComplex(input_arrays[k + function_is_handle]))
	    {
		source = mxGetPi(prhs[k + 1]);
		target = mxGetPi(input_arrays[k + function_is_handle]);
		for (r = 0; r < number_of_elements[k]; r++)
		    target[r] = source[r * number_of_points + x];
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
		if (r < mxGetNumberOfDimensions(prhs[1]))
		    dims[r] = mxGetDimensions(prhs[1])[r];
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
}
