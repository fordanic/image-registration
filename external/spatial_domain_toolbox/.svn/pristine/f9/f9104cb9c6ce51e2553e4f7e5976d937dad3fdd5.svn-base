#include "mex.h"

/*
 * r = POLYEXP_SOLVE_SYSTEM(BASIS, CONVRES_F, CONVRES_C, IS_REAL)
 * 
 * Helper function for polyexp.m, solving equations of the form (3.9)
 * in Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation
 * and Motion Estimation", when the basis functions are monomials. See
 * polyexp.m for the meaning of the parameters.
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
    int N, M;
    const int *basisdims;
    double *basis;
    const mxArray *convres_f_array;
    const mxArray *convres_c_array;
    double **G;
    double **h_r;
    double **h_i;
    int indices[3];
    int index;
    int dimensionality;
    int num_outdims;
    int outdims[4];
    int *cout_dimensions = NULL;
    mxArray *r_array;
    double *r_r;
    double *r_i;
    mxArray *cout_array = NULL;
    double *cout = NULL;
    int num_elements;
    mxArray *Qmatrix;
    double *Q;
    mxArray *qvector;
    double *q_r;
    double *q_i;
    double *p_r;
    double *p_i;
    mxArray *input[2];
    mxArray *output[1];
    int num_in, num_out;
    int is_real;
    int cout_needed = 0;
    char *cout_func_name = NULL;
    int cout_elements = 0;
    mxArray *cout_arguments[5];
    mxArray *this_cout_array;
    double *this_cout;
    int cout_num_arguments;
    
    
    /* Check the number of input and output arguments. */
    if (nrhs < 4)
	mexErrMsgTxt("Too few input arguments.");
    if (nrhs > 7)
	mexErrMsgTxt("Too many input arguments.");
    if ((nlhs > 1 && nrhs < 6) || nlhs > 2)
	mexErrMsgTxt("Too many output arguments.");
    if (nlhs < 2 && nrhs > 4)
	mexErrMsgTxt("Too few output arguments.");

    /* Check the formats of the input arguments. */
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0])
	|| mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
	mexErrMsgTxt("Unexpected format for basis.");
    }

    if (mxGetNumberOfDimensions(prhs[0]) != 2)
	mexErrMsgTxt("basis must be a matrix.");

    if (!mxIsCell(prhs[1]))
	mexErrMsgTxt("convres_f is expected to be a cell array.");
    
    if (!mxIsCell(prhs[2]))
	mexErrMsgTxt("convres_c is expected to be a cell array.");

    /* We can't check mxIsDouble(prhs[3]) since it is a logical
     * variable, which in matlab 6.5 has its own class. Neither can we
     * check mxIsNumeric(prhs[3]) since logicals are no longer
     * considered numeric in matlab 6.5.
     */
    if (mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]))
    {
	mexPrintf("%d %d %d\n",mxIsNumeric(prhs[3]), mxIsComplex(prhs[3])
		  ,mxIsSparse(prhs[3]));
	mexErrMsgTxt("Unexpected format for is_real.");
    }

    if (nrhs > 4 && !mxIsChar(prhs[4]))
	mexErrMsgTxt("cout_func is expected to be a string.");

    if (nrhs > 5 && (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5])
		     || mxIsSparse(prhs[5]) || !mxIsDouble(prhs[5])))
    {
	mexErrMsgTxt("Unexpected format for G0.");
    }

    /* We have no expectations on cout_data (prhs[6]). */

    /* Regardless whether prhs[3] is a double with logical flag
     * (matlab 5.3) or a variable of class logical (matlab 6.5),
     * mxGetScalar() automatically converts it to a double.
     */
    is_real = (int) mxGetScalar(prhs[3]);

    if (nlhs > 1)
    {
	if (nrhs <= 4)
	    mexErrMsgTxt("Too many output arguments. Cout requires a cout_func  parameter.");

	cout_needed = 1;
    }

    basisdims = mxGetDimensions(prhs[0]);
    N = basisdims[0]; /* Number of signal dimensions. */
    M = basisdims[1]; /* Number of basis functions. */
    basis = mxGetPr(prhs[0]);
    
    convres_f_array = prhs[1];
    convres_c_array = prhs[2];

    /* We want to set up a matrix and a vector with the pointers to
     * the start of the elements in the equation system.
     */
    G = mxCalloc(M * M, sizeof(*G));
    h_r = mxCalloc(M, sizeof(*h_r));
    if (!is_real)
	h_i = mxCalloc(M, sizeof(*h_i));
    
    for (i = 0; i < M; i++)
    {
	for (k = 0; k < N; k++)
	    indices[k] = (int) basis[N * i + k];

	index = mxCalcSingleSubscript(convres_f_array, N, indices);
	h_r[i] = mxGetPr(mxGetCell(convres_f_array, index));
	if (!is_real)
	    h_i[i] = mxGetPi(mxGetCell(convres_f_array, index));
	
	for (j = 0; j < M; j++)
	{
	    for (k = 0; k < N; k++)
		indices[k] = (int) (basis[N * i + k] + basis[N * j + k]);

	    index = mxCalcSingleSubscript(convres_c_array, N, indices);
	    mxGetCell(convres_c_array, index);
	    mxGetPr(mxGetCell(convres_c_array, index));
	    G[i + j * M] = mxGetPr(mxGetCell(convres_c_array, index));
	}
    }

    dimensionality = mxGetNumberOfDimensions(mxGetCell(convres_c_array,
						       index));
    num_outdims = dimensionality;
    num_elements = 1;
    for (k = 0; k < num_outdims; k++)
    {
	outdims[k] = mxGetDimensions(mxGetCell(convres_c_array, index))[k];
	num_elements *= outdims[k];
    }
    if (M > 1)
    {
	outdims[num_outdims] = M;
	num_outdims++;
    }
    
    r_array = mxCreateNumericArray(num_outdims, outdims, mxDOUBLE_CLASS,
				   is_real ? mxREAL : mxCOMPLEX);
    r_r = mxGetPr(r_array);
    if (!is_real)
	r_i = mxGetPi(r_array);
    
    Qmatrix = mxCreateDoubleMatrix(M, M, mxREAL);
    Q = mxGetPr(Qmatrix);
    input[0] = Qmatrix;
    
    qvector = mxCreateDoubleMatrix(M, 1, is_real ? mxREAL : mxCOMPLEX);
    q_r = mxGetPr(qvector);
    if (!is_real)
	q_i = mxGetPi(qvector);
    input[1] = qvector;
    
    num_in = 2;
    num_out = 1;

    if (cout_needed)
    {
	cout_num_arguments = 4;
	cout_arguments[0] = Qmatrix;
	cout_arguments[1] = (mxArray *) prhs[5];
	cout_arguments[2] = qvector;

	cout_func_name = mxCalloc(mxGetN(prhs[4]) + 1, 1);
	mxGetString(prhs[4], cout_func_name, mxGetN(prhs[4]) + 1);

	if (nrhs == 7)
	{
	    cout_arguments[4] = (mxArray *) prhs[6]; /* cout_data */
	    cout_num_arguments = 5;
	}
    }
    
    for (k = 0; k < num_elements; k++)
    {
	for (i = 0; i < M; i++)
	{
	    for (j = 0; j < M; j++)
		Q[i + j * M] = G[i + j * M][k];
	    
	    q_r[i] = h_r[i][k];
	    if (!is_real)
		q_i[i] = h_i[i][k];
	}
	mexCallMATLAB(num_out, output, num_in, input, "\\");
	if (is_real)
	{
	    p_r = mxGetPr(output[0]);
	    for (i = 0; i < M; i++)
		r_r[k + i * num_elements] = p_r[i];
	}
	else
	{
	    p_r = mxGetPr(output[0]);
	    p_i = mxGetPi(output[0]);
	    /* Although the signal has imaginary parts somewhere, the
	     * pointwise solutions may be real anywhere.
	     */
	    if (p_i) {
		for (i = 0; i < M; i++)
		{
		    r_r[k + i * num_elements] = p_r[i];
		    r_i[k + i * num_elements] = p_i[i];
		}
	    }
	    else {
		for (i = 0; i < M; i++) {
		    r_r[k + i * num_elements] = p_r[i];
		    r_i[k + i * num_elements] = 0.0;
		}
	    }
	}

	if (cout_needed)
	{
	    cout_arguments[3] = output[0];
	    mexCallMATLAB(1, &this_cout_array,
			  cout_num_arguments, cout_arguments, cout_func_name);
	    this_cout = mxGetPr(this_cout_array);

	    /* Create the output array if it doesn't already exist. */
	    if (cout == NULL)
	    {
		int cout_num_dims;
		const int *cout_dims;
		cout_num_dims = mxGetNumberOfDimensions(this_cout_array);
		cout_dims = mxGetDimensions(this_cout_array);
		cout_dimensions = (int *)mxCalloc(dimensionality
						  + cout_num_dims,
						  sizeof(int));
		for (i = 0; i < dimensionality; i++)
		    cout_dimensions[i] = outdims[i];
		for (; i < dimensionality + cout_num_dims; i++)
		    cout_dimensions[i] = cout_dims[i - dimensionality];
		
		cout_elements = (mxGetM(this_cout_array)
				 * mxGetN(this_cout_array));
		
		cout_array = mxCreateNumericArray(dimensionality
						  + cout_num_dims,
						  cout_dimensions,
						  mxDOUBLE_CLASS, mxREAL);
		cout = mxGetPr(cout_array);
	    }
	    
	    for (i = 0; i < cout_elements; i++)
		cout[k + i * num_elements] = this_cout[i];
	    
	    mxDestroyArray(this_cout_array);
	}
	
	mxDestroyArray(output[0]);
    }
    
    /* Output the computed result. */
    plhs[0] = r_array;
    if (cout_needed)
	plhs[1] = cout_array;
}
