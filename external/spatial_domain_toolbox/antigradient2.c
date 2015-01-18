/*
 * ANTIGRADIENT
 * 
 * See the help text in antigradient.m for a description.
 * 
 * Author: Gunnar Farnebäck
 *         Medical Informatics
 *         Linköping University, Sweden
 *         gunnar@imt.liu.se
 */

#include "mex.h"
#include <string.h>

#define RECURSION_SIZE_LIMIT 4

/* This can without cost be overdimensioned. */
#define MAX_LEVELS 30

/* Helper macro. */
#define VAL(x, y) ((x) ? (y) : 0)

struct data
{
  mxArray *A_array;
  double *lhs[MAX_LEVELS];
};

struct data data;

static void
clear_global_data()
{
  int k;
  data.A_array = NULL;
  for (k = 0; k < MAX_LEVELS; k++)
    data.lhs[k] = NULL;
}

/*************** 2D ****************/

static void
logmatrix2D(double *M, int m, int n, char *name, char *logfunction)
{
  mxArray *M_array;
  int dims[2];
  int i;
  mxArray *input_arrays[2];
  dims[0] = m;
  dims[1] = n;
  M_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  for (i = 0; i < m * n; i++)
    mxGetPr(M_array)[i] = M[i];
  
  input_arrays[0] = M_array;
  input_arrays[1] = mxCreateString(name);
  mexCallMATLAB(0, NULL, 2, input_arrays, logfunction);
}


static void
solve_directly2D(double *lhs, double *rhs, double *f_out, int M, int N)
{
  int s = M * N;
  int dims[2];
  mxArray *b_array;
  double *b;
  int i, j;
  mxArray *x_array;
  mxArray *input_arrays[2];

  if (data.A_array == NULL)
  {
    double *A;
    
    dims[0] = s;
    dims[1] = s;
    data.A_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    A = mxGetPr(data.A_array);
  
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	int index = j * M + i;
	A[index + s * index] = lhs[9 * index] + (lhs[9 * index] == 0.0);
	if (i > 0)
	  A[index + s * (index - 1)] = lhs[9 * index + 1];
	if (i < M-1)
	  A[index + s * (index + 1)] = lhs[9 * index + 2];
	if (j > 0)
	  A[index + s * (index - M)] = lhs[9 * index + 3];
	if (j < N-1)
	  A[index + s * (index + M)] = lhs[9 * index + 4];
	if (i > 0 && j > 0)
	  A[index + s * (index - M - 1)] = lhs[9 * index + 5];
	if (i > 0 && j < N-1)
	  A[index + s * (index + M - 1)] = lhs[9 * index + 6];
	if (i < M-1 && j > 0)
	  A[index + s * (index - M + 1)] = lhs[9 * index + 7];
	if (i < M-1 && j < N-1)
	  A[index + s * (index + M + 1)] = lhs[9 * index + 8];
      }
    
    for (i = 0; i < s*s; i++)
      A[i] += 1.0 / (s*s);
  }
  
  dims[0] = s;
  dims[1] = 1;
  b_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  b = mxGetPr(b_array);

  memcpy(b, rhs, M * N * sizeof(*rhs));
    
  input_arrays[0] = data.A_array;
  input_arrays[1] = b_array;
  mexCallMATLAB(1, &x_array, 2, input_arrays, "\\");
  memcpy(f_out, mxGetPr(x_array), s * sizeof(*f_out));
  mxDestroyArray(x_array);
  mxDestroyArray(b_array);
}


/* Gauss-Seidel smoothing iteration. Red-black ordering. */
static void
gauss_seidel2D(double *f, double *A, double *d, int M, int N)
{
  int pass;
  int i, j;
  int index;
  
  for (pass = 0; pass <= 1; pass++)
  {
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	double new_f;

	if ((i + j) % 2 != pass)
	  continue;
	
	index = i + j * M;
	if (A[9 * index] == 0.0)
	  continue;
	
	new_f = d[index];
	if (i > 0)
	  new_f -= A[9 * index + 1] * f[index - 1];

	if (i < M-1)
	  new_f -= A[9 * index + 2] * f[index + 1];

	if (j > 0)
	  new_f -= A[9 * index + 3] * f[index - M];

	if (j < N-1)
	  new_f -= A[9 * index + 4] * f[index + M];
	
	if (i > 0 && j > 0)
	  new_f -= A[9 * index + 5] * f[index - 1 - M];

	if (i > 0 && j < N-1)
	  new_f -= A[9 * index + 6] * f[index - 1 + M];

	if (i < M-1 && j > 0)
	  new_f -= A[9 * index + 7] * f[index + 1 - M];

	if (i < M-1 && j < N-1)
	  new_f -= A[9 * index + 8] * f[index + 1 + M];

	f[index] = new_f / A[9 * index];
      }
  }
}


static void
compute_residual2D(double *r, double *A, double *d, double *f, int M, int N)
{
  int i, j;

  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
      {
      int index = j * M + i;
      double residual = 0.0;
      if (A[9 * index] != 0.0)
      {
	residual = d[index] - A[9 * index] * f[index];
	if (i > 0)
	  residual -= A[9 * index + 1] * f[index - 1];
	
	if (i < M-1)
	  residual -= A[9 * index + 2] * f[index + 1];
	
	if (j > 0)
	  residual -= A[9 * index + 3] * f[index - M];
	
	if (j < N-1)
	  residual -= A[9 * index + 4] * f[index + M];
	
	if (i > 0 && j > 0)
	  residual -= A[9 * index + 5] * f[index - 1 - M];
	
	if (i > 0 && j < N-1)
	  residual -= A[9 * index + 6] * f[index - 1 + M];
	
	if (i < M-1 && j > 0)
	  residual -= A[9 * index + 7] * f[index + 1 - M];
	
	if (i < M-1 && j < N-1)
	  residual -= A[9 * index + 8] * f[index + 1 + M];
      }
      
      r[index] = residual;
    }
}


static void
downsample2D(double *rhs, int M, int N,
	     double *rhs_coarse, int Mhalf, int Nhalf,
	     double *weight, double *coarse_weight)
{
  int i, j;
  int index1;
  int index2;
  double c, n, s, w, e, nw, ne, sw, se;
  double sum;
  
  if (M % 2 == 0 && N % 2 == 0)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);

	nw = weight[index2];
	ne = weight[index2 + M];
	sw = weight[index2 + 1];
	se = weight[index2 + M + 1];
	sum = nw + ne + sw + se;
	coarse_weight[index1] = sum;

	if (sum > 0)
	{
	  rhs_coarse[index1] = 4 / sum * (nw * rhs[index2]
					  + ne * rhs[index2 + M]
					  + sw * rhs[index2 + 1]
					  + se * rhs[index2 + M + 1]);
	}
      }
  }
  
  if (M % 2 == 1 && N % 2 == 0)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	double result;
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	
	nw = 0.5 * VAL(i > 0, weight[index2 - 1]);
	ne = 0.5 * VAL(i > 0, weight[index2 + M - 1]);
	w  = weight[index2];
	e  = weight[index2 + M];
	sw = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
	se = 0.5 * VAL(i < Mhalf - 1, weight[index2 + M + 1]);
	sum = nw + ne + w + e + sw + se;
	coarse_weight[index1] = sum;

	if (sum > 0)
	{
	  result = w * rhs[index2] + e * rhs[index2 + M];
	  if (i > 0)
	    result += nw * rhs[index2 - 1] + ne * rhs[index2 + M - 1];
	  if (i < Mhalf - 1)
	    result += sw * rhs[index2 + 1] + se * rhs[index2 + M + 1];
	  
	  rhs_coarse[index1] = 4 / sum * result;
	}
      }
  }
  
  if (M % 2 == 0 && N % 2 == 1)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	double result;
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	
	nw = 0.5 * VAL(j > 0, weight[index2 - M]);
	sw = 0.5 * VAL(j > 0, weight[index2 - M + 1]);
	n  = weight[index2];
	s  = weight[index2 + 1];
	ne = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
	se = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M + 1]);
	sum = nw + sw + n + s + ne + se;
	coarse_weight[index1] = sum;

	if (sum > 0)
	{
	  result = n * rhs[index2] + s * rhs[index2 + 1];
	  if (j > 0)
	    result += nw * rhs[index2 - M] + sw * rhs[index2 - M + 1];
	  if (j < Nhalf - 1)
	    result += ne * rhs[index2 + M] + se * rhs[index2 + M + 1];
	  
	  rhs_coarse[index1] = 4 / sum * result;
	}
      }
  }

  if (M % 2 == 1 && N % 2 == 1)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	double result;
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	
	c  = weight[index2];
	n  = 0.5 * VAL(i > 0, weight[index2 - 1]);
	s  = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
	w  = 0.5 * VAL(j > 0, weight[index2 - M]);
	e  = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
	nw = 0.25 * VAL(i > 0 && j > 0, weight[index2 - M - 1]);
	ne = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + M - 1]);
	sw = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 - M + 1]);
	se = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + M + 1]);
	sum = c + n + s + w + e + nw + ne + sw + se;
	coarse_weight[index1] = sum;

	if (sum > 0)
	{
	  result = c * rhs[index2];
	  if (n > 0)
	    result += n * rhs[index2 - 1];
	  if (s > 0)
	    result += s * rhs[index2 + 1];
	  if (w > 0)
	    result += w * rhs[index2 - M];
	  if (e > 0)
	    result += e * rhs[index2 + M];
	  if (nw > 0)
	    result += nw * rhs[index2 - M - 1];
	  if (ne > 0)
	    result += ne * rhs[index2 + M - 1];
	  if (sw > 0)
	    result += sw * rhs[index2 - M + 1];
	  if (se > 0)
	    result += se * rhs[index2 + M + 1];
	
	  rhs_coarse[index1] = 4 / sum * result;
	}
      }
  }
}


static void
galerkin2D(int level, int M, int N, int Mhalf, int Nhalf,
	   double *weight, double *coarse_weight)
{
  int i, j;
  double *lhs;
  double *lhs_coarse;

  if (data.lhs[level + 1] != NULL)
    return;

  data.lhs[level + 1] = mxCalloc(9 * Mhalf * Nhalf,
				 sizeof(*data.lhs[level + 1]));
  lhs = data.lhs[level];
  lhs_coarse = data.lhs[level + 1];
  
  for (j = 0; j < Nhalf; j++)
    for (i = 0; i < Mhalf; i++)
    {
      int index1 = (j * Mhalf + i);
      int index2 = (2 * j * M + 2 * i);
      double stencil1[3][3];
      double stencil2[5][5];
      double stencil3[3][3];
      double mask1[3][3];
      int u, v;

      for (u = 0; u < 5; u++)
	for (v = 0; v < 5; v++)
	{
	  stencil2[u][v] = 0;
	  if (u < 3 && v < 3)
	  {
	    stencil1[u][v] = 0;
	    stencil3[u][v] = 0;
	  }
	}

      mask1[1][1] = coarse_weight[index1];
      mask1[0][1] = VAL(i > 0, coarse_weight[index1 - 1]);
      mask1[2][1] = VAL(i < Mhalf - 1, coarse_weight[index1 + 1]);
      mask1[1][0] = VAL(j > 0, coarse_weight[index1 - Mhalf]);
      mask1[1][2] = VAL(j < Nhalf - 1, coarse_weight[index1 + Mhalf]);
      mask1[0][0] = VAL(i > 0 && j > 0, coarse_weight[index1 - Mhalf - 1]);
      mask1[0][2] = VAL(i > 0 && j < Nhalf - 1, coarse_weight[index1 + Mhalf - 1]);
      mask1[2][0] = VAL(i < Mhalf - 1 && j > 0, coarse_weight[index1 - Mhalf + 1]);
      mask1[2][2] = VAL(i < Mhalf - 1 && j < Nhalf - 1, coarse_weight[index1 + Mhalf + 1]);
      
      if (M % 2 == 0 && N % 2 == 0)
      {
	double nw = weight[index2];
	double ne = weight[index2 + M];
	double sw = weight[index2 + 1];
	double se = weight[index2 + M + 1];

	double mean = (nw + sw + ne + se) / 4;

	/* If mean is 0 here we can short-circuit. */
	if (mean == 0)
	  continue;
	
	stencil1[1][1] = nw / mean;
	stencil1[1][2] = ne / mean;
	stencil1[2][1] = sw / mean;
	stencil1[2][2] = se / mean;
      }
      
      if (M % 2 == 1 && N % 2 == 0)
      {
	double nw = 0.5 * VAL(i > 0, weight[index2 - 1]);
	double ne = 0.5 * VAL(i > 0, weight[index2 + M - 1]);
	double w  = weight[index2];
	double e  = weight[index2 + M];
	double sw = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
	double se = 0.5 * VAL(i < Mhalf - 1, weight[index2 + M + 1]);

	double mean = (w + e + nw + ne + sw + se) / 4;

	/* If mean is 0 here we can short-circuit. */
	if (mean == 0)
	  continue;
	
	stencil1[1][1] = w / mean;
	stencil1[1][2] = e / mean;
	stencil1[2][1] = sw / mean;
	stencil1[2][2] = se / mean;
	stencil1[0][1] = nw / mean;
	stencil1[0][2] = ne / mean;
      }

      if (M % 2 == 0 && N % 2 == 1)
      {
	double nw = 0.5 * VAL(j > 0, weight[index2 - M]);
	double sw = 0.5 * VAL(j > 0, weight[index2 - M + 1]);
	double n  = weight[index2];
	double s  = weight[index2 + 1];
	double ne = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
	double se = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M + 1]);
	
	double mean = (n + s + nw + ne + sw + se) / 4;

	/* If mean is 0 here we can short-circuit. */
	if (mean == 0)
	  continue;
	
	stencil1[1][1] = n / mean;
	stencil1[2][1] = s / mean;
	stencil1[1][2] = ne / mean;
	stencil1[2][2] = se / mean;
	stencil1[1][0] = nw / mean;
	stencil1[2][0] = sw / mean;
      }
      
      if (M % 2 == 1 && N % 2 == 1)
      {
	double c  = weight[index2];
	double n  = 0.5 * VAL(i > 0, weight[index2 - 1]);
	double s  = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
	double w  = 0.5 * VAL(j > 0, weight[index2 - M]);
	double e  = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
	double nw = 0.25 * VAL(i > 0 && j > 0, weight[index2 - M - 1]);
	double ne = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + M - 1]);
	double sw = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 - M + 1]);
	double se = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + M + 1]);

	double mean = (c + s + e + w + n + se + sw + ne + nw) / 4;

	/* If mean is 0 here we can short-circuit. */
	if (mean == 0)
	  continue;
	
	stencil1[0][0] = nw / mean;
	stencil1[0][1] = n / mean;
	stencil1[0][2] = ne / mean;
	stencil1[1][0] = w / mean;
	stencil1[1][1] = c / mean;
	stencil1[1][2] = e / mean;
	stencil1[2][0] = sw / mean;
	stencil1[2][1] = s / mean;
	stencil1[2][2] = se / mean;
      }

      for (u = 0; u < 3; u++)
	for (v = 0; v < 3; v++)
	{
	  if (stencil1[u][v] != 0)
	  {
	    int index = 9 * (index2 + (u-1) + M*(v-1));
	    if (lhs[index] != 0.0)
	    {
	      stencil2[u+1][v+1] += stencil1[u][v] * lhs[index];
	      stencil2[u  ][v+1] += stencil1[u][v] * lhs[index + 1];
	      stencil2[u+2][v+1] += stencil1[u][v] * lhs[index + 2];
	      stencil2[u+1][v  ] += stencil1[u][v] * lhs[index + 3];
	      stencil2[u+1][v+2] += stencil1[u][v] * lhs[index + 4];
	      stencil2[u  ][v  ] += stencil1[u][v] * lhs[index + 5];
	      stencil2[u  ][v+2] += stencil1[u][v] * lhs[index + 6];
	      stencil2[u+2][v  ] += stencil1[u][v] * lhs[index + 7];
	      stencil2[u+2][v+2] += stencil1[u][v] * lhs[index + 8];
	    }
	  }
	}

      if (M % 2 == 0 && N % 2 == 0)
      {
	for (u = 1; u < 5; u++)
	  for (v = 1; v < 5; v++)
	  {
	    double alpha1, alpha2;
	    double nw, ne, sw, se;
	    int uu, vv;
	    double sum;
	    
	    if (u % 2 == 0)
	      alpha1 = 0.75;
	    else
	      alpha1 = 0.25;

	    if (v % 2 == 0)
	      alpha2 = 0.75;
	    else
	      alpha2 = 0.25;

	    uu = (u-1)/2;
	    vv = (v-1)/2;
	    nw = (1 - alpha1) * (1 - alpha2) * mask1[uu  ][vv  ];
	    ne = (1 - alpha1) *      alpha2  * mask1[uu  ][vv+1];
	    sw =      alpha1  * (1 - alpha2) * mask1[uu+1][vv  ];
	    se =      alpha1  *      alpha2  * mask1[uu+1][vv+1];

	    sum = nw + ne + sw + se;
	    if (sum > 0)
	    {
	      stencil3[uu  ][vv  ] += nw * stencil2[u][v] / sum;
	      stencil3[uu  ][vv+1] += ne * stencil2[u][v] / sum;
	      stencil3[uu+1][vv  ] += sw * stencil2[u][v] / sum;
	      stencil3[uu+1][vv+1] += se * stencil2[u][v] / sum;
	    }
	  }
      }

      if (M % 2 == 1 && N % 2 == 0)
      {
	for (u = 0; u < 5; u++)
	  for (v = 1; v < 5; v++)
	  {
	    double alpha1, alpha2;
	    double nw, ne, sw, se;
	    int uu, vv;
	    double sum;
	    if (u % 2 == 0)
	      alpha1 = 0;
	    else
	      alpha1 = 0.5;
	    
	    if (v % 2 == 0)
	      alpha2 = 0.75;
	    else
	      alpha2 = 0.25;
	    
	    uu = u/2;
	    vv = (v-1)/2;
	    nw = (1 - alpha1) * (1 - alpha2) * mask1[uu  ][vv  ];
	    ne = (1 - alpha1) *      alpha2  * mask1[uu  ][vv+1];
	    sw = 0;
	    se = 0;
	    if (alpha1 > 0)
	    {
	      sw = alpha1  * (1 - alpha2) * mask1[uu+1][vv  ];
	      se = alpha1  *      alpha2  * mask1[uu+1][vv+1];
	    }

	    sum = nw + ne + sw + se;
	    if (sum > 0)
	    {
	      stencil3[uu  ][vv  ] += nw * stencil2[u][v] / sum;
	      stencil3[uu  ][vv+1] += ne * stencil2[u][v] / sum;
	      if (alpha1 > 0)
	      {
		stencil3[uu+1][vv  ] += sw * stencil2[u][v] / sum;
		stencil3[uu+1][vv+1] += se * stencil2[u][v] / sum;
	      }
	    }
	  }
      }

      if (M % 2 == 0 && N % 2 == 1)
      {
	for (u = 1; u < 5; u++)
	  for (v = 0; v < 5; v++)
	  {
	    double alpha1, alpha2;
	    double nw, ne, sw, se;
	    int uu, vv;
	    double sum;
	    if (u % 2 == 0)
	      alpha1 = 0.75;
	    else
	      alpha1 = 0.25;
	    
	    if (v % 2 == 0)
	      alpha2 = 0;
	    else
	      alpha2 = 0.5;
	    
	    uu = (u-1)/2;
	    vv = v/2;
	    nw = (1 - alpha1) * (1 - alpha2) * mask1[uu  ][vv  ];
	    sw =      alpha1  * (1 - alpha2) * mask1[uu+1][vv  ];
	    ne = 0;
	    se = 0;
	    if (alpha2 > 0)
	    {
	      ne = (1 - alpha1) * alpha2 * mask1[uu  ][vv+1];
	      se =      alpha1  * alpha2 * mask1[uu+1][vv+1];
	    }

	    sum = nw + ne + sw + se;
	    if (sum > 0)
	    {
	      stencil3[uu  ][vv  ] += nw * stencil2[u][v] / sum;
	      stencil3[uu+1][vv  ] += sw * stencil2[u][v] / sum;
	      if (alpha2 > 0)
	      {
		stencil3[uu  ][vv+1] += ne * stencil2[u][v] / sum;
		stencil3[uu+1][vv+1] += se * stencil2[u][v] / sum;
	      }
	    }
	  }
      }

      if (M % 2 == 1 && N % 2 == 1)
      {
	for (u = 0; u < 5; u++)
	  for (v = 0; v < 5; v++)
	  {
	    double alpha1, alpha2;
	    double nw, ne, sw, se;
	    int uu, vv;
	    double sum;
	    if (u % 2 == 0)
	      alpha1 = 0;
	    else
	      alpha1 = 0.5;
	    
	    if (v % 2 == 0)
	      alpha2 = 0;
	    else
	      alpha2 = 0.5;

	    uu = u/2;
	    vv = v/2;
	    nw = (1 - alpha1) * (1 - alpha2) * mask1[uu  ][vv  ];
	    sw = 0;
	    if (alpha1 > 0)
	      sw =      alpha1  * (1 - alpha2) * mask1[uu+1][vv  ];
	    ne = 0;
	    if (alpha2 > 0)
	      ne = (1 - alpha1) * alpha2 * mask1[uu  ][vv+1];
	    se = 0;
	    if (alpha1 > 0 && alpha2 > 0)
	      se =      alpha1  * alpha2 * mask1[uu+1][vv+1];

	    sum = nw + ne + sw + se;
	    if (sum > 0)
	    {
	      stencil3[uu  ][vv  ] += nw * stencil2[u][v] / sum;
	      if (alpha1 > 0)
		stencil3[uu+1][vv  ] += sw * stencil2[u][v] / sum;
	      if (alpha2 > 0)
		stencil3[uu  ][vv+1] += ne * stencil2[u][v] / sum;
	      if (alpha1 > 0 && alpha2 > 0)
		stencil3[uu+1][vv+1] += se * stencil2[u][v] / sum;
	    }
	  }
      }

      lhs_coarse[9 * index1]     = stencil3[1][1];
      lhs_coarse[9 * index1 + 1] = stencil3[0][1];
      lhs_coarse[9 * index1 + 2] = stencil3[2][1];
      lhs_coarse[9 * index1 + 3] = stencil3[1][0];
      lhs_coarse[9 * index1 + 4] = stencil3[1][2];
      lhs_coarse[9 * index1 + 5] = stencil3[0][0];
      lhs_coarse[9 * index1 + 6] = stencil3[0][2];
      lhs_coarse[9 * index1 + 7] = stencil3[2][0];
      lhs_coarse[9 * index1 + 8] = stencil3[2][2];

      for (u = 1; u < 9; u++)
	if (lhs_coarse[9 * index1 + u] < 0)
	{
	  lhs_coarse[9 * index1] += lhs_coarse[9 * index1 + u];
	  lhs_coarse[9 * index1 + u] = 0;
	}
    }
}


/* Upsample and apply correction. Bilinear interpolation. */
static void
upsample2D(int M, int N,
	   double *v, int Mhalf, int Nhalf,
	   double *f_out, double *weight, double *coarse_weight)
{
  int i, j;
  int index1, index2;
  
  if (M % 2 == 0 && N % 2 == 0)
  {
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	double alpha1, alpha2;
	double nw, ne, sw, se;
	double sum;
	
	index1 = j * M + i;
	index2 = ((j + 1) / 2 - 1) * Mhalf + (i + 1) / 2 - 1;

	if (weight[index1] == 0)
	  continue;

	if (i % 2 == 0)
	  alpha1 = 0.75;
	else
	  alpha1 = 0.25;
	
	if (j % 2 == 0)
	  alpha2 = 0.75;
	else
	  alpha2 = 0.25;
	
	nw = (1 - alpha1) * (1 - alpha2) * VAL(i > 0 && j > 0, coarse_weight[index2]);
	ne = (1 - alpha1) *      alpha2  * VAL(i > 0 && j < N - 1, coarse_weight[index2 + Mhalf]);
	sw =      alpha1  * (1 - alpha2) * VAL(i < M - 1 && j > 0, coarse_weight[index2 + 1]);
	se =      alpha1  *      alpha2  * VAL(i < M - 1 && j < N - 1, coarse_weight[index2 + Mhalf + 1]);
	
	sum = nw + ne + sw + se;

	if (sum > 0)
	{
	  double contribution = 0;
	  
	  if (nw > 0)
	    contribution += nw * v[index2];
	  if (ne > 0)
	    contribution += ne * v[index2 + Mhalf];
	  if (sw > 0)
	    contribution += sw * v[index2 + 1];
	  if (se > 0)
	    contribution += se * v[index2 + Mhalf + 1];
	  
	  f_out[index1] += contribution / sum;
	}
      }
  }
    
  if (M % 2 == 1 && N % 2 == 0)
  {
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	double alpha1, alpha2;
	double nw, ne, sw, se;
	double sum;
	
	index1 = j * M + i;
	index2 = ((j + 1) / 2 - 1) * Mhalf + i / 2;

	if (weight[index1] == 0)
	  continue;

	if (i % 2 == 0)
	  alpha1 = 0.0;
	else
	  alpha1 = 0.5;
	
	if (j % 2 == 0)
	  alpha2 = 0.75;
	else
	  alpha2 = 0.25;
	
	nw = (1 - alpha1) * (1 - alpha2) * VAL(j > 0, coarse_weight[index2]);
	ne = (1 - alpha1) *      alpha2  * VAL(j < N - 1, coarse_weight[index2 + Mhalf]);
	sw =      alpha1  * (1 - alpha2) * VAL(i < M - 1 && j > 0, coarse_weight[index2 + 1]);
	se =      alpha1  *      alpha2  * VAL(i < M - 1 && j < N - 1, coarse_weight[index2 + Mhalf + 1]);
	
	sum = nw + ne + sw + se;

	if (sum > 0)
	{
	  double contribution = 0;
	  
	  if (nw > 0)
	    contribution += nw * v[index2];
	  if (ne > 0)
	    contribution += ne * v[index2 + Mhalf];
	  if (sw > 0)
	    contribution += sw * v[index2 + 1];
	  if (se > 0)
	    contribution += se * v[index2 + Mhalf + 1];
	  
	  f_out[index1] += contribution / sum;
	}
      }
  }
  
  if (M % 2 == 0 && N % 2 == 1)
  {
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	double alpha1, alpha2;
	double nw, ne, sw, se;
	double sum;
	
	index1 = j * M + i;
	index2 = (j / 2) * Mhalf + (i + 1) / 2 - 1;

	if (weight[index1] == 0)
	  continue;

	if (i % 2 == 0)
	  alpha1 = 0.75;
	else
	  alpha1 = 0.25;
	
	if (j % 2 == 0)
	  alpha2 = 0.0;
	else
	  alpha2 = 0.5;
	
	nw = (1 - alpha1) * (1 - alpha2) * VAL(i > 0, coarse_weight[index2]);
	ne = (1 - alpha1) *      alpha2  * VAL(i > 0 && j < N - 1, coarse_weight[index2 + Mhalf]);
	sw =      alpha1  * (1 - alpha2) * VAL(i < M - 1, coarse_weight[index2 + 1]);
	se =      alpha1  *      alpha2  * VAL(i < M - 1 && j < N - 1, coarse_weight[index2 + Mhalf + 1]);
	
	sum = nw + ne + sw + se;

	if (sum > 0)
	{
	  double contribution = 0;
	  
	  if (nw > 0)
	    contribution += nw * v[index2];
	  if (ne > 0)
	    contribution += ne * v[index2 + Mhalf];
	  if (sw > 0)
	    contribution += sw * v[index2 + 1];
	  if (se > 0)
	    contribution += se * v[index2 + Mhalf + 1];
	  
	  f_out[index1] += contribution / sum;
	}
      }
  }
  
  if (M % 2 == 1 && N % 2 == 1)
  {
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	double alpha1, alpha2;
	double nw, ne, sw, se;
	double sum;
	
	index1 = j * M + i;
	index2 = (j / 2) * Mhalf + i / 2;

	if (weight[index1] == 0)
	  continue;

	if (i % 2 == 0)
	  alpha1 = 0.0;
	else
	  alpha1 = 0.5;
	
	if (j % 2 == 0)
	  alpha2 = 0.0;
	else
	  alpha2 = 0.5;
	
	nw = (1 - alpha1) * (1 - alpha2) * coarse_weight[index2];
	ne = (1 - alpha1) *      alpha2  * VAL(j < N - 1, coarse_weight[index2 + Mhalf]);
	sw =      alpha1  * (1 - alpha2) * VAL(i < M - 1, coarse_weight[index2 + 1]);
	se =      alpha1  *      alpha2  * VAL(i < M - 1 && j < N - 1, coarse_weight[index2 + Mhalf + 1]);
	
	sum = nw + ne + sw + se;

	if (sum > 0)
	{
	  double contribution = 0;
	  
	  if (nw > 0)
	    contribution += nw * v[index2];
	  if (ne > 0)
	    contribution += ne * v[index2 + Mhalf];
	  if (sw > 0)
	    contribution += sw * v[index2 + 1];
	  if (se > 0)
	    contribution += se * v[index2 + Mhalf + 1];
	  
	  f_out[index1] += contribution / sum;
	}
      }
  }
}


/* Recursive multigrid function.*/
static void
poisson_multigrid2D(double *f, int level, double *rhs, double *weight,
		    int n1, int n2, int nm,
		    double *f_out,
		    int M, int N, int *directly_solved)
{
  int k;
  double *r;
  double *r_downsampled;
  double *coarse_weight;
  double *lhs = data.lhs[level];
  double *v;
  int Mhalf;
  int Nhalf;

  /* Solve a sufficiently small problem directly. */
  if (M < RECURSION_SIZE_LIMIT || N < RECURSION_SIZE_LIMIT)
  {
    solve_directly2D(lhs, rhs, f_out, M, N);
    *directly_solved = 1;
    return;
  }
  *directly_solved = 0;
  
  /* Initialize solution. */
  memcpy(f_out, f, M * N * sizeof(*f_out));
  
  /* Pre-smoothing. */
  for (k = 0; k < n1; k++)
    gauss_seidel2D(f_out, lhs, rhs, M, N);
  
  /* Compute residual. */
  r = mxCalloc(M * N, sizeof(*r));
  compute_residual2D(r, lhs, rhs, f_out, M, N);

  /* Downsample residual. */
  Mhalf = (M + 1) / 2;
  Nhalf = (N + 1) / 2;
  r_downsampled = mxCalloc(Mhalf * Nhalf, sizeof(*r_downsampled));
  coarse_weight = mxCalloc(Mhalf * Nhalf, sizeof(*coarse_weight));
  downsample2D(r, M, N, r_downsampled, Mhalf, Nhalf, weight, coarse_weight);
  galerkin2D(level, M, N, Mhalf, Nhalf, weight, coarse_weight);
  
  /* Recurse to compute a correction. */
  v = mxCalloc(Mhalf * Nhalf, sizeof(*v));
  for (k = 0; k < nm; k++)
  {
    int directly_solved;
    poisson_multigrid2D(v, level + 1, r_downsampled, coarse_weight,
			n1, n2, nm, v, Mhalf, Nhalf, &directly_solved);
    if (directly_solved)
      break;
  }
  
  upsample2D(M, N, v, Mhalf, Nhalf, f_out, weight, coarse_weight);
  
  /* Post-smoothing. */
  for (k = 0; k < n2; k++)
    gauss_seidel2D(f_out, lhs, rhs, M, N);

  /* Set the mean value to zero.
   *
   * FIXME: This should not be needed (I believe) and might indicate
   * some bug elsewhere.
   */
  if (1)
  {
    double sum = 0.0;
    int num_samples_in_mask = 0;
    double mean;
    int i;
    
    for (i = 0; i < M * N; i++)
      if (weight[i] != 0)
      {
	sum += f_out[i];
	num_samples_in_mask++;
      }
    
    mean = sum / num_samples_in_mask;
    for (i = 0; i < M * N; i++)
      if (weight[i] != 0)
	f_out[i] -= mean;
  }
  
  mxFree(r);
  mxFree(r_downsampled);
  mxFree(coarse_weight);
  mxFree(v);
}


/* It is assumed that f_out is initialized to zero when called. */
static void
poisson_full_multigrid2D(int level, double *rhs, double *weight,
			 int number_of_iterations, int M, int N, double *f_out)
{
  double *rhs_downsampled;
  double *coarse_weight;
  double *f_coarse;
  int k;
  
  /* Unless already coarsest scale, first recurse to coarser scale. */
  if (M >= RECURSION_SIZE_LIMIT && N >= RECURSION_SIZE_LIMIT)
  {
    /* Downsample right hand side. */
    int Mhalf = (M + 1) / 2;
    int Nhalf = (N + 1) / 2;
    rhs_downsampled = mxCalloc(Mhalf * Nhalf, sizeof(*rhs_downsampled));
    coarse_weight = mxCalloc(Mhalf * Nhalf, sizeof(*coarse_weight));
    downsample2D(rhs, M, N, rhs_downsampled, Mhalf, Nhalf,
		 weight, coarse_weight);
    galerkin2D(level, M, N, Mhalf, Nhalf, weight, coarse_weight);
    f_coarse = mxCalloc(Mhalf * Nhalf, sizeof(*f_coarse));
    poisson_full_multigrid2D(level + 1, rhs_downsampled, coarse_weight,
			     number_of_iterations, Mhalf, Nhalf, f_coarse);
    
    /* Upsample the coarse result. */
    upsample2D(M, N, f_coarse, Mhalf, Nhalf, f_out,
	       weight, coarse_weight);

    mxFree(f_coarse);
    mxFree(coarse_weight);
    mxFree(rhs_downsampled);
  }
  
  /* Perform number_of_iterations standard multigrid cycles. */
  for (k = 0; k < number_of_iterations; k++)
  {
    int directly_solved;
    poisson_multigrid2D(f_out, level, rhs, weight, 2, 2, 2, f_out, M, N,
			&directly_solved);
    if (directly_solved)
      break;
  }
}


static void
antigradient2D(double *g, double *mask, double mu, int number_of_iterations,
	       int M, int N, double *f_out)
{
  double *rhs;
  double *lhs;
  double *weight;
  double sum;
  double mean;
  int i, j;
  int num_samples_in_mask;

  clear_global_data();

  /* Compute left and right hand sides of Poisson problem with Neumann
   * boundary conditions, discretized by finite differences.
   */
  rhs = mxCalloc(M * N, sizeof(*rhs));
  lhs = mxCalloc(9 * M * N, sizeof(*lhs));
  data.lhs[0] = lhs;
  weight = mxCalloc(M * N, sizeof(*weight));
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
    {
      int index1 = j * M + i;
      int index2 = index1 + M * N;
      double d = 0.0;
      int N_missing = 0;
      int S_missing = 0;
      int W_missing = 0;
      int E_missing = 0;

      if (mask && mask[index1] == 0)
	continue;

      weight[index1] = 1;
      
      if (i == 0 || (mask && mask[index1 - 1] == 0))
	N_missing = 1;

      if (i == M - 1 || (mask && mask[index1 + 1] == 0))
	S_missing = 1;
      
      if (j == 0 || (mask && mask[index1 - M] == 0))
	W_missing = 1;

      if (j == N - 1 || (mask && mask[index1 + M] == 0))
	E_missing = 1;
      
      if (N_missing && !S_missing)
	d = g[index1 + 1] + g[index1];
      else if (!N_missing && S_missing)
	d = - g[index1] - g[index1 - 1];
      else if (!N_missing && !S_missing)
	d = 0.5 * (g[index1 + 1] - g[index1 - 1]);
      
      if (W_missing && !E_missing)
	d += g[index2 + M] + g[index2];
      else if (!W_missing && E_missing)
	d += - g[index2] - g[index2 - M];
      else if (!W_missing && !E_missing)
	d += 0.5 * (g[index2 + M] - g[index2 - M]);
      
      rhs[index1] = d;

      lhs[9 * index1] = -2 * ((!N_missing || !S_missing)
			      + (!W_missing || !E_missing));
      if (!N_missing)
	lhs[9 * index1 + 1] = 1 + (S_missing);
      if (!S_missing)
	lhs[9 * index1 + 2] = 1 + (N_missing);
      if (!W_missing)
	lhs[9 * index1 + 3] = 1 + (E_missing);
      if (!E_missing)
	lhs[9 * index1 + 4] = 1 + (W_missing);
    }
  
  /* Solve the equation system with the full multigrid algorithm.
   * Use W cycles and 2 presmoothing and 2 postsmoothing
   * Gauss-Seidel iterations.
   */
  poisson_full_multigrid2D(0, rhs, weight, number_of_iterations, M, N, f_out);
  
  /* Fix the mean value. */
  sum = 0.0;
  num_samples_in_mask = 0;
  for (i = 0; i < M * N; i++)
    if (weight[i])
    {
      sum += f_out[i];
      num_samples_in_mask++;
    }
  
  mean = sum / num_samples_in_mask;
  for (i = 0; i < M * N; i++)
    if (weight[i])
    {
      f_out[i] -= mean;
      f_out[i] += mu;
    }

  mxFree(rhs);
  mxFree(weight);
  for (i = 0; i < MAX_LEVELS; i++)
    if (data.lhs[i] != NULL)
      mxFree(data.lhs[i]);
  mxDestroyArray(data.A_array);
}


/*************** 3D ****************/

static void
logmatrix3D(double *M, int m, int n, int p, char *name, char *logfunction)
{
  mxArray *M_array;
  int dims[3];
  int i;
  mxArray *input_arrays[2];
  dims[0] = m;
  dims[1] = n;
  dims[2] = p;
  M_array = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  for (i = 0; i < m * n * p; i++)
    mxGetPr(M_array)[i] = M[i];
  
  input_arrays[0] = M_array;
  input_arrays[1] = mxCreateString(name);
  mexCallMATLAB(0, NULL, 2, input_arrays, logfunction);
}


static void
solve_directly3D(double *lhs, double *rhs, double *f_out, int M, int N, int P)
{
  int s = M * N * P;
  int MN = M * N;
  int dims[2];
  mxArray *b_array;
  double *b;
  int i, j, p;
  mxArray *x_array;
  mxArray *input_arrays[2];

  if (data.A_array == NULL)
  {
    double *A;
    
    dims[0] = s;
    dims[1] = s;
    data.A_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    A = mxGetPr(data.A_array);
  
    for (p = 0; p < P; p++)
      for (j = 0; j < N; j++)
	for (i = 0; i < M; i++)
	{
	  int index = (p * N + j) * M + i;
	  int k;
	  if (lhs[27 * index + 13] == 0)
	    A[index + s * index] = 1;
	  else
	  {
	    for (k = 0; k < 27; k++)
	    {
	      int u = (k % 3) - 1;
	      int v = ((k / 3) % 3) - 1;
	      int w = ((k / 9) % 3) - 1;
	      if (lhs[27 * index + k] != 0)
		A[index + s * (index + u + v * M + w * MN)] = lhs[27 * index + k];
	    }
	  }
	}
    
    for (i = 0; i < s*s; i++)
      A[i] += 1.0 / (s*s);
  }
  
  dims[0] = s;
  dims[1] = 1;
  b_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  b = mxGetPr(b_array);
  memcpy(b, rhs, s * sizeof(*rhs));
    
  input_arrays[0] = data.A_array;
  input_arrays[1] = b_array;
  mexCallMATLAB(1, &x_array, 2, input_arrays, "\\");
  memcpy(f_out, mxGetPr(x_array), s * sizeof(*f_out));
  mxDestroyArray(x_array);
  mxDestroyArray(b_array);
}


/* Gauss-Seidel smoothing iteration. Red-black ordering. */
static void
gauss_seidel3D(double *f, double *A, double *d, int M, int N, int P)
{
  int pass;
  int i, j, p;
  int index;
  int MN = M * N;
  
  for (pass = 0; pass <= 1; pass++)
  {
    for (p = 0; p < P; p++)
      for (j = 0; j < N; j++)
	for (i = 0; i < M; i++)
	{
	  double new_f;
	  int k;
	  
	  if ((i + j + p) % 2 != pass)
	    continue;
	  
	  index = (p * N + j) * M + i;

	  if (A[27 * index + 13] == 0)
	    continue;
	  
	  new_f = d[index];
	  for (k = 0; k < 27; k++)
	  {
	    int u = (k % 3) - 1;
	    int v = ((k / 3) % 3) - 1;
	    int w = ((k / 9) % 3) - 1;

	    if (k != 13 && A[27 * index + k] != 0)
	      new_f -= A[27 * index + k] * f[index + u + v * M + w * MN];
	  }

	  f[index] = new_f / A[27 * index + 13];
	}
  }
}


static void
compute_residual3D(double *r, double *A, double *d, double *f,
                   int M, int N, int P)
{
  int i, j, p;
  int MN = M * N;

  for (p = 0; p < P; p++)
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	int index = (p * N + j) * M + i;
	double residual = 0;
	if (A[27 * index + 13] != 0)
	{
	  int k;
	  residual = d[index];
	  for (k = 0; k < 27; k++)
	  {
	    int u = (k % 3) - 1;
	    int v = ((k / 3) % 3) - 1;
	    int w = ((k / 9) % 3) - 1;

	    if (A[27 * index + k] != 0)
	      residual -= A[27 * index + k] * f[index + u + v * M + w * MN];
	  }
	}

	r[index] = residual;
      }
}

static void
downsample3D(double *rhs, int M, int N, int P,
             double *rhs_coarse, int Mhalf, int Nhalf, int Phalf,
             double *weight, double *coarse_weight)
{
  int i, j, p;
  int index1;
  int index2;
  int MN = M * N;
  double w[3][3][3];
  double sum;
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[1][1][1] = weight[index2];
          w[1][1][2] = weight[index2 + MN];
          w[1][2][1] = weight[index2 + M];
          w[1][2][2] = weight[index2 + MN + M];
          w[2][1][1] = weight[index2 + 1];
          w[2][1][2] = weight[index2 + MN + 1];
          w[2][2][1] = weight[index2 + M + 1];
          w[2][2][2] = weight[index2 + MN + M + 1];
          
          sum = w[1][1][1] + w[1][1][2] + w[1][2][1] + w[1][2][2] + w[2][1][1] + w[2][1][2] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          w[0][1][2] = 0.5 * VAL(i > 0, weight[index2 + MN - 1]);
          w[0][2][1] = 0.5 * VAL(i > 0, weight[index2 + M - 1]);
          w[0][2][2] = 0.5 * VAL(i > 0, weight[index2 + MN + M - 1]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = weight[index2 + MN];
          w[1][2][1] = weight[index2 + M];
          w[1][2][2] = weight[index2 + MN + M];
          w[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          w[2][1][2] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + MN + 1]);
          w[2][2][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + M + 1]);
          w[2][2][2] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[0][1][1] + w[0][1][2] + w[0][2][1] + w[0][2][2] + w[1][1][1] + w[1][1][2] + w[1][2][1] + w[1][2][2] + w[2][1][1] + w[2][1][2] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[0][1][1] > 0)
              result += w[0][1][1] * rhs[index2 - 1];
            if (w[0][1][2] > 0)
              result += w[0][1][2] * rhs[index2 + MN - 1];
            if (w[0][2][1] > 0)
              result += w[0][2][1] * rhs[index2 + M - 1];
            if (w[0][2][2] > 0)
              result += w[0][2][2] * rhs[index2 + MN + M - 1];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          w[1][0][2] = 0.5 * VAL(j > 0, weight[index2 + MN - M]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = weight[index2 + MN];
          w[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          w[1][2][2] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + MN + M]);
          w[2][0][1] = 0.5 * VAL(j > 0, weight[index2 - M + 1]);
          w[2][0][2] = 0.5 * VAL(j > 0, weight[index2 + MN - M + 1]);
          w[2][1][1] = weight[index2 + 1];
          w[2][1][2] = weight[index2 + MN + 1];
          w[2][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M + 1]);
          w[2][2][2] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[1][0][1] + w[1][0][2] + w[1][1][1] + w[1][1][2] + w[1][2][1] + w[1][2][2] + w[2][0][1] + w[2][0][2] + w[2][1][1] + w[2][1][2] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[1][0][1] > 0)
              result += w[1][0][1] * rhs[index2 - M];
            if (w[1][0][2] > 0)
              result += w[1][0][2] * rhs[index2 + MN - M];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][0][1] > 0)
              result += w[2][0][1] * rhs[index2 - M + 1];
            if (w[2][0][2] > 0)
              result += w[2][0][2] * rhs[index2 + MN - M + 1];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[0][0][1] = 0.25 * VAL(i > 0 && j > 0, weight[index2 - M - 1]);
          w[0][0][2] = 0.25 * VAL(i > 0 && j > 0, weight[index2 + MN - M - 1]);
          w[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          w[0][1][2] = 0.5 * VAL(i > 0, weight[index2 + MN - 1]);
          w[0][2][1] = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + M - 1]);
          w[0][2][2] = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + MN + M - 1]);
          w[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          w[1][0][2] = 0.5 * VAL(j > 0, weight[index2 + MN - M]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = weight[index2 + MN];
          w[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          w[1][2][2] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + MN + M]);
          w[2][0][1] = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 - M + 1]);
          w[2][0][2] = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 + MN - M + 1]);
          w[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          w[2][1][2] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + MN + 1]);
          w[2][2][1] = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + M + 1]);
          w[2][2][2] = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[0][0][1] + w[0][0][2] + w[0][1][1] + w[0][1][2] + w[0][2][1] + w[0][2][2] + w[1][0][1] + w[1][0][2] + w[1][1][1] + w[1][1][2] + w[1][2][1] + w[1][2][2] + w[2][0][1] + w[2][0][2] + w[2][1][1] + w[2][1][2] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[0][0][1] > 0)
              result += w[0][0][1] * rhs[index2 - M - 1];
            if (w[0][0][2] > 0)
              result += w[0][0][2] * rhs[index2 + MN - M - 1];
            if (w[0][1][1] > 0)
              result += w[0][1][1] * rhs[index2 - 1];
            if (w[0][1][2] > 0)
              result += w[0][1][2] * rhs[index2 + MN - 1];
            if (w[0][2][1] > 0)
              result += w[0][2][1] * rhs[index2 + M - 1];
            if (w[0][2][2] > 0)
              result += w[0][2][2] * rhs[index2 + MN + M - 1];
            if (w[1][0][1] > 0)
              result += w[1][0][1] * rhs[index2 - M];
            if (w[1][0][2] > 0)
              result += w[1][0][2] * rhs[index2 + MN - M];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][0][1] > 0)
              result += w[2][0][1] * rhs[index2 - M + 1];
            if (w[2][0][2] > 0)
              result += w[2][0][2] * rhs[index2 + MN - M + 1];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          w[1][2][0] = 0.5 * VAL(p > 0, weight[index2 - MN + M]);
          w[1][2][1] = weight[index2 + M];
          w[1][2][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + M]);
          w[2][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN + 1]);
          w[2][1][1] = weight[index2 + 1];
          w[2][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + 1]);
          w[2][2][0] = 0.5 * VAL(p > 0, weight[index2 - MN + M + 1]);
          w[2][2][1] = weight[index2 + M + 1];
          w[2][2][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[1][1][0] + w[1][1][1] + w[1][1][2] + w[1][2][0] + w[1][2][1] + w[1][2][2] + w[2][1][0] + w[2][1][1] + w[2][1][2] + w[2][2][0] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[1][1][0] > 0)
              result += w[1][1][0] * rhs[index2 - MN];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][0] > 0)
              result += w[1][2][0] * rhs[index2 - MN + M];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][1][0] > 0)
              result += w[2][1][0] * rhs[index2 - MN + 1];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][0] > 0)
              result += w[2][2][0] * rhs[index2 - MN + M + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[0][1][0] = 0.25 * VAL(i > 0 && p > 0, weight[index2 - MN - 1]);
          w[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          w[0][1][2] = 0.25 * VAL(i > 0 && p < Phalf - 1, weight[index2 + MN - 1]);
          w[0][2][0] = 0.25 * VAL(i > 0 && p > 0, weight[index2 - MN + M - 1]);
          w[0][2][1] = 0.5 * VAL(i > 0, weight[index2 + M - 1]);
          w[0][2][2] = 0.25 * VAL(i > 0 && p < Phalf - 1, weight[index2 + MN + M - 1]);
          w[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          w[1][2][0] = 0.5 * VAL(p > 0, weight[index2 - MN + M]);
          w[1][2][1] = weight[index2 + M];
          w[1][2][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + M]);
          w[2][1][0] = 0.25 * VAL(i < Mhalf - 1 && p > 0, weight[index2 - MN + 1]);
          w[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          w[2][1][2] = 0.25 * VAL(i < Mhalf - 1 && p < Phalf - 1, weight[index2 + MN + 1]);
          w[2][2][0] = 0.25 * VAL(i < Mhalf - 1 && p > 0, weight[index2 - MN + M + 1]);
          w[2][2][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + M + 1]);
          w[2][2][2] = 0.25 * VAL(i < Mhalf - 1 && p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[0][1][0] + w[0][1][1] + w[0][1][2] + w[0][2][0] + w[0][2][1] + w[0][2][2] + w[1][1][0] + w[1][1][1] + w[1][1][2] + w[1][2][0] + w[1][2][1] + w[1][2][2] + w[2][1][0] + w[2][1][1] + w[2][1][2] + w[2][2][0] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[0][1][0] > 0)
              result += w[0][1][0] * rhs[index2 - MN - 1];
            if (w[0][1][1] > 0)
              result += w[0][1][1] * rhs[index2 - 1];
            if (w[0][1][2] > 0)
              result += w[0][1][2] * rhs[index2 + MN - 1];
            if (w[0][2][0] > 0)
              result += w[0][2][0] * rhs[index2 - MN + M - 1];
            if (w[0][2][1] > 0)
              result += w[0][2][1] * rhs[index2 + M - 1];
            if (w[0][2][2] > 0)
              result += w[0][2][2] * rhs[index2 + MN + M - 1];
            if (w[1][1][0] > 0)
              result += w[1][1][0] * rhs[index2 - MN];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][0] > 0)
              result += w[1][2][0] * rhs[index2 - MN + M];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][1][0] > 0)
              result += w[2][1][0] * rhs[index2 - MN + 1];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][0] > 0)
              result += w[2][2][0] * rhs[index2 - MN + M + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[1][0][0] = 0.25 * VAL(j > 0 && p > 0, weight[index2 - MN - M]);
          w[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          w[1][0][2] = 0.25 * VAL(j > 0 && p < Phalf - 1, weight[index2 + MN - M]);
          w[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          w[1][2][0] = 0.25 * VAL(j < Nhalf - 1 && p > 0, weight[index2 - MN + M]);
          w[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          w[1][2][2] = 0.25 * VAL(j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M]);
          w[2][0][0] = 0.25 * VAL(j > 0 && p > 0, weight[index2 - MN - M + 1]);
          w[2][0][1] = 0.5 * VAL(j > 0, weight[index2 - M + 1]);
          w[2][0][2] = 0.25 * VAL(j > 0 && p < Phalf - 1, weight[index2 + MN - M + 1]);
          w[2][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN + 1]);
          w[2][1][1] = weight[index2 + 1];
          w[2][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + 1]);
          w[2][2][0] = 0.25 * VAL(j < Nhalf - 1 && p > 0, weight[index2 - MN + M + 1]);
          w[2][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M + 1]);
          w[2][2][2] = 0.25 * VAL(j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[1][0][0] + w[1][0][1] + w[1][0][2] + w[1][1][0] + w[1][1][1] + w[1][1][2] + w[1][2][0] + w[1][2][1] + w[1][2][2] + w[2][0][0] + w[2][0][1] + w[2][0][2] + w[2][1][0] + w[2][1][1] + w[2][1][2] + w[2][2][0] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[1][0][0] > 0)
              result += w[1][0][0] * rhs[index2 - MN - M];
            if (w[1][0][1] > 0)
              result += w[1][0][1] * rhs[index2 - M];
            if (w[1][0][2] > 0)
              result += w[1][0][2] * rhs[index2 + MN - M];
            if (w[1][1][0] > 0)
              result += w[1][1][0] * rhs[index2 - MN];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][0] > 0)
              result += w[1][2][0] * rhs[index2 - MN + M];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][0][0] > 0)
              result += w[2][0][0] * rhs[index2 - MN - M + 1];
            if (w[2][0][1] > 0)
              result += w[2][0][1] * rhs[index2 - M + 1];
            if (w[2][0][2] > 0)
              result += w[2][0][2] * rhs[index2 + MN - M + 1];
            if (w[2][1][0] > 0)
              result += w[2][1][0] * rhs[index2 - MN + 1];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][0] > 0)
              result += w[2][2][0] * rhs[index2 - MN + M + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          
          w[0][0][0] = 0.125 * VAL(i > 0 && j > 0 && p > 0, weight[index2 - MN - M - 1]);
          w[0][0][1] = 0.25 * VAL(i > 0 && j > 0, weight[index2 - M - 1]);
          w[0][0][2] = 0.125 * VAL(i > 0 && j > 0 && p < Phalf - 1, weight[index2 + MN - M - 1]);
          w[0][1][0] = 0.25 * VAL(i > 0 && p > 0, weight[index2 - MN - 1]);
          w[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          w[0][1][2] = 0.25 * VAL(i > 0 && p < Phalf - 1, weight[index2 + MN - 1]);
          w[0][2][0] = 0.125 * VAL(i > 0 && j < Nhalf - 1 && p > 0, weight[index2 - MN + M - 1]);
          w[0][2][1] = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + M - 1]);
          w[0][2][2] = 0.125 * VAL(i > 0 && j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M - 1]);
          w[1][0][0] = 0.25 * VAL(j > 0 && p > 0, weight[index2 - MN - M]);
          w[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          w[1][0][2] = 0.25 * VAL(j > 0 && p < Phalf - 1, weight[index2 + MN - M]);
          w[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          w[1][1][1] = weight[index2];
          w[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          w[1][2][0] = 0.25 * VAL(j < Nhalf - 1 && p > 0, weight[index2 - MN + M]);
          w[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          w[1][2][2] = 0.25 * VAL(j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M]);
          w[2][0][0] = 0.125 * VAL(i < Mhalf - 1 && j > 0 && p > 0, weight[index2 - MN - M + 1]);
          w[2][0][1] = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 - M + 1]);
          w[2][0][2] = 0.125 * VAL(i < Mhalf - 1 && j > 0 && p < Phalf - 1, weight[index2 + MN - M + 1]);
          w[2][1][0] = 0.25 * VAL(i < Mhalf - 1 && p > 0, weight[index2 - MN + 1]);
          w[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          w[2][1][2] = 0.25 * VAL(i < Mhalf - 1 && p < Phalf - 1, weight[index2 + MN + 1]);
          w[2][2][0] = 0.125 * VAL(i < Mhalf - 1 && j < Nhalf - 1 && p > 0, weight[index2 - MN + M + 1]);
          w[2][2][1] = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + M + 1]);
          w[2][2][2] = 0.125 * VAL(i < Mhalf - 1 && j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          sum = w[0][0][0] + w[0][0][1] + w[0][0][2] + w[0][1][0] + w[0][1][1] + w[0][1][2] + w[0][2][0] + w[0][2][1] + w[0][2][2] + w[1][0][0] + w[1][0][1] + w[1][0][2] + w[1][1][0] + w[1][1][1] + w[1][1][2] + w[1][2][0] + w[1][2][1] + w[1][2][2] + w[2][0][0] + w[2][0][1] + w[2][0][2] + w[2][1][0] + w[2][1][1] + w[2][1][2] + w[2][2][0] + w[2][2][1] + w[2][2][2];
          coarse_weight[index1] = sum;
          
          if (sum > 0)
          {
            double result = 0;
            if (w[0][0][0] > 0)
              result += w[0][0][0] * rhs[index2 - MN - M - 1];
            if (w[0][0][1] > 0)
              result += w[0][0][1] * rhs[index2 - M - 1];
            if (w[0][0][2] > 0)
              result += w[0][0][2] * rhs[index2 + MN - M - 1];
            if (w[0][1][0] > 0)
              result += w[0][1][0] * rhs[index2 - MN - 1];
            if (w[0][1][1] > 0)
              result += w[0][1][1] * rhs[index2 - 1];
            if (w[0][1][2] > 0)
              result += w[0][1][2] * rhs[index2 + MN - 1];
            if (w[0][2][0] > 0)
              result += w[0][2][0] * rhs[index2 - MN + M - 1];
            if (w[0][2][1] > 0)
              result += w[0][2][1] * rhs[index2 + M - 1];
            if (w[0][2][2] > 0)
              result += w[0][2][2] * rhs[index2 + MN + M - 1];
            if (w[1][0][0] > 0)
              result += w[1][0][0] * rhs[index2 - MN - M];
            if (w[1][0][1] > 0)
              result += w[1][0][1] * rhs[index2 - M];
            if (w[1][0][2] > 0)
              result += w[1][0][2] * rhs[index2 + MN - M];
            if (w[1][1][0] > 0)
              result += w[1][1][0] * rhs[index2 - MN];
            if (w[1][1][1] > 0)
              result += w[1][1][1] * rhs[index2];
            if (w[1][1][2] > 0)
              result += w[1][1][2] * rhs[index2 + MN];
            if (w[1][2][0] > 0)
              result += w[1][2][0] * rhs[index2 - MN + M];
            if (w[1][2][1] > 0)
              result += w[1][2][1] * rhs[index2 + M];
            if (w[1][2][2] > 0)
              result += w[1][2][2] * rhs[index2 + MN + M];
            if (w[2][0][0] > 0)
              result += w[2][0][0] * rhs[index2 - MN - M + 1];
            if (w[2][0][1] > 0)
              result += w[2][0][1] * rhs[index2 - M + 1];
            if (w[2][0][2] > 0)
              result += w[2][0][2] * rhs[index2 + MN - M + 1];
            if (w[2][1][0] > 0)
              result += w[2][1][0] * rhs[index2 - MN + 1];
            if (w[2][1][1] > 0)
              result += w[2][1][1] * rhs[index2 + 1];
            if (w[2][1][2] > 0)
              result += w[2][1][2] * rhs[index2 + MN + 1];
            if (w[2][2][0] > 0)
              result += w[2][2][0] * rhs[index2 - MN + M + 1];
            if (w[2][2][1] > 0)
              result += w[2][2][1] * rhs[index2 + M + 1];
            if (w[2][2][2] > 0)
              result += w[2][2][2] * rhs[index2 + MN + M + 1];
            
            rhs_coarse[index1] = 8 / sum * result;
          }
        }
      }
    }
  }
}


static void
galerkin3D(int level, int M, int N, int P, int Mhalf, int Nhalf, int Phalf,
           double *weight, double *coarse_weight)
{
  int i, j, p;
  int k;
  double *lhs;
  double *lhs_coarse;
  int MN = M * N;
  int MNhalf = Mhalf * Nhalf;
  double lw[3][3][3];
  double mean;
  
  if (data.lhs[level + 1] != NULL)
    return;
  
  data.lhs[level + 1] = mxCalloc(27 * Mhalf * Nhalf * Phalf,
                                 sizeof(*data.lhs[level + 1]));
  lhs = data.lhs[level];
  lhs_coarse = data.lhs[level + 1];
  
  for (p = 0; p < Phalf; p++)
  {
    for (j = 0; j < Nhalf; j++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        int index1 = ((p * Nhalf + j) * Mhalf + i);
        int index2 = ((2 * p * N + 2 * j) * M + 2 * i);
        double stencil1[3][3][3];
        double stencil2[5][5][5];
        double stencil3[3][3][3];
        double mask1[3][3][3];
        int u, v, w;
        
        for (u = 0; u < 5; u++)
        {
          for (v = 0; v < 5; v++)
          {
            for (w = 0; w < 5; w++)
            {
              stencil2[u][v][w] = 0;
              if (u < 3 && v < 3 && w < 3)
              {
                stencil1[u][v][w] = 0;
                stencil3[u][v][w] = 0;
              }
            }
          }
        }
        mask1[0][0][0] = VAL(i > 0 && j > 0 && p > 0,
                             coarse_weight[index1 - MNhalf - Mhalf - 1]);
        mask1[0][0][1] = VAL(i > 0 && j > 0,
                             coarse_weight[index1 - Mhalf - 1]);
        mask1[0][0][2] = VAL(i > 0 && j > 0 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf - Mhalf - 1]);
        mask1[0][1][0] = VAL(i > 0 && p > 0,
                             coarse_weight[index1 - MNhalf - 1]);
        mask1[0][1][1] = VAL(i > 0,
                             coarse_weight[index1 - 1]);
        mask1[0][1][2] = VAL(i > 0 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf - 1]);
        mask1[0][2][0] = VAL(i > 0 && j < Nhalf - 1 && p > 0,
                             coarse_weight[index1 - MNhalf + Mhalf - 1]);
        mask1[0][2][1] = VAL(i > 0 && j < Nhalf - 1,
                             coarse_weight[index1 + Mhalf - 1]);
        mask1[0][2][2] = VAL(i > 0 && j < Nhalf - 1 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf + Mhalf - 1]);
        mask1[1][0][0] = VAL(j > 0 && p > 0,
                             coarse_weight[index1 - MNhalf - Mhalf]);
        mask1[1][0][1] = VAL(j > 0,
                             coarse_weight[index1 - Mhalf]);
        mask1[1][0][2] = VAL(j > 0 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf - Mhalf]);
        mask1[1][1][0] = VAL(p > 0,
                             coarse_weight[index1 - MNhalf]);
        mask1[1][1][1] = coarse_weight[index1];
        mask1[1][1][2] = VAL(p < Phalf - 1,
                             coarse_weight[index1 + MNhalf]);
        mask1[1][2][0] = VAL(j < Nhalf - 1 && p > 0,
                             coarse_weight[index1 - MNhalf + Mhalf]);
        mask1[1][2][1] = VAL(j < Nhalf - 1,
                             coarse_weight[index1 + Mhalf]);
        mask1[1][2][2] = VAL(j < Nhalf - 1 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf + Mhalf]);
        mask1[2][0][0] = VAL(i < Mhalf - 1 && j > 0 && p > 0,
                             coarse_weight[index1 - MNhalf - Mhalf + 1]);
        mask1[2][0][1] = VAL(i < Mhalf - 1 && j > 0,
                             coarse_weight[index1 - Mhalf + 1]);
        mask1[2][0][2] = VAL(i < Mhalf - 1 && j > 0 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf - Mhalf + 1]);
        mask1[2][1][0] = VAL(i < Mhalf - 1 && p > 0,
                             coarse_weight[index1 - MNhalf + 1]);
        mask1[2][1][1] = VAL(i < Mhalf - 1,
                             coarse_weight[index1 + 1]);
        mask1[2][1][2] = VAL(i < Mhalf - 1 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf + 1]);
        mask1[2][2][0] = VAL(i < Mhalf - 1 && j < Nhalf - 1 && p > 0,
                             coarse_weight[index1 - MNhalf + Mhalf + 1]);
        mask1[2][2][1] = VAL(i < Mhalf - 1 && j < Nhalf - 1,
                             coarse_weight[index1 + Mhalf + 1]);
        mask1[2][2][2] = VAL(i < Mhalf - 1 && j < Nhalf - 1 && p < Phalf - 1,
                             coarse_weight[index1 + MNhalf + Mhalf + 1]);
        
        if (M % 2 == 0 && N % 2 == 0 && P % 2 == 0)
        {
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = weight[index2 + MN];
          lw[1][2][1] = weight[index2 + M];
          lw[1][2][2] = weight[index2 + MN + M];
          lw[2][1][1] = weight[index2 + 1];
          lw[2][1][2] = weight[index2 + MN + 1];
          lw[2][2][1] = weight[index2 + M + 1];
          lw[2][2][2] = weight[index2 + MN + M + 1];
          
          mean = (lw[1][1][1] + lw[1][1][2] + lw[1][2][1] + lw[1][2][2] + lw[2][1][1] + lw[2][1][2] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 1 && N % 2 == 0 && P % 2 == 0)
        {
          lw[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          lw[0][1][2] = 0.5 * VAL(i > 0, weight[index2 + MN - 1]);
          lw[0][2][1] = 0.5 * VAL(i > 0, weight[index2 + M - 1]);
          lw[0][2][2] = 0.5 * VAL(i > 0, weight[index2 + MN + M - 1]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = weight[index2 + MN];
          lw[1][2][1] = weight[index2 + M];
          lw[1][2][2] = weight[index2 + MN + M];
          lw[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          lw[2][1][2] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + MN + 1]);
          lw[2][2][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + M + 1]);
          lw[2][2][2] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[0][1][1] + lw[0][1][2] + lw[0][2][1] + lw[0][2][2] + lw[1][1][1] + lw[1][1][2] + lw[1][2][1] + lw[1][2][2] + lw[2][1][1] + lw[2][1][2] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[0][1][1] = lw[0][1][1] / mean;
          stencil1[0][1][2] = lw[0][1][2] / mean;
          stencil1[0][2][1] = lw[0][2][1] / mean;
          stencil1[0][2][2] = lw[0][2][2] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 0 && N % 2 == 1 && P % 2 == 0)
        {
          lw[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          lw[1][0][2] = 0.5 * VAL(j > 0, weight[index2 + MN - M]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = weight[index2 + MN];
          lw[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          lw[1][2][2] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + MN + M]);
          lw[2][0][1] = 0.5 * VAL(j > 0, weight[index2 - M + 1]);
          lw[2][0][2] = 0.5 * VAL(j > 0, weight[index2 + MN - M + 1]);
          lw[2][1][1] = weight[index2 + 1];
          lw[2][1][2] = weight[index2 + MN + 1];
          lw[2][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M + 1]);
          lw[2][2][2] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[1][0][1] + lw[1][0][2] + lw[1][1][1] + lw[1][1][2] + lw[1][2][1] + lw[1][2][2] + lw[2][0][1] + lw[2][0][2] + lw[2][1][1] + lw[2][1][2] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[1][0][1] = lw[1][0][1] / mean;
          stencil1[1][0][2] = lw[1][0][2] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][0][1] = lw[2][0][1] / mean;
          stencil1[2][0][2] = lw[2][0][2] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 1 && N % 2 == 1 && P % 2 == 0)
        {
          lw[0][0][1] = 0.25 * VAL(i > 0 && j > 0, weight[index2 - M - 1]);
          lw[0][0][2] = 0.25 * VAL(i > 0 && j > 0, weight[index2 + MN - M - 1]);
          lw[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          lw[0][1][2] = 0.5 * VAL(i > 0, weight[index2 + MN - 1]);
          lw[0][2][1] = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + M - 1]);
          lw[0][2][2] = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + MN + M - 1]);
          lw[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          lw[1][0][2] = 0.5 * VAL(j > 0, weight[index2 + MN - M]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = weight[index2 + MN];
          lw[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          lw[1][2][2] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + MN + M]);
          lw[2][0][1] = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 - M + 1]);
          lw[2][0][2] = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 + MN - M + 1]);
          lw[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          lw[2][1][2] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + MN + 1]);
          lw[2][2][1] = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + M + 1]);
          lw[2][2][2] = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[0][0][1] + lw[0][0][2] + lw[0][1][1] + lw[0][1][2] + lw[0][2][1] + lw[0][2][2] + lw[1][0][1] + lw[1][0][2] + lw[1][1][1] + lw[1][1][2] + lw[1][2][1] + lw[1][2][2] + lw[2][0][1] + lw[2][0][2] + lw[2][1][1] + lw[2][1][2] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[0][0][1] = lw[0][0][1] / mean;
          stencil1[0][0][2] = lw[0][0][2] / mean;
          stencil1[0][1][1] = lw[0][1][1] / mean;
          stencil1[0][1][2] = lw[0][1][2] / mean;
          stencil1[0][2][1] = lw[0][2][1] / mean;
          stencil1[0][2][2] = lw[0][2][2] / mean;
          stencil1[1][0][1] = lw[1][0][1] / mean;
          stencil1[1][0][2] = lw[1][0][2] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][0][1] = lw[2][0][1] / mean;
          stencil1[2][0][2] = lw[2][0][2] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 0 && N % 2 == 0 && P % 2 == 1)
        {
          lw[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          lw[1][2][0] = 0.5 * VAL(p > 0, weight[index2 - MN + M]);
          lw[1][2][1] = weight[index2 + M];
          lw[1][2][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + M]);
          lw[2][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN + 1]);
          lw[2][1][1] = weight[index2 + 1];
          lw[2][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + 1]);
          lw[2][2][0] = 0.5 * VAL(p > 0, weight[index2 - MN + M + 1]);
          lw[2][2][1] = weight[index2 + M + 1];
          lw[2][2][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[1][1][0] + lw[1][1][1] + lw[1][1][2] + lw[1][2][0] + lw[1][2][1] + lw[1][2][2] + lw[2][1][0] + lw[2][1][1] + lw[2][1][2] + lw[2][2][0] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[1][1][0] = lw[1][1][0] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][0] = lw[1][2][0] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][1][0] = lw[2][1][0] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][0] = lw[2][2][0] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 1 && N % 2 == 0 && P % 2 == 1)
        {
          lw[0][1][0] = 0.25 * VAL(i > 0 && p > 0, weight[index2 - MN - 1]);
          lw[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          lw[0][1][2] = 0.25 * VAL(i > 0 && p < Phalf - 1, weight[index2 + MN - 1]);
          lw[0][2][0] = 0.25 * VAL(i > 0 && p > 0, weight[index2 - MN + M - 1]);
          lw[0][2][1] = 0.5 * VAL(i > 0, weight[index2 + M - 1]);
          lw[0][2][2] = 0.25 * VAL(i > 0 && p < Phalf - 1, weight[index2 + MN + M - 1]);
          lw[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          lw[1][2][0] = 0.5 * VAL(p > 0, weight[index2 - MN + M]);
          lw[1][2][1] = weight[index2 + M];
          lw[1][2][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + M]);
          lw[2][1][0] = 0.25 * VAL(i < Mhalf - 1 && p > 0, weight[index2 - MN + 1]);
          lw[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          lw[2][1][2] = 0.25 * VAL(i < Mhalf - 1 && p < Phalf - 1, weight[index2 + MN + 1]);
          lw[2][2][0] = 0.25 * VAL(i < Mhalf - 1 && p > 0, weight[index2 - MN + M + 1]);
          lw[2][2][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + M + 1]);
          lw[2][2][2] = 0.25 * VAL(i < Mhalf - 1 && p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[0][1][0] + lw[0][1][1] + lw[0][1][2] + lw[0][2][0] + lw[0][2][1] + lw[0][2][2] + lw[1][1][0] + lw[1][1][1] + lw[1][1][2] + lw[1][2][0] + lw[1][2][1] + lw[1][2][2] + lw[2][1][0] + lw[2][1][1] + lw[2][1][2] + lw[2][2][0] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[0][1][0] = lw[0][1][0] / mean;
          stencil1[0][1][1] = lw[0][1][1] / mean;
          stencil1[0][1][2] = lw[0][1][2] / mean;
          stencil1[0][2][0] = lw[0][2][0] / mean;
          stencil1[0][2][1] = lw[0][2][1] / mean;
          stencil1[0][2][2] = lw[0][2][2] / mean;
          stencil1[1][1][0] = lw[1][1][0] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][0] = lw[1][2][0] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][1][0] = lw[2][1][0] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][0] = lw[2][2][0] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 0 && N % 2 == 1 && P % 2 == 1)
        {
          lw[1][0][0] = 0.25 * VAL(j > 0 && p > 0, weight[index2 - MN - M]);
          lw[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          lw[1][0][2] = 0.25 * VAL(j > 0 && p < Phalf - 1, weight[index2 + MN - M]);
          lw[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          lw[1][2][0] = 0.25 * VAL(j < Nhalf - 1 && p > 0, weight[index2 - MN + M]);
          lw[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          lw[1][2][2] = 0.25 * VAL(j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M]);
          lw[2][0][0] = 0.25 * VAL(j > 0 && p > 0, weight[index2 - MN - M + 1]);
          lw[2][0][1] = 0.5 * VAL(j > 0, weight[index2 - M + 1]);
          lw[2][0][2] = 0.25 * VAL(j > 0 && p < Phalf - 1, weight[index2 + MN - M + 1]);
          lw[2][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN + 1]);
          lw[2][1][1] = weight[index2 + 1];
          lw[2][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN + 1]);
          lw[2][2][0] = 0.25 * VAL(j < Nhalf - 1 && p > 0, weight[index2 - MN + M + 1]);
          lw[2][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M + 1]);
          lw[2][2][2] = 0.25 * VAL(j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[1][0][0] + lw[1][0][1] + lw[1][0][2] + lw[1][1][0] + lw[1][1][1] + lw[1][1][2] + lw[1][2][0] + lw[1][2][1] + lw[1][2][2] + lw[2][0][0] + lw[2][0][1] + lw[2][0][2] + lw[2][1][0] + lw[2][1][1] + lw[2][1][2] + lw[2][2][0] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[1][0][0] = lw[1][0][0] / mean;
          stencil1[1][0][1] = lw[1][0][1] / mean;
          stencil1[1][0][2] = lw[1][0][2] / mean;
          stencil1[1][1][0] = lw[1][1][0] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][0] = lw[1][2][0] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][0][0] = lw[2][0][0] / mean;
          stencil1[2][0][1] = lw[2][0][1] / mean;
          stencil1[2][0][2] = lw[2][0][2] / mean;
          stencil1[2][1][0] = lw[2][1][0] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][0] = lw[2][2][0] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        if (M % 2 == 1 && N % 2 == 1 && P % 2 == 1)
        {
          lw[0][0][0] = 0.125 * VAL(i > 0 && j > 0 && p > 0, weight[index2 - MN - M - 1]);
          lw[0][0][1] = 0.25 * VAL(i > 0 && j > 0, weight[index2 - M - 1]);
          lw[0][0][2] = 0.125 * VAL(i > 0 && j > 0 && p < Phalf - 1, weight[index2 + MN - M - 1]);
          lw[0][1][0] = 0.25 * VAL(i > 0 && p > 0, weight[index2 - MN - 1]);
          lw[0][1][1] = 0.5 * VAL(i > 0, weight[index2 - 1]);
          lw[0][1][2] = 0.25 * VAL(i > 0 && p < Phalf - 1, weight[index2 + MN - 1]);
          lw[0][2][0] = 0.125 * VAL(i > 0 && j < Nhalf - 1 && p > 0, weight[index2 - MN + M - 1]);
          lw[0][2][1] = 0.25 * VAL(i > 0 && j < Nhalf - 1, weight[index2 + M - 1]);
          lw[0][2][2] = 0.125 * VAL(i > 0 && j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M - 1]);
          lw[1][0][0] = 0.25 * VAL(j > 0 && p > 0, weight[index2 - MN - M]);
          lw[1][0][1] = 0.5 * VAL(j > 0, weight[index2 - M]);
          lw[1][0][2] = 0.25 * VAL(j > 0 && p < Phalf - 1, weight[index2 + MN - M]);
          lw[1][1][0] = 0.5 * VAL(p > 0, weight[index2 - MN]);
          lw[1][1][1] = weight[index2];
          lw[1][1][2] = 0.5 * VAL(p < Phalf - 1, weight[index2 + MN]);
          lw[1][2][0] = 0.25 * VAL(j < Nhalf - 1 && p > 0, weight[index2 - MN + M]);
          lw[1][2][1] = 0.5 * VAL(j < Nhalf - 1, weight[index2 + M]);
          lw[1][2][2] = 0.25 * VAL(j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M]);
          lw[2][0][0] = 0.125 * VAL(i < Mhalf - 1 && j > 0 && p > 0, weight[index2 - MN - M + 1]);
          lw[2][0][1] = 0.25 * VAL(i < Mhalf - 1 && j > 0, weight[index2 - M + 1]);
          lw[2][0][2] = 0.125 * VAL(i < Mhalf - 1 && j > 0 && p < Phalf - 1, weight[index2 + MN - M + 1]);
          lw[2][1][0] = 0.25 * VAL(i < Mhalf - 1 && p > 0, weight[index2 - MN + 1]);
          lw[2][1][1] = 0.5 * VAL(i < Mhalf - 1, weight[index2 + 1]);
          lw[2][1][2] = 0.25 * VAL(i < Mhalf - 1 && p < Phalf - 1, weight[index2 + MN + 1]);
          lw[2][2][0] = 0.125 * VAL(i < Mhalf - 1 && j < Nhalf - 1 && p > 0, weight[index2 - MN + M + 1]);
          lw[2][2][1] = 0.25 * VAL(i < Mhalf - 1 && j < Nhalf - 1, weight[index2 + M + 1]);
          lw[2][2][2] = 0.125 * VAL(i < Mhalf - 1 && j < Nhalf - 1 && p < Phalf - 1, weight[index2 + MN + M + 1]);
          
          mean = (lw[0][0][0] + lw[0][0][1] + lw[0][0][2] + lw[0][1][0] + lw[0][1][1] + lw[0][1][2] + lw[0][2][0] + lw[0][2][1] + lw[0][2][2] + lw[1][0][0] + lw[1][0][1] + lw[1][0][2] + lw[1][1][0] + lw[1][1][1] + lw[1][1][2] + lw[1][2][0] + lw[1][2][1] + lw[1][2][2] + lw[2][0][0] + lw[2][0][1] + lw[2][0][2] + lw[2][1][0] + lw[2][1][1] + lw[2][1][2] + lw[2][2][0] + lw[2][2][1] + lw[2][2][2]) / 8;
          
          if (mean == 0)
            continue;
          
          stencil1[0][0][0] = lw[0][0][0] / mean;
          stencil1[0][0][1] = lw[0][0][1] / mean;
          stencil1[0][0][2] = lw[0][0][2] / mean;
          stencil1[0][1][0] = lw[0][1][0] / mean;
          stencil1[0][1][1] = lw[0][1][1] / mean;
          stencil1[0][1][2] = lw[0][1][2] / mean;
          stencil1[0][2][0] = lw[0][2][0] / mean;
          stencil1[0][2][1] = lw[0][2][1] / mean;
          stencil1[0][2][2] = lw[0][2][2] / mean;
          stencil1[1][0][0] = lw[1][0][0] / mean;
          stencil1[1][0][1] = lw[1][0][1] / mean;
          stencil1[1][0][2] = lw[1][0][2] / mean;
          stencil1[1][1][0] = lw[1][1][0] / mean;
          stencil1[1][1][1] = lw[1][1][1] / mean;
          stencil1[1][1][2] = lw[1][1][2] / mean;
          stencil1[1][2][0] = lw[1][2][0] / mean;
          stencil1[1][2][1] = lw[1][2][1] / mean;
          stencil1[1][2][2] = lw[1][2][2] / mean;
          stencil1[2][0][0] = lw[2][0][0] / mean;
          stencil1[2][0][1] = lw[2][0][1] / mean;
          stencil1[2][0][2] = lw[2][0][2] / mean;
          stencil1[2][1][0] = lw[2][1][0] / mean;
          stencil1[2][1][1] = lw[2][1][1] / mean;
          stencil1[2][1][2] = lw[2][1][2] / mean;
          stencil1[2][2][0] = lw[2][2][0] / mean;
          stencil1[2][2][1] = lw[2][2][1] / mean;
          stencil1[2][2][2] = lw[2][2][2] / mean;
        }
        
        for (u = 0; u < 3; u++)
        {
          for (v = 0; v < 3; v++)
          {
            for (w = 0; w < 3; w++)
            {
              if (stencil1[u][v][w] != 0)
              {
                int index = 27 * (index2 + ((w - 1) * N + v - 1) * M + u - 1);
                if (lhs[index + 13] != 0)
                {
                  int k;
                  for (k = 0; k < 27; k++)
                  {
                    int a = (k % 3);
                    int b = ((k / 3) % 3);
                    int c = ((k / 9) % 3);
                    stencil2[u + a][v + b][w + c] += stencil1[u][v][w] * lhs[index + k];
                  }
                }
              }
            }
          }
        }
        
        if (M % 2 == 0 && N % 2 == 0 && P % 2 == 0)
        {
          for (u = 1; u < 5; u++)
          {
            for (v = 1; v < 5; v++)
            {
              for (w = 1; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0.75;
                else
                  alpha1 = 0.25;
                
                if (v % 2 == 0)
                  alpha2 = 0.75;
                else
                  alpha2 = 0.25;
                
                if (w % 2 == 0)
                  alpha3 = 0.75;
                else
                  alpha3 = 0.25;
                
                uu = (u - 1) / 2;
                vv = (v - 1) / 2;
                ww = (w - 1) / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 1 && N % 2 == 0 && P % 2 == 0)
        {
          for (u = 0; u < 5; u++)
          {
            for (v = 1; v < 5; v++)
            {
              for (w = 1; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0;
                else
                  alpha1 = 0.5;
                
                if (v % 2 == 0)
                  alpha2 = 0.75;
                else
                  alpha2 = 0.25;
                
                if (w % 2 == 0)
                  alpha3 = 0.75;
                else
                  alpha3 = 0.25;
                
                uu = u / 2;
                vv = (v - 1) / 2;
                ww = (w - 1) / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 0 && N % 2 == 1 && P % 2 == 0)
        {
          for (u = 1; u < 5; u++)
          {
            for (v = 0; v < 5; v++)
            {
              for (w = 1; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0.75;
                else
                  alpha1 = 0.25;
                
                if (v % 2 == 0)
                  alpha2 = 0;
                else
                  alpha2 = 0.5;
                
                if (w % 2 == 0)
                  alpha3 = 0.75;
                else
                  alpha3 = 0.25;
                
                uu = (u - 1) / 2;
                vv = v / 2;
                ww = (w - 1) / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 1 && N % 2 == 1 && P % 2 == 0)
        {
          for (u = 0; u < 5; u++)
          {
            for (v = 0; v < 5; v++)
            {
              for (w = 1; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0;
                else
                  alpha1 = 0.5;
                
                if (v % 2 == 0)
                  alpha2 = 0;
                else
                  alpha2 = 0.5;
                
                if (w % 2 == 0)
                  alpha3 = 0.75;
                else
                  alpha3 = 0.25;
                
                uu = u / 2;
                vv = v / 2;
                ww = (w - 1) / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 0 && N % 2 == 0 && P % 2 == 1)
        {
          for (u = 1; u < 5; u++)
          {
            for (v = 1; v < 5; v++)
            {
              for (w = 0; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0.75;
                else
                  alpha1 = 0.25;
                
                if (v % 2 == 0)
                  alpha2 = 0.75;
                else
                  alpha2 = 0.25;
                
                if (w % 2 == 0)
                  alpha3 = 0;
                else
                  alpha3 = 0.5;
                
                uu = (u - 1) / 2;
                vv = (v - 1) / 2;
                ww = w / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 1 && N % 2 == 0 && P % 2 == 1)
        {
          for (u = 0; u < 5; u++)
          {
            for (v = 1; v < 5; v++)
            {
              for (w = 0; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0;
                else
                  alpha1 = 0.5;
                
                if (v % 2 == 0)
                  alpha2 = 0.75;
                else
                  alpha2 = 0.25;
                
                if (w % 2 == 0)
                  alpha3 = 0;
                else
                  alpha3 = 0.5;
                
                uu = u / 2;
                vv = (v - 1) / 2;
                ww = w / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 0 && N % 2 == 1 && P % 2 == 1)
        {
          for (u = 1; u < 5; u++)
          {
            for (v = 0; v < 5; v++)
            {
              for (w = 0; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0.75;
                else
                  alpha1 = 0.25;
                
                if (v % 2 == 0)
                  alpha2 = 0;
                else
                  alpha2 = 0.5;
                
                if (w % 2 == 0)
                  alpha3 = 0;
                else
                  alpha3 = 0.5;
                
                uu = (u - 1) / 2;
                vv = v / 2;
                ww = w / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        if (M % 2 == 1 && N % 2 == 1 && P % 2 == 1)
        {
          for (u = 0; u < 5; u++)
          {
            for (v = 0; v < 5; v++)
            {
              for (w = 0; w < 5; w++)
              {
                double alpha1, alpha2, alpha3;
                double unw, dnw, une, dne, usw, dsw, use, dse;
                int uu, vv, ww;
                double sum;
                
                if (u % 2 == 0)
                  alpha1 = 0;
                else
                  alpha1 = 0.5;
                
                if (v % 2 == 0)
                  alpha2 = 0;
                else
                  alpha2 = 0.5;
                
                if (w % 2 == 0)
                  alpha3 = 0;
                else
                  alpha3 = 0.5;
                
                uu = u / 2;
                vv = v / 2;
                ww = w / 2;
                unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                       mask1[uu][vv][ww]);
                dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                       VAL(alpha3 > 0,
                           mask1[uu][vv][ww + 1]));
                une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                       VAL(alpha2 > 0,
                           mask1[uu][vv + 1][ww]));
                dne = ((1 - alpha1) * alpha2 * alpha3 *
                       VAL(alpha2 > 0 && alpha3 > 0,
                           mask1[uu][vv + 1][ww + 1]));
                usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                       VAL(alpha1 > 0,
                           mask1[uu + 1][vv][ww]));
                dsw = (alpha1 * (1 - alpha2) * alpha3 *
                       VAL(alpha1 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv][ww + 1]));
                use = (alpha1 * alpha2 * (1 - alpha3) *
                       VAL(alpha1 > 0 && alpha2 > 0,
                           mask1[uu + 1][vv + 1][ww]));
                dse = (alpha1 * alpha2 * alpha3 *
                       VAL(alpha1 > 0 && alpha2 > 0 && alpha3 > 0,
                           mask1[uu + 1][vv + 1][ww + 1]));
                
                sum = unw + dnw + une + dne + usw + dsw + use + dse;
                
                if (sum > 0)
                {
                  if (unw > 0)
                    stencil3[uu][vv][ww] += unw * stencil2[u][v][w] / sum;
                  if (dnw > 0)
                    stencil3[uu][vv][ww + 1] += dnw * stencil2[u][v][w] / sum;
                  if (une > 0)
                    stencil3[uu][vv + 1][ww] += une * stencil2[u][v][w] / sum;
                  if (dne > 0)
                    stencil3[uu][vv + 1][ww + 1] += dne * stencil2[u][v][w] / sum;
                  if (usw > 0)
                    stencil3[uu + 1][vv][ww] += usw * stencil2[u][v][w] / sum;
                  if (dsw > 0)
                    stencil3[uu + 1][vv][ww + 1] += dsw * stencil2[u][v][w] / sum;
                  if (use > 0)
                    stencil3[uu + 1][vv + 1][ww] += use * stencil2[u][v][w] / sum;
                  if (dse > 0)
                    stencil3[uu + 1][vv + 1][ww + 1] += dse * stencil2[u][v][w] / sum;
                }
              }
            }
          }
        }
        
        for (k = 0; k < 27; k++)
        {
          int a = (k % 3);
          int b = ((k / 3) % 3);
          int c = ((k / 9) % 3);
          lhs_coarse[27 * index1 + k] = stencil3[a][b][c];
        }
        
        for (u = 0; u < 27; u++)
        {
          if (u != 13 && lhs_coarse[27 * index1 + u] < 0)
          {
            lhs_coarse[27 * index1 + 13] += lhs_coarse[27 * index1 + u];
            lhs_coarse[27 * index1 + u] = 0;
          }
        }
      }
    }
  }
}


static void
upsample3D(int M, int N, int P,
           double *v, int Mhalf, int Nhalf, int Phalf,
           double *f_out, double *weight, double *coarse_weight)
{
  int i, j, p;
  int index1;
  int index2;
  int MNhalf = Mhalf * Nhalf;
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = (((p + 1) / 2 - 1) * Nhalf + (j + 1) / 2 - 1) * Mhalf + (i + 1) / 2 - 1;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0.75;
          else
            alpha1 = 0.25;
          
          if (j % 2 == 0)
            alpha2 = 0.75;
          else
            alpha2 = 0.25;
          
          if (p % 2 == 0)
            alpha3 = 0.75;
          else
            alpha3 = 0.25;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(i > 0 && j > 0 && p > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(i > 0 && j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(i > 0 && j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(i > 0 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1 && j > 0 && p > 0,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = (((p + 1) / 2 - 1) * Nhalf + (j + 1) / 2 - 1) * Mhalf + i / 2;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0;
          else
            alpha1 = 0.5;
          
          if (j % 2 == 0)
            alpha2 = 0.75;
          else
            alpha2 = 0.25;
          
          if (p % 2 == 0)
            alpha3 = 0.75;
          else
            alpha3 = 0.25;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(j > 0 && p > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1 && j > 0 && p > 0,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = (((p + 1) / 2 - 1) * Nhalf + j / 2) * Mhalf + (i + 1) / 2 - 1;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0.75;
          else
            alpha1 = 0.25;
          
          if (j % 2 == 0)
            alpha2 = 0;
          else
            alpha2 = 0.5;
          
          if (p % 2 == 0)
            alpha3 = 0.75;
          else
            alpha3 = 0.25;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(i > 0 && p > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(i > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(i > 0 && j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(i > 0 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1 && p > 0,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = (((p + 1) / 2 - 1) * Nhalf + j / 2) * Mhalf + i / 2;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0;
          else
            alpha1 = 0.5;
          
          if (j % 2 == 0)
            alpha2 = 0;
          else
            alpha2 = 0.5;
          
          if (p % 2 == 0)
            alpha3 = 0.75;
          else
            alpha3 = 0.25;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(p > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1 && p > 0,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1 && p > 0,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 1)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = ((p / 2) * Nhalf + (j + 1) / 2 - 1) * Mhalf + (i + 1) / 2 - 1;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0.75;
          else
            alpha1 = 0.25;
          
          if (j % 2 == 0)
            alpha2 = 0.75;
          else
            alpha2 = 0.25;
          
          if (p % 2 == 0)
            alpha3 = 0;
          else
            alpha3 = 0.5;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(i > 0 && j > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(i > 0 && j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(i > 0 && j < N - 1,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(i > 0 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1 && j > 0,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 1)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = ((p / 2) * Nhalf + (j + 1) / 2 - 1) * Mhalf + i / 2;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0;
          else
            alpha1 = 0.5;
          
          if (j % 2 == 0)
            alpha2 = 0.75;
          else
            alpha2 = 0.25;
          
          if (p % 2 == 0)
            alpha3 = 0;
          else
            alpha3 = 0.5;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(j > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(j < N - 1,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1 && j > 0,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && j > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 1)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = ((p / 2) * Nhalf + j / 2) * Mhalf + (i + 1) / 2 - 1;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0.75;
          else
            alpha1 = 0.25;
          
          if (j % 2 == 0)
            alpha2 = 0;
          else
            alpha2 = 0.5;
          
          if (p % 2 == 0)
            alpha3 = 0;
          else
            alpha3 = 0.5;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 VAL(i > 0,
                     coarse_weight[index2]));
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(i > 0 && p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(i > 0 && j < N - 1,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(i > 0 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 1)
  {
    for (p = 0; p < P; p++)
    {
      for (j = 0; j < N; j++)
      {
        for (i = 0; i < M; i++)
        {
          double alpha1, alpha2, alpha3;
          double unw, dnw, une, dne, usw, dsw, use, dse;
          double sum;
          
          index1 = (p * N + j) * M + i;
          index2 = ((p / 2) * Nhalf + j / 2) * Mhalf + i / 2;
          
          if (weight[index1] == 0)
            continue;
          
          if (i % 2 == 0)
            alpha1 = 0;
          else
            alpha1 = 0.5;
          
          if (j % 2 == 0)
            alpha2 = 0;
          else
            alpha2 = 0.5;
          
          if (p % 2 == 0)
            alpha3 = 0;
          else
            alpha3 = 0.5;
          
          unw = ((1 - alpha1) * (1 - alpha2) * (1 - alpha3) *
                 coarse_weight[index2]);
          dnw = ((1 - alpha1) * (1 - alpha2) * alpha3 *
                 VAL(p < P - 1,
                     coarse_weight[index2 + MNhalf]));
          une = ((1 - alpha1) * alpha2 * (1 - alpha3) *
                 VAL(j < N - 1,
                     coarse_weight[index2 + Mhalf]));
          dne = ((1 - alpha1) * alpha2 * alpha3 *
                 VAL(j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf]));
          usw = (alpha1 * (1 - alpha2) * (1 - alpha3) *
                 VAL(i < M - 1,
                     coarse_weight[index2 + 1]));
          dsw = (alpha1 * (1 - alpha2) * alpha3 *
                 VAL(i < M - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + 1]));
          use = (alpha1 * alpha2 * (1 - alpha3) *
                 VAL(i < M - 1 && j < N - 1,
                     coarse_weight[index2 + Mhalf + 1]));
          dse = (alpha1 * alpha2 * alpha3 *
                 VAL(i < M - 1 && j < N - 1 && p < P - 1,
                     coarse_weight[index2 + MNhalf + Mhalf + 1]));
          
          sum = unw + dnw + une + dne + usw + dsw + use + dse;
          
          if (sum > 0)
          {
            double contribution = 0;
            
            if (unw > 0)
              contribution += unw * v[index2];
            if (dnw > 0)
              contribution += dnw * v[index2 + MNhalf];
            if (une > 0)
              contribution += une * v[index2 + Mhalf];
            if (dne > 0)
              contribution += dne * v[index2 + MNhalf + Mhalf];
            if (usw > 0)
              contribution += usw * v[index2 + 1];
            if (dsw > 0)
              contribution += dsw * v[index2 + MNhalf + 1];
            if (use > 0)
              contribution += use * v[index2 + Mhalf + 1];
            if (dse > 0)
              contribution += dse * v[index2 + MNhalf + Mhalf + 1];
            
            f_out[index1] += contribution / sum;
          }
        }
      }
    }
  }
}


/* Recursive multigrid function.*/
static void
poisson_multigrid3D(int level, double *rhs, double *weight,
		    int n1, int n2, int nm,
		    double *f_out,
		    int M, int N, int P, int *directly_solved)
{
  int k;
  double *r;
  double *r_downsampled;
  double *coarse_weight;
  double *lhs = data.lhs[level];
  double *v;
  int Mhalf;
  int Nhalf;
  int Phalf;
  int MN = M * N;
  
  /* Solve a sufficiently small problem directly. */
  if (M < RECURSION_SIZE_LIMIT
      || N < RECURSION_SIZE_LIMIT
      || P < RECURSION_SIZE_LIMIT)
  {
    solve_directly3D(lhs, rhs, f_out, M, N, P);
    *directly_solved = 1;
    return;
  }
  *directly_solved = 0;
  
  /* Pre-smoothing. */
  for (k = 0; k < n1; k++)
    gauss_seidel3D(f_out, lhs, rhs, M, N, P);
  
  /* Compute residual. */
  r = mxCalloc(M * N * P, sizeof(*r));
  compute_residual3D(r, lhs, rhs, f_out, M, N, P);
  
  /* Downsample residual. */
  Mhalf = (M + 1) / 2;
  Nhalf = (N + 1) / 2;
  Phalf = (P + 1) / 2;
  r_downsampled = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*r_downsampled));
  coarse_weight = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*coarse_weight));
  downsample3D(r, M, N, P, r_downsampled, Mhalf, Nhalf, Phalf,
	       weight, coarse_weight);
  galerkin3D(level, M, N, P, Mhalf, Nhalf, Phalf, weight, coarse_weight);

  /* Recurse to compute a correction. */
  v = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*v));
  for (k = 0; k < nm; k++)
  {
    int directly_solved;
    poisson_multigrid3D(level + 1, r_downsampled, coarse_weight,
			n1, n2, nm, v, Mhalf, Nhalf, Phalf, &directly_solved);
    if (directly_solved)
      break;
  }
  
  upsample3D(M, N, P, v, Mhalf, Nhalf, Phalf, f_out, weight, coarse_weight);
  
  /* Post-smoothing. */
  for (k = 0; k < n2; k++)
    gauss_seidel3D(f_out, lhs, rhs, M, N, P);
  
  /* Set the mean value to zero.
   *
   * FIXME: This should not be needed (I believe) and might indicate
   * some bug elsewhere.
   */
  if (0)
  {
    double sum = 0.0;
    int num_samples_in_mask = 0;
    double mean;
    int i;
    
    for (i = 0; i < M * N * P; i++)
      if (weight[i] != 0)
      {
	sum += f_out[i];
	num_samples_in_mask++;
      }
    
    mean = sum / num_samples_in_mask;
    for (i = 0; i < M * N * P; i++)
      if (weight[i] != 0)
	f_out[i] -= mean;
  }
  
  mxFree(r);
  mxFree(r_downsampled);
  mxFree(coarse_weight);
  mxFree(v);
}


/* It is assumed that f_out is initialized to zero when called. */
static void
poisson_full_multigrid3D(int level, double *rhs, double *weight,
			 int number_of_iterations, int M, int N, int P,
			 double *f_out)
{
  double *rhs_downsampled;
  double *coarse_weight;
  double *f_coarse;
  int k;
  
  /* Unless already coarsest scale, first recurse to coarser scale. */
  if (M >= RECURSION_SIZE_LIMIT
      && N >= RECURSION_SIZE_LIMIT
      && P >= RECURSION_SIZE_LIMIT)
  {
    /* Downsample right hand side. */
    int Mhalf = (M + 1) / 2;
    int Nhalf = (N + 1) / 2;
    int Phalf = (P + 1) / 2;
    rhs_downsampled = mxCalloc(Mhalf * Nhalf * Phalf,
			       sizeof(*rhs_downsampled));
    coarse_weight = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*coarse_weight));
    downsample3D(rhs, M, N, P, rhs_downsampled, Mhalf, Nhalf, Phalf,
		 weight, coarse_weight);
    galerkin3D(level, M, N, P, Mhalf, Nhalf, Phalf, weight, coarse_weight);

    f_coarse = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*f_coarse));
    poisson_full_multigrid3D(level + 1, rhs_downsampled, coarse_weight,
			     number_of_iterations,
			     Mhalf, Nhalf, Phalf, f_coarse);
    /* Upsample the coarse result. */
    upsample3D(M, N, P, f_coarse, Mhalf, Nhalf, Phalf, f_out,
	       weight, coarse_weight);

    mxFree(f_coarse);
    mxFree(coarse_weight);
    mxFree(rhs_downsampled);
  }
  
  /* Perform number_of_iterations standard multigrid cycles. */
  for (k = 0; k < number_of_iterations; k++)
  {
    int directly_solved;
    poisson_multigrid3D(level, rhs, weight, 2, 2, 2, f_out, M, N, P,
			&directly_solved);
    if (directly_solved)
      break;
  }
}


static void
antigradient3D(double *g, double *mask, double mu, int number_of_iterations,
	       int M, int N, int P, double *f_out)
{
  double *rhs;
  double *lhs;
  double *weight;
  double sum;
  double mean;
  int i, j, p;
  int num_samples_in_mask;
  int MN = M * N;
  
  clear_global_data();

  /* Compute left and right hand sides of Poisson problem with Neumann
   * boundary conditions, discretized by finite differences.
   */
  rhs = mxCalloc(M * N * P, sizeof(*rhs));
  lhs = mxCalloc(27 * M * N * P, sizeof(*lhs));
  data.lhs[0] = lhs;
  weight = mxCalloc(M * N * P, sizeof(*weight));
  for (p = 0; p < P; p++)
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	int index1 = (p * N + j) * M + i;
	int index2 = index1 + M * N * P;
	int index3 = index1 + 2 * M * N * P;
	double d = 0.0;
	int N_missing = 0;
	int S_missing = 0;
	int W_missing = 0;
	int E_missing = 0;
	int U_missing = 0;
	int D_missing = 0;
	
	if (mask && mask[index1] == 0)
	  continue;
	
	weight[index1] = 1;
	
	if (i == 0 || (mask && mask[index1 - 1] == 0))
	  N_missing = 1;
	
	if (i == M - 1 || (mask && mask[index1 + 1] == 0))
	  S_missing = 1;
	
	if (j == 0 || (mask && mask[index1 - M] == 0))
	  W_missing = 1;
	
	if (j == N - 1 || (mask && mask[index1 + M] == 0))
	  E_missing = 1;
	
	if (p == 0 || (mask && mask[index1 - MN] == 0))
	  U_missing = 1;
	
	if (p == P - 1 || (mask && mask[index1 + MN] == 0))
	  D_missing = 1;
	
	if (N_missing && !S_missing)
	  d = g[index1 + 1] + g[index1];
	else if (!N_missing && S_missing)
	  d = - g[index1] - g[index1 - 1];
	else if (!N_missing && !S_missing)
	  d = 0.5 * (g[index1 + 1] - g[index1 - 1]);
	
	if (W_missing && !E_missing)
	  d += g[index2 + M] + g[index2];
	else if (!W_missing && E_missing)
	  d += - g[index2] - g[index2 - M];
	else if (!W_missing && !E_missing)
	  d += 0.5 * (g[index2 + M] - g[index2 - M]);
	
	if (U_missing && !D_missing)
	  d += g[index3 + MN] + g[index3];
	else if (!U_missing && D_missing)
	  d += - g[index3] - g[index3 - MN];
	else if (!U_missing && !D_missing)
	  d += 0.5 * (g[index3 + MN] - g[index3 - MN]);
	
	rhs[index1] = d;

	lhs[27 * index1 + 13] = -2 * ((!N_missing || !S_missing)
				      + (!W_missing || !E_missing)
				      + (!U_missing || !D_missing));
	if (!N_missing)
	  lhs[27 * index1 + 12] = 1 + (S_missing);
	if (!S_missing)
	  lhs[27 * index1 + 14] = 1 + (N_missing);
	if (!W_missing)
	  lhs[27 * index1 + 10] = 1 + (E_missing);
	if (!E_missing)
	  lhs[27 * index1 + 16] = 1 + (W_missing);
	if (!U_missing)
	  lhs[27 * index1 + 4]  = 1 + (D_missing);
	if (!D_missing)
	  lhs[27 * index1 + 22] = 1 + (U_missing);
      }

  /* Solve the equation system with the full multigrid algorithm.
   * Use W cycles and 2 presmoothing and 2 postsmoothing
   * Gauss-Seidel iterations.
   */
  poisson_full_multigrid3D(0, rhs, weight,
			   number_of_iterations, M, N, P, f_out);
  
  /* Fix the mean value. */
  sum = 0.0;
  num_samples_in_mask = 0;
  for (i = 0; i < M * N * P; i++)
    if (weight[i])
    {
      sum += f_out[i];
      num_samples_in_mask++;
    }
  
  mean = sum / num_samples_in_mask;
  for (i = 0; i < M * N * P; i++)
    if (weight[i])
    {
      f_out[i] -= mean;
      f_out[i] += mu;
    }

  mxFree(rhs);
  mxFree(weight);
  for (i = 0; i < MAX_LEVELS; i++)
    if (data.lhs[i] != NULL)
      mxFree(data.lhs[i]);
  mxDestroyArray(data.A_array);
}


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int M, N, P;
  double *g;
  double *mask = NULL;
  double *f_out;
  double mu;
  int number_of_iterations;
  int dim;
  int argno;
  
  /* Check the input and output arguments. */
  
  /* First we expect an MxNx2 or MxNxPx3 array. */
  dim = mxGetNumberOfDimensions(prhs[0]) - 1;
  if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0])
      || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0])
      || dim < 2 || dim > 3
      || mxGetDimensions(prhs[0])[dim] != dim)
  {
    mexErrMsgTxt("g is expected to be an MxNx2 or MxNxPx3 array.");
  }
  
  M = mxGetDimensions(prhs[0])[0];
  N = mxGetDimensions(prhs[0])[1];
  if (dim > 2)
    P = mxGetDimensions(prhs[0])[2];
  else
    P = 0;
  
  g = mxGetPr(prhs[0]);

  /* Then an optional mask. */
  argno = 1;
  if (argno >= nrhs)
    mask = NULL;
  else if (mxGetNumberOfElements(prhs[1]) > 1)
  {
    if (!mxIsNumeric(prhs [1]) || mxIsComplex(prhs [1])
	|| mxIsSparse(prhs [1]) || !mxIsDouble(prhs [1]))
    {
      mexErrMsgTxt("mask is expected to be a numeric array.");
    }

    if (mxGetNumberOfDimensions(prhs[1]) != dim
	|| mxGetDimensions(prhs[1])[0] != M
	|| mxGetDimensions(prhs[1])[1] != N
	|| (dim > 2 && mxGetDimensions(prhs[1])[2] != P))
    {
      mexErrMsgTxt("g and mask have incompatible sizes.");
    }
    
    mask = mxGetPr(prhs[1]);
    argno++;
  }
  else if (mxIsEmpty(prhs[1]))
  {
    mask = NULL;
    argno++;
  }

  /* Next two scalars. */
  if (argno >= nrhs)
    mu = 0.0;
  else
  {
    if (!mxIsNumeric(prhs[argno]) || mxIsComplex(prhs[argno])
	|| mxIsSparse(prhs[argno]) || !mxIsDouble(prhs[argno])
	|| mxGetNumberOfElements(prhs[argno]) != 1)
    {
      mexErrMsgTxt("mu is expected to be a scalar.");
    }
    mu = mxGetScalar(prhs[argno]);
  }
  argno++;

  if (argno >= nrhs)
  {
    number_of_iterations = 2;
  }
  else
  {
    if (!mxIsNumeric(prhs[argno]) || mxIsComplex(prhs[argno])
	|| mxIsSparse(prhs[argno]) || !mxIsDouble(prhs[argno])
	|| mxGetNumberOfElements(prhs[argno]) != 1)
    {
      mexErrMsgTxt("N is expected to be a scalar.");
    }
    number_of_iterations = (int) mxGetScalar(prhs[argno]);
    if (number_of_iterations < 0
	|| (double) number_of_iterations != mxGetScalar(prhs[argno]))
    {
      mexErrMsgTxt("N is expected to be a positive integer.");
    }
  }
  
  plhs[0] = mxCreateNumericArray(dim, mxGetDimensions(prhs[0]),
				 mxDOUBLE_CLASS, mxREAL);
  f_out = mxGetPr(plhs[0]);
  
  if (dim == 2)
    antigradient2D(g, mask, mu, number_of_iterations, M, N, f_out);
  else
    antigradient3D(g, mask, mu, number_of_iterations, M, N, P, f_out);
}

