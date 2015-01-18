/* Please see velocity_segment.m for documentation.
 *
 *  Author: Gunnar Farnebäck
 *          Computer Vision Laboratory
 *          Linköping University, Sweden
 *          gf@isy.liu.se
 */
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <float.h>
#include <stdio.h> /* For DumpCandidateHeap */
#include <string.h>

#define MAX_PARAMS 8

#define DEFAULT_MINSIZE 500
#define DEFAULT_MAXSIZE 500
#define DEFAULT_KERNELDIST 4
#define DEFAULT_LAMBDA 0.1 /* Factor to make it harder for new regions. */
#define DEFAULT_COVERAGE 1.0 /* How large part of the image to segment. */
/* Whether to recompute the parameters for a region at the end. */
#define DEFAULT_PARAMETER_RECOMPUTATION 0
/* How many times the initial blocks are regrown. */
#define DEFAULT_NUMBER_OF_INITIAL_ITERATIONS 2
/* Which motion models to use. Default is affine only. */
#define DEFAULT_MOTION_MODELS 2
#define DEFAULT_VERBOSE 0

/* Motion models */
#define CONSTANT_MOTION 1
#define AFFINE_MOTION 2
#define EIGHT_PARAMETER_MOTION 4

struct region
{
    double param[MAX_PARAMS];
    int motion_model;
    int startx;
    int starty;
    int size;
};

struct candidate
{
    double param[MAX_PARAMS];
    int motion_model;
    double cost;
    int x;
    int y;
    int last_rebuild;
};



static int
number_of_parameters(int motion_model)
{
    int n;
    switch (motion_model)
    {
      case CONSTANT_MOTION:
	n = 2;
	break;
      case AFFINE_MOTION:
	n = 6;
	break;
      case EIGHT_PARAMETER_MOTION:
	n = 8;
	break;
    }
    return n;
}

static void
clear_quadratic_form(double *Q, int motion_model)
{
    int i, j;
    int n = number_of_parameters(motion_model);
    for (j = 0; j < n+1; j++)
	for (i = 0; i < n+1; i++)
	    Q[j * (MAX_PARAMS+1) + i] = 0.0;
}

static void
add_to_quadratic_form(double *Q, double *T, int x, int y, int xsize,
		      int ysize, int motion_model, int debug)
{
    int index = x + y * xsize;
    int N = xsize * ysize;
    double T11 = T[index];
    double T12 = T[index + N];
    double T13 = T[index + 2*N];
    double T22 = T[index + 4*N];
    double T23 = T[index + 5*N];
    double T33 = T[index + 8*N];
    double x2, xy, y2;
    double x3, x2y, xy2, y3;
    double x4, x3y, x2y2, xy3, y4;
    if (debug)
	mexPrintf("%f %f %f %f %f %f\n", T11, T12, T13, T22, T23, T33);
    switch (motion_model)
    {
      case CONSTANT_MOTION:
	Q[0 + 0 * (MAX_PARAMS+1)] += T11;
	Q[1 + 1 * (MAX_PARAMS+1)] += T22;
	Q[2 + 2 * (MAX_PARAMS+1)] += T33;

	Q[0 + 1 * (MAX_PARAMS+1)] += T12;
	Q[0 + 2 * (MAX_PARAMS+1)] += T13;
	Q[1 + 2 * (MAX_PARAMS+1)] += T23;
	      		           
	Q[1 + 0 * (MAX_PARAMS+1)] += T12;
	Q[2 + 0 * (MAX_PARAMS+1)] += T13;
	Q[2 + 1 * (MAX_PARAMS+1)] += T23;
	break;

      case AFFINE_MOTION:
	x2 = x * x;
	xy = x * y;
	y2 = y * y;
	Q[0 + 0 * (MAX_PARAMS+1)] += x2 * T11;
	Q[1 + 1 * (MAX_PARAMS+1)] += y2 * T11;
	Q[2 + 2 * (MAX_PARAMS+1)] += T11;
	Q[3 + 3 * (MAX_PARAMS+1)] += x2 * T22;
	Q[4 + 4 * (MAX_PARAMS+1)] += y2 * T22;
	Q[5 + 5 * (MAX_PARAMS+1)] += T22;
	Q[6 + 6 * (MAX_PARAMS+1)] += T33;
	      
	Q[0 + 1 * (MAX_PARAMS+1)] += xy * T11;
	Q[0 + 2 * (MAX_PARAMS+1)] += x * T11;
	Q[0 + 3 * (MAX_PARAMS+1)] += x2 * T12;
	Q[0 + 4 * (MAX_PARAMS+1)] += xy * T12;
	Q[0 + 5 * (MAX_PARAMS+1)] += x * T12;
	Q[0 + 6 * (MAX_PARAMS+1)] += x * T13;
	Q[1 + 2 * (MAX_PARAMS+1)] += y * T11;
	Q[1 + 3 * (MAX_PARAMS+1)] += xy * T12;
	Q[1 + 4 * (MAX_PARAMS+1)] += y2 * T12;
	Q[1 + 5 * (MAX_PARAMS+1)] += y * T12;
	Q[1 + 6 * (MAX_PARAMS+1)] += y * T13;
	Q[2 + 3 * (MAX_PARAMS+1)] += x * T12;
	Q[2 + 4 * (MAX_PARAMS+1)] += y * T12;
	Q[2 + 5 * (MAX_PARAMS+1)] += T12;
	Q[2 + 6 * (MAX_PARAMS+1)] += T13;
	Q[3 + 4 * (MAX_PARAMS+1)] += xy * T22;
	Q[3 + 5 * (MAX_PARAMS+1)] += x * T22;
	Q[3 + 6 * (MAX_PARAMS+1)] += x * T23;
	Q[4 + 5 * (MAX_PARAMS+1)] += y * T22;
	Q[4 + 6 * (MAX_PARAMS+1)] += y * T23;
	Q[5 + 6 * (MAX_PARAMS+1)] += T23;
	      
	Q[1 + 0 * (MAX_PARAMS+1)] += xy * T11;
	Q[2 + 0 * (MAX_PARAMS+1)] += x * T11;
	Q[3 + 0 * (MAX_PARAMS+1)] += x2 * T12;
	Q[4 + 0 * (MAX_PARAMS+1)] += xy * T12;
	Q[5 + 0 * (MAX_PARAMS+1)] += x * T12;
	Q[6 + 0 * (MAX_PARAMS+1)] += x * T13;
	Q[2 + 1 * (MAX_PARAMS+1)] += y * T11;
	Q[3 + 1 * (MAX_PARAMS+1)] += xy * T12;
	Q[4 + 1 * (MAX_PARAMS+1)] += y2 * T12;
	Q[5 + 1 * (MAX_PARAMS+1)] += y * T12;
	Q[6 + 1 * (MAX_PARAMS+1)] += y * T13;
	Q[3 + 2 * (MAX_PARAMS+1)] += x * T12;
	Q[4 + 2 * (MAX_PARAMS+1)] += y * T12;
	Q[5 + 2 * (MAX_PARAMS+1)] += T12;
	Q[6 + 2 * (MAX_PARAMS+1)] += T13;
	Q[4 + 3 * (MAX_PARAMS+1)] += xy * T22;
	Q[5 + 3 * (MAX_PARAMS+1)] += x * T22;
	Q[6 + 3 * (MAX_PARAMS+1)] += x * T23;
	Q[5 + 4 * (MAX_PARAMS+1)] += y * T22;
	Q[6 + 4 * (MAX_PARAMS+1)] += y * T23;
	Q[6 + 5 * (MAX_PARAMS+1)] += T23;
	break;

      case EIGHT_PARAMETER_MOTION:
	x2   = x * x;
	xy   = x * y;
	y2   = y * y;
	x3   = x * x2;
	x2y  = x * xy;
	xy2  = x * y2;
	y3   = y * y2;
	x4   = x * x3;
	x3y  = x * x2y;
	x2y2 = x * xy2;
	xy3  = x * y3;
	y4   = y * y3;
	Q[0 + 0 * (MAX_PARAMS+1)] += T11;
	Q[1 + 1 * (MAX_PARAMS+1)] += x2 * T11;
	Q[2 + 2 * (MAX_PARAMS+1)] += y2 * T11;
	Q[3 + 3 * (MAX_PARAMS+1)] += T22;
	Q[4 + 4 * (MAX_PARAMS+1)] += x2 * T22;
	Q[5 + 5 * (MAX_PARAMS+1)] += y2 * T22;
	Q[6 + 6 * (MAX_PARAMS+1)] += x4 * T11 + 2.0 * x3y * T12 + x2y2 * T22;
	Q[7 + 7 * (MAX_PARAMS+1)] += x2y2 * T11 + 2.0 * xy3 * T12 + y4 * T22;
	Q[8 + 8 * (MAX_PARAMS+1)] += T33;
	      
	Q[0 + 1 * (MAX_PARAMS+1)] += x * T11;
	Q[0 + 2 * (MAX_PARAMS+1)] += y * T11;
	Q[0 + 3 * (MAX_PARAMS+1)] += T12;
	Q[0 + 4 * (MAX_PARAMS+1)] += x * T12;
	Q[0 + 5 * (MAX_PARAMS+1)] += y * T12;
	Q[0 + 6 * (MAX_PARAMS+1)] += x2 * T11 + xy * T12;
	Q[0 + 7 * (MAX_PARAMS+1)] += xy * T11 + y2 * T12;
	Q[0 + 8 * (MAX_PARAMS+1)] += T13;
	Q[1 + 2 * (MAX_PARAMS+1)] += xy * T11;
	Q[1 + 3 * (MAX_PARAMS+1)] += x * T12;
	Q[1 + 4 * (MAX_PARAMS+1)] += x2 * T12;
	Q[1 + 5 * (MAX_PARAMS+1)] += xy * T12;
	Q[1 + 6 * (MAX_PARAMS+1)] += x3 * T11 + x2y * T12;
	Q[1 + 7 * (MAX_PARAMS+1)] += x2y * T11 + xy2 * T12;
	Q[1 + 8 * (MAX_PARAMS+1)] += x * T13;
	Q[2 + 3 * (MAX_PARAMS+1)] += y * T12;
	Q[2 + 4 * (MAX_PARAMS+1)] += xy * T12;
	Q[2 + 5 * (MAX_PARAMS+1)] += y2 * T12;
	Q[2 + 6 * (MAX_PARAMS+1)] += x2y * T11 + xy2 * T12;
	Q[2 + 7 * (MAX_PARAMS+1)] += xy2 * T11 + y3 * T12;
	Q[2 + 8 * (MAX_PARAMS+1)] += y * T13;
	Q[3 + 4 * (MAX_PARAMS+1)] += x * T22;
	Q[3 + 5 * (MAX_PARAMS+1)] += y * T22;
	Q[3 + 6 * (MAX_PARAMS+1)] += x2 * T12 + xy * T22;
	Q[3 + 7 * (MAX_PARAMS+1)] += xy * T12 + y2 * T22;
	Q[3 + 8 * (MAX_PARAMS+1)] += T23;
	Q[4 + 5 * (MAX_PARAMS+1)] += xy * T22;
	Q[4 + 6 * (MAX_PARAMS+1)] += x3 * T12 + x2y * T22;
	Q[4 + 7 * (MAX_PARAMS+1)] += x2y * T12 + xy2 * T22;
	Q[4 + 8 * (MAX_PARAMS+1)] += x * T23;
	Q[5 + 6 * (MAX_PARAMS+1)] += x2y * T12 + xy2 * T22;
	Q[5 + 7 * (MAX_PARAMS+1)] += xy2 * T12 + y3 * T22;
	Q[5 + 8 * (MAX_PARAMS+1)] += y * T23;
	Q[6 + 7 * (MAX_PARAMS+1)] += x3y * T11 + 2.0 * x2y2 * T12 + xy3 * T22;
	Q[6 + 8 * (MAX_PARAMS+1)] += x2 * T13 + xy * T23;
	Q[7 + 8 * (MAX_PARAMS+1)] += xy * T13 + y2 * T23;
	      
	Q[1 + 0 * (MAX_PARAMS+1)] += x * T11;
	Q[2 + 0 * (MAX_PARAMS+1)] += y * T11;
	Q[3 + 0 * (MAX_PARAMS+1)] += T12;
	Q[4 + 0 * (MAX_PARAMS+1)] += x * T12;
	Q[5 + 0 * (MAX_PARAMS+1)] += y * T12;
	Q[6 + 0 * (MAX_PARAMS+1)] += x2 * T11 + xy * T12;
	Q[7 + 0 * (MAX_PARAMS+1)] += xy * T11 + y2 * T12;
	Q[8 + 0 * (MAX_PARAMS+1)] += T13;
	Q[2 + 1 * (MAX_PARAMS+1)] += xy * T11;
	Q[3 + 1 * (MAX_PARAMS+1)] += x * T12;
	Q[4 + 1 * (MAX_PARAMS+1)] += x2 * T12;
	Q[5 + 1 * (MAX_PARAMS+1)] += xy * T12;
	Q[6 + 1 * (MAX_PARAMS+1)] += x3 * T11 + x2y * T12;
	Q[7 + 1 * (MAX_PARAMS+1)] += x2y * T11 + xy2 * T12;
	Q[8 + 1 * (MAX_PARAMS+1)] += x * T13;
	Q[3 + 2 * (MAX_PARAMS+1)] += y * T12;
	Q[4 + 2 * (MAX_PARAMS+1)] += xy * T12;
	Q[5 + 2 * (MAX_PARAMS+1)] += y2 * T12;
	Q[6 + 2 * (MAX_PARAMS+1)] += x2y * T11 + xy2 * T12;
	Q[7 + 2 * (MAX_PARAMS+1)] += xy2 * T11 + y3 * T12;
	Q[8 + 2 * (MAX_PARAMS+1)] += y * T13;
	Q[4 + 3 * (MAX_PARAMS+1)] += x * T22;
	Q[5 + 3 * (MAX_PARAMS+1)] += y * T22;
	Q[6 + 3 * (MAX_PARAMS+1)] += x2 * T12 + xy * T22;
	Q[7 + 3 * (MAX_PARAMS+1)] += xy * T12 + y2 * T22;
	Q[8 + 3 * (MAX_PARAMS+1)] += T23;
	Q[5 + 4 * (MAX_PARAMS+1)] += xy * T22;
	Q[6 + 4 * (MAX_PARAMS+1)] += x3 * T12 + x2y * T22;
	Q[7 + 4 * (MAX_PARAMS+1)] += x2y * T12 + xy2 * T22;
	Q[8 + 4 * (MAX_PARAMS+1)] += x * T23;
	Q[6 + 5 * (MAX_PARAMS+1)] += x2y * T12 + xy2 * T22;
	Q[7 + 5 * (MAX_PARAMS+1)] += xy2 * T12 + y3 * T22;
	Q[8 + 5 * (MAX_PARAMS+1)] += y * T23;
	Q[7 + 6 * (MAX_PARAMS+1)] += x3y * T11 + 2.0 * x2y2 * T12 + xy3 * T22;
	Q[8 + 6 * (MAX_PARAMS+1)] += x2 * T13 + xy * T23;
	Q[8 + 7 * (MAX_PARAMS+1)] += xy * T13 + y2 * T23;
	
	break;
      default:
	mexErrMsgTxt("Unknown motion model.");
	break;
    }
}
    
static void
compute_optimal_param(double *Q, double *param, int motion_model)
{
    /* Here a linear equation system needs to be solved. */
    /* Solved by call back to matlab. */
    int n = number_of_parameters(motion_model);
    mxArray *Qbarmat = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Qbar = mxGetPr(Qbarmat);
    mxArray *qmat = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *q = mxGetPr(qmat);
    mxArray *leftargs[1];
    mxArray *rightargs[2];
    int i, j;

    for (j = 0; j < n; j++)
    {
	for (i = 0; i < n; i++)
	    Qbar[i + j * n] = Q[i + j * (MAX_PARAMS+1)];
	q[j] = -Q[n + j * (MAX_PARAMS+1)];
    }

    rightargs[0] = Qbarmat;
    rightargs[1] = qmat;
    mexCallMATLAB(1, leftargs, 2, rightargs, "\\");

    for (i = 0; i < n; i++)
	param[i] = mxGetPr(leftargs[0])[i];

    mxDestroyArray(Qbarmat);
    mxDestroyArray(qmat);
    mxDestroyArray(leftargs[0]);
}

static double
trace_from_q(double *Q, int motion_model)
{
    switch (motion_model)
    {
      case CONSTANT_MOTION:
	return (Q[0 + 0 * (MAX_PARAMS+1)] + Q[1 + 1 * (MAX_PARAMS+1)] +
		Q[2 + 2 * (MAX_PARAMS+1)]);
	break;
      case AFFINE_MOTION:
	return (Q[2 + 2 * (MAX_PARAMS+1)] + Q[5 + 5 * (MAX_PARAMS+1)] +
		Q[6 + 6 * (MAX_PARAMS+1)]);
	break;
      case EIGHT_PARAMETER_MOTION:
	return (Q[0 + 0 * (MAX_PARAMS+1)] + Q[3 + 3 * (MAX_PARAMS+1)] +
		Q[8 + 8 * (MAX_PARAMS+1)]);
	break;
    }
    return 1.0; /* Never reached */
}

static double
compute_normalized_cost(double *Q, double *param, int motion_model)
{
    double cost = 0.0;
    int i, j;
    int n = number_of_parameters(motion_model);
    for (j = 0; j < n+1; j++)
	for (i = 0; i < n+1; i++)
	    cost += param[i] * param[j] * Q[i + j * (MAX_PARAMS+1)];
    /* Obs, obs, obs! */
    cost /= trace_from_q(Q, motion_model);
    return cost;
}

static void
compute_velocity(double *vx, double *vy, int x, int y, double *param,
		 int motion_model)
{
    switch (motion_model)
    {
      case CONSTANT_MOTION:
	*vx = param[0];
	*vy = param[1];
	break;
      case AFFINE_MOTION:
	*vx = x*param[0] + y*param[1] + param[2];
	*vy = x*param[3] + y*param[4] + param[5];
	break;
      case EIGHT_PARAMETER_MOTION:
	*vx = param[0] + x*param[1] + y*param[2] + x*x*param[6] + x*y*param[7];
	*vy = param[3] + x*param[4] + y*param[5] + x*y*param[6] + y*y*param[7];
	break;
    }	
}

static double
normalized_distance(double *T, int x, int y, int xsize, int ysize,
		    double *param, int motion_model)
{
    int index = x + y * xsize;
    int N = xsize * ysize;
    double T11 = T[index];
    double T12 = T[index + N];
    double T13 = T[index + 2*N];
    double T22 = T[index + 4*N];
    double T23 = T[index + 5*N];
    double T33 = T[index + 8*N];
    double trace = T11 + T22 + T33;
    double vx;
    double vy;
    double distance;
    double speedsquare;
    compute_velocity(&vx, &vy, x, y, param, motion_model);
    distance = (vx*T11*vx + 2.0*vx*T12*vy + 2.0*vx*T13 + vy*T22*vy +
		2.0*vy*T23 + T33);
    speedsquare = vx*vx + vy*vy + 1;

    if (trace <= 0.0)
	return 1.0 / 3.0; /* Same as for multiple of I. */
    else
	return distance / (trace * speedsquare);
}

static void
add_to_heap(int *heap, int *heapsize, int index, double *costs)
{
    int i = *heapsize;
    heap[*heapsize] = index;
    (*heapsize)++;
    /* Walk up */
    while (i > 0)
    {
	/* This can be improved on by not dropping the moving value
	/* until has reached its position.
	 */
	if (costs[heap[i]] < costs[heap[(i+1)/2-1]])
	{
	    int j = (i+1)/2-1;
	    int tmp = heap[j];
	    heap[j] = heap[i];
	    heap[i] = tmp;
	    i = j;
	}
	else
	    break;
    }
}

static int
get_top_of_heap(int *heap, int *heapsize, double *costs)
{
    int first = heap[0];
    int last;
    int i, j;
    (*heapsize)--;
    last = heap[*heapsize];
    i = 0;
    while (2*i+1 < *heapsize)
    {
	j = 2*i+1;
	if (j+1 < *heapsize && costs[heap[j+1]] < costs[heap[j]])
	    j++;
	if (costs[heap[j]] < costs[last])
	{
	    heap[i] = heap[j];
	    i = j;
	}
	else
	    break;
    }
    heap[i] = last;
    return first;
}

static void
order_candidate_heap(int *heap, int n, struct candidate *candidates)
{
    int i, j, k;
    int current;
    for (k = (n+1)/2-1; k >= 0; k--)
    {
	current = heap[k];
	i = k;
	while (2*i+1 < n)
	{
	    j = 2*i+1;
	    if (j+1 < n
		&& candidates[heap[j+1]].cost < candidates[heap[j]].cost)
		j++;
	    if (candidates[heap[j]].cost < candidates[current].cost)
	    {
		heap[i] = heap[j];
		i = j;
	    }
	    else
		break;
	}
	heap[i] = current;
    }
}

static void
dump_candidate_heap(int *heap, int n, struct candidate *candidates,
		    char *filename)
{
    FILE *fd = fopen(filename, "w");
    int i;
    for (i = 0; i < n; i++)
	fprintf(fd, "%5d %5d %g\n", i+1, heap[i], candidates[heap[i]].cost);
    fclose(fd);
}

static void
reorder_candidate_heap(int *heap, int n, struct candidate *candidates)
{
    int current = heap[0];
    int i, j;
    i = 0;
    while (2*i+1 < n)
    {
	j = 2*i+1;
	if (j+1 < n && candidates[heap[j+1]].cost < candidates[heap[j]].cost)
	    j++;
	if (candidates[heap[j]].cost < candidates[current].cost)
	{
	    heap[i] = heap[j];
	    i = j;
	}
	else
	    break;
    }
    heap[i] = current;
}

static int
get_first_candidate(int *heap, int *n, struct candidate* candidates)
{
    int first = heap[0];
    (*n)--;
    heap[0] = heap[*n];
    reorder_candidate_heap(heap, *n, candidates);
    return first;
}

static int
rebuild_candidate(struct candidate *candidate, int *map,
		  int minsize, int maxsize, int xsize, int ysize,
		  int *heap, double *costs, int *pixels, double *T,
		  double *Q, int debug)
{
    int size = 0;
    int heapsize = 0;
    int x = candidate->x;
    int y = candidate->y;
    int index = x + y*xsize;
    double largest_cost = 0.0;
    int i;

    if (map[index] != 0)
    {
	candidate->cost = 0.0; /* Should never be used in this case. */
	return 0;
    }
    costs[index] = normalized_distance(T, x, y, xsize, ysize, candidate->param,
				       candidate->motion_model);
    add_to_heap(heap, &heapsize, index, costs);
    while (size < maxsize && heapsize > 0)
    {
	index = get_top_of_heap(heap, &heapsize, costs);
	
	if (map[index] != 0)
	    continue;
	
	pixels[size] = index;
	size++;
	map[index] = -1;
	if (largest_cost < costs[index])
	    largest_cost = costs[index];
	x = index % xsize;
	y = index / xsize;
	if (x+1 < xsize) /* South */
	{
	    int newindex = index + 1;
	    if (map[newindex] == 0)
	    {
		costs[newindex] = normalized_distance(T, x+1, y, xsize, ysize,
						      candidate->param,
						      candidate->motion_model);
		add_to_heap(heap, &heapsize, newindex, costs);
	    }
	}
	if (x-1 >= 0) /* North */
	{
	    int newindex = index - 1;
	    if (map[newindex] == 0)
	    {
		costs[newindex] = normalized_distance(T, x-1, y, xsize, ysize,
						      candidate->param,
						      candidate->motion_model);
		add_to_heap(heap, &heapsize, newindex, costs);
	    }
	}
	if (y+1 < ysize) /* East */
	{
	    int newindex = index + xsize;
	    if (map[newindex] == 0)
	    {
		costs[newindex] = normalized_distance(T, x, y+1, xsize, ysize,
						      candidate->param,
						      candidate->motion_model);
		add_to_heap(heap, &heapsize, newindex, costs);
	    }
	}
	if (y-1 >= 0) /* West */
	{
	    int newindex = index - xsize;
	    if (map[newindex] == 0)
	    {
		costs[newindex] = normalized_distance(T, x, y-1, xsize, ysize,
						      candidate->param,
						      candidate->motion_model);
		add_to_heap(heap, &heapsize, newindex, costs);
	    }
	}
    }
    
    if (debug)
	mexPrintf("Startpoint: %d %d\n", candidate->x, candidate->y);
    
    candidate->cost = largest_cost;
    for (i = 0; i < size; i++)
    {
	map[pixels[i]] = 0;
	if (debug)
	    mexPrintf("%d\n", pixels[i]);
    }
    
    if (size >= minsize)
    {
	int motion_model = candidate->motion_model;
	clear_quadratic_form(Q, motion_model);
	for (i = 0; i < size; i++)
	    add_to_quadratic_form(Q, T, pixels[i] % xsize, pixels[i] / xsize,
				  xsize, ysize, motion_model, 0);
	compute_optimal_param(Q, candidate->param, motion_model);
	if (debug)
	{
	    for (i=0; i<number_of_parameters(motion_model); i++)
		mexPrintf("%f ", candidate->param[i]);
	    mexPrintf("\n");
	}
    }
    return size;
}

static int
get_cheapest_pixel(int *heap, int *heapsize, double *costs,
		   int *heap_position)
{
    int first = heap[0];
    int last;
    int i, j;
    (*heapsize)--;
    last = heap[*heapsize];
    i = 0;
    while (2*i+1 < *heapsize)
    {
	j = 2*i+1;
	if (j+1 < *heapsize && costs[heap[j+1]] < costs[heap[j]])
	    j++;
	if (costs[heap[j]] < costs[last])
	{
	    heap[i] = heap[j];
	    heap_position[heap[i]] = i;
	    i = j;
	}
	else
	    break;
    }
    heap[i] = last;
    heap_position[heap[i]] = i;
    return first;
}

static void
remove_from_pixel_heap(int *heap, int *heapsize, double *costs,
		       int *heap_position, int index)
{
    int last;
    int i, j;
    (*heapsize)--;
    last = heap[*heapsize];
    i = heap_position[index];
    while (2*i+1 < *heapsize)
    {
	j = 2*i+1;
	if (j+1 < *heapsize && costs[heap[j+1]] < costs[heap[j]])
	    j++;
	if (costs[heap[j]] < costs[last])
	{
	    heap[i] = heap[j];
	    heap_position[heap[i]] = i;
	    i = j;
	}
	else
	    break;
    }
    heap[i] = last;
    heap_position[heap[i]] = i;
}

static void
add_to_pixel_heap(int *heap, int *heapsize, double *costs,
		  int *heap_position, int index)
{
    int i = *heapsize;
    heap[*heapsize] = index;
    heap_position[index] = *heapsize;
    (*heapsize)++;
    /* Walk up */
    while (i > 0)
    {
	/* This can be improved on by not dropping the moving value */
	/* until has reached its position. */
	if (costs[heap[i]] < costs[heap[(i+1)/2-1]])
	{
	    int j = (i+1)/2-1;
	    int tmp = heap[j];
	    heap[j] = heap[i];
	    heap[i] = tmp;
	    heap_position[heap[i]] = i;
	    heap_position[heap[j]] = j;
	    i = j;
	}
	else
	    break;
    }
}

static void
ascend_in_pixel_heap(int *heap, int heapsize, double *costs,
		     int *heap_position, int index)
{
    int i = heap_position[index];
    /* Walk up */
    while (i > 0)
    {
	/* This can be improved on by not dropping the moving value */
	/* until has reached its position. */
	if (costs[heap[i]] < costs[heap[(i+1)/2-1]])
	{
	    int j = (i+1)/2-1;
	    int tmp = heap[j];
	    heap[j] = heap[i];
	    heap[i] = tmp;
	    heap_position[heap[i]] = i;
	    heap_position[heap[j]] = j;
	    i = j;
	}
	else
	    break;
    }
}

static void
maybe_add_boundary_pixel(double *T, int x, int y, int xsize, int ysize,
			 struct region *regions, int n,
			 int *pixels_heap, int *pixels_heap_size,
			 double *pixel_costs, int *pixel_heap_position,
			 int newindex, int *pixel_claimed_by)
{
    double cost = normalized_distance(T, x, y, xsize, ysize, regions[n].param,
				      regions[n].motion_model);
    int new_boundary_pixel = (pixel_claimed_by[newindex] == 0);
    if (new_boundary_pixel || cost < pixel_costs[newindex])
    {
	pixel_costs[newindex] = cost;
	pixel_claimed_by[newindex] = n;
    }
    else
	return;
    
    if (new_boundary_pixel)
	add_to_pixel_heap(pixels_heap, pixels_heap_size, pixel_costs,
			  pixel_heap_position, newindex);
    else
	ascend_in_pixel_heap(pixels_heap, *pixels_heap_size, pixel_costs,
			     pixel_heap_position, newindex);
}

static void
recompute_all_parameters(int *map, struct region *regions, int
			 number_of_regions, double *T, int xsize,
			 int ysize) 
{
    double *Qs;
    int x, y;
    int index;
    int region;
    Qs = mxCalloc((MAX_PARAMS+1) * (MAX_PARAMS+1) * number_of_regions,
		  sizeof(*Qs));
    for (y = 0; y < ysize; y++)
	for (x = 0; x < xsize; x++)
	{
	    index = x + y * xsize;
	    region = map[index];
	    if (region > 0)
		add_to_quadratic_form(&(Qs[((MAX_PARAMS+1) *
					    (MAX_PARAMS+1) *
					    (region-1))]),
				      T, x, y, xsize, ysize,
				      regions[region].motion_model, 0);
	}
    for (region = 1; region <= number_of_regions; region++)
    {
	compute_optimal_param(&(Qs[(MAX_PARAMS+1)*(MAX_PARAMS+1)*(region-1)]),
			      regions[region].param,
			      regions[region].motion_model);
    }
    mxFree(Qs);
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int minsize = DEFAULT_MINSIZE;
    int maxsize = DEFAULT_MAXSIZE;
    int kerneldist = DEFAULT_KERNELDIST;
    double lambda = DEFAULT_LAMBDA;
    double coverage = DEFAULT_COVERAGE;
    int parameter_recomputation = DEFAULT_PARAMETER_RECOMPUTATION;
    int number_of_initial_iterations = DEFAULT_NUMBER_OF_INITIAL_ITERATIONS;
    int motion_models = DEFAULT_MOTION_MODELS;
    int verbose = DEFAULT_VERBOSE;
    
    int number_of_candidates = 0;
    struct candidate *candidates;
    
    const int *dims;
    int xsize;
    int ysize;
    int i, j;
    int x, y;
    
    double Q[(MAX_PARAMS+1) * (MAX_PARAMS+1)];
    double *T;
    int *map;
    int *candidate_build_heap;
    double *candidate_build_costs;
    int *candidate_build_pixels;
    int *candidate_heap;
    int *pixels_heap;
    double *pixel_costs;
    int *pixel_heap_position;
    int *pixel_claimed_by;
    int pixels_heap_size;
    int number_of_regions;
    int number_of_occupied_pixels;
    struct region *regions;
    int number_of_models;
    
    /* Check the number of input and output arguments. */
    if (nrhs < 1)
	mexErrMsgTxt("Too few input arguments.");
    if (nrhs > 2)
	mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 4)
	mexErrMsgTxt("Too many output arguments.");

    /* Check the formats of the input arguments. */
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
	mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) ||
	mxGetNumberOfDimensions(prhs[0]) != 4)
    {
	mexErrMsgTxt("T must be a real and full four-dimensional array, stored as doubles.");
    }
    
    if (nrhs == 2)
	if (!mxIsStruct(prhs[1]) || mxGetM(prhs[1]) * mxGetN(prhs[1]) > 1)
	    mexErrMsgTxt("OPTIONS must be a scalar structure.");
    
    /* We won't deal with an empty tensor field. */
    if (mxIsEmpty(prhs[0]))
	mexErrMsgTxt("T must not be empty.");
    
    dims = mxGetDimensions(prhs[0]);
    if (dims[2] != 3 || dims[3] != 3)
	mexErrMsgTxt("T must be M x N x 3 x 3.");
    xsize = dims[0];
    ysize = dims[1];
    T = mxGetPr(prhs[0]);
    
    /* Check if there are any options set. */
    if (nrhs == 2)
    {
	for (i = 0; i < mxGetNumberOfFields(prhs[1]); i++)
	{
	    const mxArray *fieldvalue = mxGetFieldByNumber(prhs[1], 0, i);
	    const char *fieldname = mxGetFieldNameByNumber(prhs[1], i);

	    if (strcmp(fieldname, "lambda") == 0)
		lambda = mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "minsize") == 0)
		minsize = (int) mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "maxsize") == 0)
		maxsize = (int) mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "kerneldist") == 0)
		kerneldist = (int) mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "coverage") == 0)
		coverage = mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "parameter_recomputation") == 0)
		parameter_recomputation = (int) mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "number_of_initial_iterations") == 0)
		number_of_initial_iterations = (int) mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "motion_models") == 0)
		motion_models = (int) mxGetScalar(fieldvalue);
	    else if (strcmp(fieldname, "verbose") == 0)
		verbose = (int) mxGetScalar(fieldvalue);
	    else
	    {
		const char *prefix = "Unknown option: ";
		char *buf = mxCalloc(strlen(prefix) + strlen(fieldname) + 1,
				     sizeof(*buf));
		sprintf(buf, "%s%s", prefix, fieldname);
		mexWarnMsgTxt(buf);
		mxFree(buf);
	    }

	    /* We can't have minsize larger than maxsize. */
	    if (minsize > maxsize)
		minsize = maxsize;

	    /* Verify that the field value really is a scalar. */
	    if (!mxIsNumeric(fieldvalue) || mxIsComplex(fieldvalue) ||
		mxIsSparse(fieldvalue) || !mxIsDouble(fieldvalue) ||
		mxGetM(fieldvalue)*mxGetN(fieldvalue) > 1)
	    {
		mexErrMsgTxt("The value for the options must be real scalars, stored as a doubles and not in sparse arrays.");
	    }
	}
    }

    map = mxCalloc(xsize * ysize, sizeof(*map));
    candidate_build_heap = mxCalloc(xsize * ysize,
				    sizeof(*candidate_build_heap));
    candidate_build_costs = mxCalloc(xsize * ysize,
				     sizeof(*candidate_build_costs));
    candidate_build_pixels = mxCalloc(xsize * ysize,
				      sizeof(*candidate_build_pixels));
    
    /* Count the number of initial candidate regions. */
    number_of_candidates = 0;
    for (y = 7; y < ysize-7; y += kerneldist)
	for (x = 7; x < xsize-7; x += kerneldist)
	    number_of_candidates++;

    if (motion_models < 1 || motion_models > 7)
	motion_models = DEFAULT_MOTION_MODELS;
    number_of_models = (((motion_models&1) != 0) + ((motion_models&2) != 0)
			+ ((motion_models&4) != 0));
    number_of_candidates *= number_of_models;
    candidates = mxCalloc(number_of_candidates, sizeof(*candidates));

    if (verbose)
	mexPrintf("%d\n", number_of_candidates);

    i = 0;
    for (y = 7; y < ysize-7; y += kerneldist)
	for (x = 7; x < xsize-7; x += kerneldist)
	    for (j = 0; j < 3; j++)
		if (motion_models & (1<<j))
		{
		    candidates[i].x = x;
		    candidates[i].y = y;
		    candidates[i].motion_model = 1<<j;
		    i++;
		}
    
    for (i = 0; i < number_of_candidates; i++)
    {
	int xx = candidates[i].x;
	int yy = candidates[i].y;
	int motion_model = candidates[i].motion_model;
	clear_quadratic_form(Q, motion_model);

	for (y = yy-10; y <= yy+10; y++)
	    for (x = xx-10; x <= xx+10; x++)
		if (x >= 0 && x < xsize && y >= 0 && y < ysize)
		    add_to_quadratic_form(Q, T, x, y, xsize, ysize,
					  motion_model, 0);
	compute_optimal_param(Q, candidates[i].param, motion_model);

	/* This value is actually not used, unless the subsequent loop is */
	/* removed. */
	if (number_of_initial_iterations == 0)
	    candidates[i].cost = compute_normalized_cost(Q,
							 candidates[i].param,
							 motion_model);

	/* Improve the parameters a specified number of times */
	/* (default 2) by regrowing the region and recomputing the */
	/* parameters. */
	for (j = 0; j < number_of_initial_iterations; j++)
	{
	    (void)rebuild_candidate(&candidates[i], map, minsize, maxsize,
				    xsize, ysize, candidate_build_heap,
				    candidate_build_costs,
				    candidate_build_pixels, T, Q, 0);
	}
	candidates[i].last_rebuild = 0;
	if (i % 200 == 0)
	{
	    if (verbose)
	    {
		mexPrintf("%d\n", i);
#if 0
		for (j = 0; j < number_of_parameters(motion_model); j++)
		    mexPrintf("%f ", candidates[i].param[j]);
		mexPrintf("\n");
#endif
	    }
	}
    }
    
    candidate_heap = mxCalloc(number_of_candidates, sizeof(*candidate_heap));
    for (i = 0; i < number_of_candidates; i++)
    {
	candidate_heap[i] = i;
    }
    order_candidate_heap(candidate_heap, number_of_candidates, candidates);

#if 0
    dump_candidate_heap(candidate_heap, number_of_candidates, candidates,
			"/tmp/candidate_dump2");
#endif
    
    pixels_heap = mxCalloc(xsize * ysize, sizeof(*pixels_heap));
    pixel_costs = mxCalloc(xsize * ysize, sizeof(*pixel_costs));
    pixel_heap_position = mxCalloc(xsize * ysize,
				   sizeof(*pixel_heap_position));
    pixel_claimed_by = mxCalloc(xsize * ysize, sizeof(*pixel_claimed_by));
    regions = mxCalloc(number_of_candidates + 1, sizeof(*regions));
    
    pixels_heap_size = 0;
    number_of_regions = 0;
    number_of_occupied_pixels = 0;

    while (number_of_occupied_pixels < (int)(coverage * xsize * ysize))
    {
	double add_pixel_cost;
	double add_candidate_cost;
	if (number_of_occupied_pixels % 1000 == 0 && verbose)
	{
	    mexPrintf("Regionsizes at %d pixels:\n",
		      number_of_occupied_pixels);
	    for (i = 1; i <= number_of_regions; i++)
		mexPrintf("%2d %5d\n", i, regions[i].size);
	}
	if (pixels_heap_size > 0)
	    add_pixel_cost = pixel_costs[pixels_heap[0]];
	else
	    add_pixel_cost = DBL_MAX;
	if (number_of_candidates > 0)
	    add_candidate_cost = candidates[candidate_heap[0]].cost;
	else
	    add_candidate_cost = DBL_MAX;

	if (add_pixel_cost < lambda * add_candidate_cost)
	{ /* Add the least expensive pixel to an existing region. */
	    int index = get_cheapest_pixel(pixels_heap, &pixels_heap_size,
					   pixel_costs, pixel_heap_position);
	    int x = index % xsize;
	    int y = index / xsize;
	    int n = pixel_claimed_by[index];
	    map[index] = n;
	    number_of_occupied_pixels++;
	    regions[n].size++;
	    if (x+1 < xsize) /* South */
	    {
		int newindex = index + 1;
		if (map[newindex] == 0)
		    maybe_add_boundary_pixel(T, x+1, y, xsize, ysize, regions,
					     n, pixels_heap, &pixels_heap_size,
					     pixel_costs, pixel_heap_position,
					     newindex, pixel_claimed_by);
	    }
	    if (x-1 >= 0) /* North */
	    {
		int newindex = index - 1;
		if (map[newindex] == 0)
		    maybe_add_boundary_pixel(T, x-1, y, xsize, ysize, regions,
					     n, pixels_heap, &pixels_heap_size,
					     pixel_costs, pixel_heap_position,
					     newindex, pixel_claimed_by);

	    }
	    if (y+1 < ysize) /* East */
	    {
		int newindex = index + xsize;
		if (map[newindex] == 0)
		    maybe_add_boundary_pixel(T, x, y+1, xsize, ysize, regions,
					     n, pixels_heap, &pixels_heap_size,
					     pixel_costs, pixel_heap_position,
					     newindex, pixel_claimed_by);
	    }
	    if (y-1 >= 0) /* West */
	    {
		int newindex = index - xsize;
		if (map[newindex] == 0)
		    maybe_add_boundary_pixel(T, x, y-1, xsize, ysize, regions,
					     n, pixels_heap, &pixels_heap_size,
					     pixel_costs, pixel_heap_position,
					     newindex, pixel_claimed_by);
	    }
	}
	else
	{ /* Consider raising a candidate to region status. */
	    int n = candidate_heap[0];

	    if (candidates[n].last_rebuild < number_of_occupied_pixels)
	    { /* Too old info */
		int size = rebuild_candidate(&candidates[n], map, minsize,
					     maxsize, xsize, ysize,
					     candidate_build_heap,
					     candidate_build_costs, 
					     candidate_build_pixels, T, Q, 0);
		
		candidates[n].last_rebuild = number_of_occupied_pixels;
		if (size < minsize) /* Remove completely. */
		{
		    (void) get_first_candidate(candidate_heap,
					       &number_of_candidates,
					       candidates);
		}
		else
		    reorder_candidate_heap(candidate_heap,
					   number_of_candidates,
					   candidates);
	    }
	    else
	    { /* Ok, let's raise. */
		int size = rebuild_candidate(&candidates[n], map, minsize,
					     maxsize, xsize, ysize,
					     candidate_build_heap,
					     candidate_build_costs,
					     candidate_build_pixels, T, Q, 0);
		int m;
		number_of_regions++;
		m = number_of_regions;
		if (verbose)
		    mexPrintf("Region %d(%d) has emerged with %d pixels.\n",
			      m, n, size);
		regions[m].startx = candidates[n].x;
		regions[m].starty = candidates[n].y;
		regions[m].motion_model = candidates[n].motion_model;
		regions[m].size = 0;
		for (i=0; i<number_of_parameters(regions[m].motion_model); i++)
		    regions[m].param[i] = candidates[n].param[i];
		for (i = 0; i < size; i++)
		{
		    int index = candidate_build_pixels[i];
		    int x, y;
		    map[index] = m;
		    if (pixel_claimed_by[index] != 0)
			remove_from_pixel_heap(pixels_heap, &pixels_heap_size,
					       pixel_costs,
					       pixel_heap_position, index);
		    x = index % xsize;
		    y = index / xsize;
		    number_of_occupied_pixels++;
		    regions[m].size++;
		    if (x+1 < xsize) /* South */
		    {
			int newindex = index + 1;
			if (map[newindex] == 0)
			    maybe_add_boundary_pixel(T, x+1, y, xsize, ysize,
						     regions, m, pixels_heap,
						     &pixels_heap_size, 
						     pixel_costs,
						     pixel_heap_position,
						     newindex,
						     pixel_claimed_by);
		    }
		    if (x-1 >= 0) /* North */
		    {
			int newindex = index - 1;
			if (map[newindex] == 0)
			    maybe_add_boundary_pixel(T, x-1, y, xsize, ysize,
						     regions, m, pixels_heap,
						     &pixels_heap_size,
						     pixel_costs,
						     pixel_heap_position,
						     newindex,
						     pixel_claimed_by);
			
		    }
		    if (y+1 < ysize) /* East */
		    {
			int newindex = index + xsize;
			if (map[newindex] == 0)
			    maybe_add_boundary_pixel(T, x, y+1, xsize, ysize,
						     regions, m, pixels_heap,
						     &pixels_heap_size,
						     pixel_costs,
						     pixel_heap_position,
						     newindex,
						     pixel_claimed_by);
		    }
		    if (y-1 >= 0) /* West */
		    {
			int newindex = index - xsize;
			if (map[newindex] == 0)
			    maybe_add_boundary_pixel(T, x, y-1, xsize, ysize,
						     regions, m, pixels_heap,
						     &pixels_heap_size,
						     pixel_costs,
						     pixel_heap_position,
						     newindex,
						     pixel_claimed_by);
		    }
		}
	    }
	}
    }

    /* Recompute parameters if this is enabled. */
    if (parameter_recomputation)
	recompute_all_parameters(map, regions, number_of_regions, T,
				 xsize, ysize);
    
    /* Return the velocity field. */
    if (nlhs > 0)
    {
	int dims[3];
	double *param;
	double vx;
	double vy;
	double *velocities;
	int motion_model;
	int index;
	dims[0] = xsize;
	dims[1] = ysize;
	dims[2] = 2;
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	velocities = mxGetPr(plhs[0]);
	for (y = 0; y < ysize; y++)
	    for (x = 0; x < xsize; x++)
	    {
		index = x + y * xsize;
		if (map[index] != 0)
		{
		    param = regions[map[index]].param;
		    motion_model = regions[map[index]].motion_model;
		    compute_velocity(&vx, &vy, x, y, param, motion_model);
		}
		else
		{
		    vx = 0.0;
		    vy = 0.0;
		}
		velocities[index] = vx;
		velocities[index + xsize * ysize] = vy;
	    }
    }
    
    /* Return the region map. */
    if (nlhs > 1)
    {
	plhs[1] = mxCreateDoubleMatrix(xsize, ysize, mxREAL);
	for (i = 0; i < xsize * ysize; i++)
	    mxGetPr(plhs[1])[i] = (double)map[i];
    }

    /* Return the parameters. */
    if (nlhs > 2)
    {
	double *new_params;
	mxArray *param_matrix;
	int n;
	plhs[2] = mxCreateCellMatrix(number_of_regions, 1);
	for (i = 0; i < number_of_regions; i++)
	{
	    n = number_of_parameters(regions[i+1].motion_model);
	    param_matrix = mxCreateDoubleMatrix(n, 1, mxREAL);
	    new_params = mxGetPr(param_matrix);
	    for (j=0; j<n; j++)
		new_params[j] = regions[i+1].param[j];
	    mxSetCell(plhs[2], i, param_matrix);
	}
    }

    /* Return the startpoints for the candidate regions which were raised */
    /* to real regions. (I.e., one point inside each real region.) */
    if (nlhs > 3)
    {
	double *centers;
	plhs[3] = mxCreateDoubleMatrix(2, number_of_regions, mxREAL);
	centers = mxGetPr(plhs[3]);
	for (i = 0; i < number_of_regions; i++)
	{
	    centers[2*i] = regions[i+1].startx;
	    centers[2*i + 1] = regions[i+1].starty;
	}
    }
}
