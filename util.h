#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include <sys/time.h>
#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "def.h"
#include "pri_queue.h"

struct Result;
class  MaxK_List;

extern timeval g_start_time;		// global param: start time
extern timeval g_end_time;			// global param: end time

extern float   g_memory;			// global param: estimated memory usage
extern float   g_indextime;			// global param: indexing time

extern float   g_runtime;			// global param: running time
extern float   g_ratio;				// global param: overall ratio
extern float   g_recall;			// global param: recall

// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_txt_data(					// read data (text) from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	const char *fname,					// address of data set
	float **data,						// data objects (return)
	float **norm_d);					// l2-norm of data objects (return)

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	const char *fname,					// address of data
	float **data,						// data objects (return)
	float **norm_d);					// l2-norm of data objects (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int    qn,							// number of query objects
	const  char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

// -----------------------------------------------------------------------------
void write_pre_recall(				// write precision-recall curves to disk
	FILE  *fp,							// file pointer
	const float **pre,					// precision 
	const float **recall);				// recall

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *norm1,					// l2-norm of 1st point
	const float *p2,					// 2nd point
	const float *norm2);				// l2-norm of 2nd point

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc L2 square distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_ratio(					// calc overall ratio
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
int get_hits(						// get the number of hits between two ID list
	int   k,							// top-k value
	int   t,							// top-t value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
int norm_distribution(				// analyse norm distribution of data
	int   n,							// number of data objects
	int   d,							// dimensionality
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const char  *out_path);				// output path

// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth MIP results
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// dimensionality
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const char  *truth_set);			// address of truth set
