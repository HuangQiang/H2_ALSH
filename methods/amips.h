#pragma once

#include <iostream>
#include <algorithm>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "h2_alsh.h"
#include "l2_alsh.h"
#include "l2_alsh2.h"
#include "xbox.h"
#include "sign_alsh.h"
#include "simple_lsh.h"

namespace mips {

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
	
// -----------------------------------------------------------------------------
int linear_scan(					// k-MIP search by linear scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

// -----------------------------------------------------------------------------
int l2_alsh(						// k-MIP search by l2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   m,							// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for ANN search
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

// -----------------------------------------------------------------------------
int l2_alsh2(						// k-MIP search by l2_alsh2
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   m,							// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for ANN search
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

// -----------------------------------------------------------------------------
int xbox(							// k-MIP search by xbox
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	const char *method_name1,			// name of method
	const char *method_name2,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

// -----------------------------------------------------------------------------
int sign_alsh(						// k-MIP search by sign_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

// -----------------------------------------------------------------------------
int simple_lsh(						// k-MIP search by simple_lsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

// -----------------------------------------------------------------------------
int h2_alsh(						// k-MIP search by h2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R);					// MIP ground truth results

} // end namespace mips
