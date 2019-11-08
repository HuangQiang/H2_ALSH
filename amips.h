#ifndef __AMIPS_H
#define __AMIPS_H


// -----------------------------------------------------------------------------
int linear_scan(					// k-MIP search by linear scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int l2_alsh(						// k-MIP search by l2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   m,							// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int l2_alsh2(						// k-MIP search by l2_alsh2
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   m,							// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int xbox(							// k-MIP search by xbox
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int sign_alsh(						// k-MIP search by sign_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int simple_lsh(						// k-MIP search by simple_lsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int h2_alsh(						// k-MIP search by h2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

#endif // __AMIPS_H
