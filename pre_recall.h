#ifndef __PRE_RECALL_H
#define __PRE_RECALL_H


// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision-recall curve of h2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int sign_alsh_precision_recall(		// precision-recall curve of sign_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision-recall curve of simple_lsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path);				// output path


#endif // __PRE_RECALL_H
