#ifndef __AMIP_H
#define __AMIP_H

// -----------------------------------------------------------------------------
//  interface of this package
// -----------------------------------------------------------------------------
int h2_alsh(						// c-AMIP search via h2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int l2_alsh(						// c-AMIP search via l2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int l2_alsh2(						// c-AMIP search via l2_alsh2
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int xbox(							// c-AMIP search via xbox
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int sign_alsh(						// c-AMIP search via sign_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int simple_lsh(						// c-AMIP search via simple_lsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear_scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision recall curve of h2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int sign_alsh_precision_recall(		// precision recall curve of sign_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for ANN search
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision recall curve of simple_lsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	float nn_ratio,						// approximation ratio for ANN search
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder);			// output folder


#endif // __AMIP_H
