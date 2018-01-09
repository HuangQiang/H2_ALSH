#ifndef __AMIP_H
#define __AMIP_H


// -----------------------------------------------------------------------------
// assistant functions of main function
// -----------------------------------------------------------------------------
int ground_truth(					// output the ground truth results
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float **data,						// data set
	float **query,						// query set
	char *truth_set);					// address of ground truth file

// -----------------------------------------------------------------------------
int l2_alsh(						// mip search via l2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int l2_alsh2(						// mip search via l2_alsh2
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int xbox(							// mip search via xbox
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int h2_alsh(						// mip search via h2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int sign_alsh(						// mip search via sign_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	int m,								// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int simple_lsh(						// mip search via simple_lsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear_scan
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimensionality
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int l2_alsh_precision_recall(		// precision recall curve of l2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int l2_alsh2_precision_recall(		// precision recall curve of l2_alsh2
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int xbox_precision_recall(			// precision recall curve of via xbox
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision recall curve of h2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int sign_alsh_precision_recall(		// precision recall curve of sign_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	int m,								// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision recall curve of simple_lsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder);				// output folder

#endif
