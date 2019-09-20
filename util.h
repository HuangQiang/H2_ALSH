#ifndef __UTIL_H
#define __UTIL_H

class MaxK_List;

extern timeval g_start_time;		// global parameter: start time
extern timeval g_end_time;			// global parameter: end time

extern float   g_runtime;			// global parameter: running time
extern float   g_ratio;				// global parameter: overall ratio
extern float   g_recall;			// global parameter: recall

// -----------------------------------------------------------------------------
//  struct Result
// -----------------------------------------------------------------------------
struct Result {						// structure for furthest neighbor / hash value
	float key_;							// distance / random projection value
	int   id_;							// object id
};

// -----------------------------------------------------------------------------
inline int cmp(						// cmp func for lower_bound (ascending)
	Result a, 							// 1st element
	Result b)							// 2nd element
{
	return a.key_ < b.key_;
}

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
//  uitlity functions
// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_data(						// read data from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	const char *fname,					// address of data set
	float **data,						// data objects (return)
	float **norm_d);					// l2-norm of data objects (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int    qn,							// number of query objects
	const  char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

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
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	const Result *result);				// MIP results

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
void k_mip_search(					// k-MIP search
	int   n, 							// number of data objects
	int   qn,							// number of query objects
	int   d, 							// dimensionality
	int   k,							// top-k value
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	Result **result);					// k-MIP results (return)

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

#endif // __UTIL_H
