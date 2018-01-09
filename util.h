#ifndef __UTIL_H
#define __UTIL_H

#include <map>
using std::map;
class MaxK_List;

// -----------------------------------------------------------------------------
// global functions
// -----------------------------------------------------------------------------
void error(							// an error message
	char* msg,							// an message
	bool is_exit);						// whether exit the program

// -----------------------------------------------------------------------------
int read_data(						// read data set from disk
	int n,								// number of data points
	int d,								// dimensionality
	char* fname,						// address of data
	float** data);						// data (return)

// -----------------------------------------------------------------------------
int check_path(						// check input path whether exist
	char* path);						// input path

// -----------------------------------------------------------------------------
float uniform(						// r.v. from Uniform(min, max)
	float min,							// min value
	float max);							// max value

// -----------------------------------------------------------------------------
float gaussian(						// r.v. from Gaussian(mean, sigma)
	float mean,							// mean value
	float sigma);						// std value

// -----------------------------------------------------------------------------
float normal_pdf(					// pdf of Guassian(mean, std)
	float x,							// variable
	float u,							// mean
	float sigma);						// standard error

// -----------------------------------------------------------------------------
float normal_cdf(					// cdf of N(0, 1) in (-inf, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

// -----------------------------------------------------------------------------
float new_cdf(						// cdf of N(0, 1) in [-x, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim);							// dimension

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc L2 square distance
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim);							// dimension

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L2 distance
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim);							// dimension

// -----------------------------------------------------------------------------
float calc_l1_dist(					// calc L1 distance
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim);							// dimension

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall between two ID list
	int n,								// length of ID list
	MaxK_List* list,					// ID list of k-nn method
	map<int, int>& id_rank);			// ground truth ID list

// -----------------------------------------------------------------------------
int get_hits(						// get the number of hits between two ID list
	int k,								// length of ID list
	int t,								// top-t
	MaxK_List* list,					// ID list of k-nn method
	map<int, int>& id_rank);			// ground truth ID list


#endif
