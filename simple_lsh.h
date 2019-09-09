#ifndef __SIMPLE_LSH_H
#define __SIMPLE_LSH_H

class SRP_LSH;

// -----------------------------------------------------------------------------
//  Simple-LSH is used to solve the problem of c-Approximate Maximum Inner 
//  Product (c-AMIP) search.
//
//  the idea was introduced by Behnam Neyshabur and Nathan Srebro in their 
//  paper "On Symmetric and Asymmetric LSHs for Inner Product Search", In 
//  Proceedings of the 32nd International Conference on International Conference 
//  on Machine Learning (ICML), pages 1926–1934, 2015.
// -----------------------------------------------------------------------------
class Simple_LSH {
public:
	Simple_LSH(						// constructor
		int   n,						// number of data
		int   d,						// dimension of data
		int   K,						// number of hash tables
		float ratio,					// approximation ratio
		FILE  *fp,						// output file pointer
		const float** data);			// data objects

	~Simple_LSH();					// destructor

	// -------------------------------------------------------------------------
	int kmip(						// c-k-AMIP search
		int   top_k,					// top-k value
		const float* query,				// input query
		MaxK_List* list);				// top-k mip results

protected:
	int   n_pts_;					// number of data points
	int   dim_;						// dimension of data
	int   K_;						// number of hash tables
	float appr_ratio_;				// approximation ratio for AMC search
	const float **data_;			// data objects
	
	float M_;						// max norm of data objects
	int   simple_lsh_dim_;			// dimension of simple_lsh data
	float **simple_lsh_data_;		// simple_lsh data
	SRP_LSH *lsh_;					// SRP_LSH

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading
};

#endif // __SIMPLE_LSH_H
