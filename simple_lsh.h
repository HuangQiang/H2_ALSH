#ifndef __SIMPLE_LSH_H
#define __SIMPLE_LSH_H

class SRP_LSH;


// -----------------------------------------------------------------------------
class Simple_LSH {
public:
	Simple_LSH();					// constructor
	~Simple_LSH();					// destructor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int n,							// number of data
		int d,							// dimension of data
		int K,							// number of hash tables
		float ratio,					// approximation ratio
		float** data);					// input data

	// -------------------------------------------------------------------------
	int kmip(						// top-k approximate mip search
		float* query,					// input query
		int top_k,						// top-k value
		MaxK_List* list);				// top-k mip results

private:
	int   n_pts_;					// number of data points
	int   dim_;						// dimension of data
	int   K_;						// number of hash tables
	float appr_ratio_;				// approximation ratio
	float **data_;					// original data objects
	
	float M_;						// max norm of original data objects
	int   simple_lsh_dim_;			// dimension of simple_lsh data
	float **simple_lsh_data_;		// simple_lsh data
	
	SRP_LSH *lsh_;					// SRP_LSH

	// -------------------------------------------------------------------------
	int pre_processing();			// pre-processing of data

	// -------------------------------------------------------------------------
	void display_params();			// display parameters

	// -------------------------------------------------------------------------
	int indexing();					// indexing the new data

};

#endif