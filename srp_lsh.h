#ifndef __SRP_LSH_H
#define __SRP_LSH_H

// -----------------------------------------------------------------------------
//  Sign-Random Projection LSH (SRP_LSH): SRP_LSH is used to solve the problem 
//  of approximate Maximum Cosine Search (MCS).
// -----------------------------------------------------------------------------
class SRP_LSH {
public:
	SRP_LSH();						// constructor
	~SRP_LSH();						// destructor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int K,							// number of hash tables
		int d,							// dimensionality of dataset
		int n,							// cardinality of dataset
		float** data);					// input data objects

	// -------------------------------------------------------------------------
	int kmc(						// top-k approximate maximum cosine search
		float* query,					// input query
		int top_k,						// top-k value
		MaxK_List* list);				// top-k mip results

private:
	int K_;							// number of hash tables
	int dim_;						// dimensionality of dataset
	int n_pts_;						// cardinality of dataset

	bool  **hash_code_;				// hash code of data objects
	float **proj_;					// random projection vectors
	float **data_;					// data objects

	// -------------------------------------------------------------------------
	void pre_processing();			// pre-processing of data

	// -------------------------------------------------------------------------
	void generate();				// generate random projection vectors

	// -------------------------------------------------------------------------
	void get_projection_vector(		// get vector after random projection
		float* data,					// input data
		bool* hash_code);				// hash code of input data (return)
};

#endif
