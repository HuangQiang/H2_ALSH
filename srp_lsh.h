#ifndef __SRP_LSH_H
#define __SRP_LSH_H

// -----------------------------------------------------------------------------
//  Sign-Random Projection LSH (SRP_LSH) is used to solve the problem of 
//  c-Approximate Maximum Cosine (c-AMC) search
// 
//  the idea was introduced by Moses S Charikar in his paper "Similarity 
//  estimation techniques from rounding algorithms", In Proceedings of the 
//  thiry-fourth annual ACM symposium on Theory of computing (STOC), pages 
//  380–388, 2002.
// -----------------------------------------------------------------------------
class SRP_LSH {
public:
	SRP_LSH(
		int n,							// cardinality of dataset
		int d,							// dimensionality of dataset
		int K,							// number of hash tables
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~SRP_LSH();						// destructor

	// -------------------------------------------------------------------------
	int kmc(						// c-k-AMC search
		int   top_k,					// top-k value
		const float *query,				// input query
		MaxK_List *list);				// top-k MC results  (return)

protected:
	int   n_pts_;					// cardinality of dataset
	int   dim_;						// dimensionality of dataset
	int   K_;						// number of hash tables
	const float **data_;			// data objects

	bool  **hash_code_;				// hash code of data objects
	float **proj_;					// random projection vectors

	// -------------------------------------------------------------------------
	void gen_random_vectors();		// generate random projection vectors

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading

	// -------------------------------------------------------------------------
	void get_proj_vector(			// get vector after random projection
		const float *data,				// input data
		bool *hash_code);				// hash code of input data (return)
};

#endif // __SRP_LSH_H
