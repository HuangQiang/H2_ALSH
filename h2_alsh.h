#ifndef __H2_ALSH_H
#define __H2_ALSH_H

class QALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  Assistant Data Structure for H2-ALSH
// -----------------------------------------------------------------------------
struct Block {
	int   n_pts_;
	float M_;
	int   *index_;
	QALSH *lsh_;

	Block() { n_pts_ = 0; M_ = 0; index_ = NULL; lsh_ = NULL; }
	~Block() {
		if (index_ != NULL) { delete[] index_; index_ = NULL; }
		if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
	}
};

// -----------------------------------------------------------------------------
//  Asymmetric Locality-Sensitive Hashing based on Homocentric Hypersphere 
//  partition (H2_ALSH) is used to solve the problem of c-Approximate Maximum 
//  Inner Product (c-AMIP) search
// -----------------------------------------------------------------------------
class H2_ALSH {
public:
	H2_ALSH(						// constructor
		int   n,						// number of data objects
		int   d,						// dimension of data objects
		float nn_ratio,					// approximation ratio for NN
		float mip_ratio,				// approximation ratio for MIP
		const float **data, 			// input data
		const float **norm_d);			// l2-norm of data objects

	// -------------------------------------------------------------------------
	~H2_ALSH();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kmip(						// k-MIP search
		int   top_k,					// top-k value
		const float *query,				// input query
		const float *norm_q,			// l2-norm of query
		MaxK_List *list);				// top-k MIP results (return) 

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimension of data objects
	float nn_ratio_;				// approximation ratio for NN
	float mip_ratio_;				// approximation ratio for MIP
	const float **data_;			// original data objects
	const float **norm_d_;			// l2-norm of data objects
	
	float b_;						// compression ratio
	float M_;						// max norm of the data objects
	float **h2_alsh_data_;			// h2_alsh data
	int   num_blocks_;				// number of blocks
	std::vector<Block*> blocks_;	// blocks
	
	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading
};

#endif // __H2_ALSH_H
