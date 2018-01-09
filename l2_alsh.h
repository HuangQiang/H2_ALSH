#ifndef __L2_ALSH_H
#define __L2_ALSH_H

class QALSH_Col;
class MaxK_List;

// -----------------------------------------------------------------------------
//  L2_ALSH: L2_ALSH is an Asymmetric LSH scheme based on Euclidean distance 
//  which is used to solve the problem of approximate Maximum Inner Product 
//  (MIP) search.
// -----------------------------------------------------------------------------
class L2_ALSH {
public:
	L2_ALSH();						// constructor
	~L2_ALSH();						// destructor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int n,							// number of data
		int d,							// dimension of data
		int m,							// additional dimension of data
		float U,						// scale factor for data
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
	int   m_;						// additional dimension of data
	float U_;						// scale factor
	float appr_ratio_;				// approximation ratio
	float **data_;					// original data 

	float M_;						// max norm of data
	int   l2_alsh_dim_;				// dimension of l2_alsh data (dim_ + m_)
	float **l2_alsh_data_;			// l2_alsh data

	QALSH_Col* lsh_;				// qalsh

	// -------------------------------------------------------------------------
	int pre_processing();			// pre-processing of data

	// -------------------------------------------------------------------------
	void display_params();			// display parameters

	// -------------------------------------------------------------------------
	int indexing();					// indexing the new data
};

#endif
