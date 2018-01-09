#ifndef __L2_ALSH2_H
#define __L2_ALSH2_H

class QALSH_Col;
class MaxK_List;

// -----------------------------------------------------------------------------
//  L2_ALSH2: L2_ALSH2 is an Asymmetric LSH scheme based on Euclidean distance 
//  which is used to solve the problem of approximate Maximum Inner Product 
//  (MIP) search.
// -----------------------------------------------------------------------------
class L2_ALSH2 {
public:
	L2_ALSH2();						// constructor
	~L2_ALSH2();					// destructor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int n,							// number of data
		int qn,							// number of queries
		int d,							// dimension of data
		int m,							// additional dimension of data
		float U,						// scale factor for data
		float ratio,					// approximation ratio
		float** data,					// input data
		float** query);					// input query

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

	float M_;						// max norm of data and query
	int   l2_alsh2_dim_;			// dim of l2_alsh2 data (dim_ + 2 * m_)
	float **l2_alsh2_data_;			// l2_alsh2 data

	QALSH_Col* lsh_;				// qalsh

	// -------------------------------------------------------------------------
	int pre_processing(				// pre-processing of data
		int qn,							// number of queries
		float** query);					// input query

	// -------------------------------------------------------------------------
	void display_params();			// display parameters

	// -------------------------------------------------------------------------
	int indexing();					// indexing the new data
};

#endif
