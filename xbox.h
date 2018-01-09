#ifndef __XBOX_H
#define __XBOX_H

class QALSH_Col;
class MaxK_List;

// -----------------------------------------------------------------------------
//  XBox: XBox is used to solve the problem of approximate Maximum Inner 
//  Product (MIP) search.
// -----------------------------------------------------------------------------
class XBox {
public:
	XBox();							// constructor
	~XBox();						// destructor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int n,							// number of data objects
		int d,							// dimension of data objects
		float ratio,					// approximation ratio
		float** data);					// original data objects

	// -------------------------------------------------------------------------
	int kmip(						// top-k approximate mip search
		float* query,					// input query
		int top_k,						// top-k value
		bool used_new_transform,		// used new transformation
		MaxK_List* list);				// top-k mip results

private:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimension of data objects
	float appr_ratio_;				// approximation ratio

	int   xbox_dim_;				// dimension of xbox data (dim_ + 1)
	float **xbox_data_;				// xbox data
	float **data_;					// original data objects
	float M_;						// max norm of data objects
	
	QALSH_Col* lsh_;				// qalsh

	// -------------------------------------------------------------------------
	void display_params();			// display parameters

	// -------------------------------------------------------------------------
	int pre_processing();			// pre-processing

	// -------------------------------------------------------------------------
	int indexing();					// indexing new transformation of data
};


#endif
