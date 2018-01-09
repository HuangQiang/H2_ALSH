#ifndef __SIGN_ALSH_H
#define __SIGN_ALSH_H

class SRP_LSH;

// -------------------------------------------------------------------------
class Sign_ALSH {
public:
	Sign_ALSH();					// constructor
	~Sign_ALSH();					// destrcutor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int n,							// number of data
		int d,							// dimension of data
		int K,							// number of hash tables
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
	int   K_;						// number of hash tables
	int   m_;						// additional dimension of data
	float U_;						// scale factor
	float appr_ratio_;				// approximation ratio
	float **data_;					// original data objects

	float M_;						// max norm of original data objects
	int   sign_alsh_dim_;			// dimension of sign_alsh data
	float **sign_alsh_data_;		// sign_alsh data
	
	SRP_LSH *lsh_;					// SRP_LSH

	// -------------------------------------------------------------------------
	int pre_processing();			// pre-processing of data

	// -------------------------------------------------------------------------
	void display_params();			// display parameters

	// -------------------------------------------------------------------------
	int indexing();					// indexing the new data

};

#endif