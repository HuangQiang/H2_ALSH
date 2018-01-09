#ifndef __H2_ALSH_H
#define __H2_ALSH_H

#include <vector>
using namespace std;

class QALSH_Col;
class MaxK_List;

// -----------------------------------------------------------------------------
//  We divide the dataset into a set of blocks which are the rings 
//  with same center
// -----------------------------------------------------------------------------
struct Block {
	int n_pts_;
	float M_;
	int* index_;
	QALSH_Col* lsh_;

	Block() { n_pts_ = 0; M_ = 0; index_ = NULL; lsh_ = NULL; }
	~Block() {
		if (index_ != NULL) { delete[] index_; index_ = NULL; }
		if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
	}
};

// -----------------------------------------------------------------------------
//  H2_ALSH: H2_ALSH is used to solve the problem of approximate Maximum Inner 
//  Product (MIP) search.
// -----------------------------------------------------------------------------
class H2_ALSH {
public:
	H2_ALSH();						// constructor
	virtual ~H2_ALSH();				// destructor

	// -------------------------------------------------------------------------
	void init(						// init the parameters
		int n,							// number of data
		int d,							// dimension of data
		float nn_ratio,					// approximation ratio for NN
		float mip_ratio,				// approximation ratio for MIP
		float** data);					// input data

	// -------------------------------------------------------------------------
	int kmip(						// top-k approximate mip search
		float* query,					// input query
		int top_k,						// top-k value
		MaxK_List* list);				// top-k mip results

private:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimension of data objects
	float nn_ratio_;				// approximation ratio for NN
	float mip_ratio_;				// approximation ratio for MIP
	float **data_;					// original data objects

	int   N_;						// number of sufficient condition
	float cmpr_ratio_;				// compression ratio
	float M_;						// max norm of the data objects

	int   h2_alsh_dim_;				// dimension of h2_alsh data (dim_ + 1)
	float **h2_alsh_data_;			// h2_alsh data
	int   block_count_;				// number of blocks
	vector<Block*> blocks_;			// blocks
	
	// -------------------------------------------------------------------------
	int pre_processing();			// pre-processing of data

	// -------------------------------------------------------------------------
	void display_params();			// display parameters

	// -------------------------------------------------------------------------
	void display_block_params(		// display parameters for block
		int block_index,				// block id
		Block* block);					// input block

	// -------------------------------------------------------------------------
	int indexing(					// indexing the new data
		int start_index,
		int n_pts,
		int block_index,
		Block* block);
};


#endif
