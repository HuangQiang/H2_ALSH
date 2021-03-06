#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"

namespace mips {

// -----------------------------------------------------------------------------
//  Assistant Data Structure for H2-ALSH
// -----------------------------------------------------------------------------
struct Block {
	int   n_pts_;
	float M_;
	int   *index_;
	QALSH *lsh_;

	Block() { n_pts_ = 0; M_ = 0; index_ = NULL; lsh_ = NULL; }
	~Block() { if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; } }
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

	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		ret += SIZEINT * n_pts_; 	// for h2_alsh_id_
		for (auto block : blocks_) { // for blocks_
			ret += sizeof(*block);
			if (block->lsh_ != NULL) ret += block->lsh_->get_memory_usage();
		}
		return ret;
	}

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimension of data objects
	float ratio_;					// approximation ratio for MIP
	float b_;						// compression ratio
	float M_;						// max norm of the data objects
	const float **data_;			// original data objects
	const float **norm_d_;			// l2-norm of data objects
	
	int *h2_alsh_id_;				// data id after h2_alsh transformation
	std::vector<Block*> blocks_;	// blocks
};

} // end namespace mips
