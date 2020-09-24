#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "srp_lsh.h"

namespace mips {

// -----------------------------------------------------------------------------
//  Sign-LSH is used to solve the problem of c-Approximate Maximum Inner 
//  Product (c-AMIP) search
//
//  the idea was introduced by Anshumali Shrivastava and Ping Li in their paper 
//  "Improved asymmetric locality sensitive hashing (ALSH) for Maximum Inner 
//  Product Search (MIPS)", In Proceedings of the Thirty-First Conference on 
//  Uncertainty in Artificial Intelligence (UAI), pages 812–821, 2015.
// -----------------------------------------------------------------------------
class Sign_ALSH {
public:
	Sign_ALSH(						// constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		int   K,						// number of hash tables
		int   m,						// additional dimension of data
		float U,						// scale factor for data
		const float **data, 			// input data
		const float **norm_d);			// l2-norm of data objects

	// -------------------------------------------------------------------------	
	~Sign_ALSH();					// destrcutor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kmip(						// c-k-AMIP search
		int   top_k,					// top-k value
		const float *query,				// input query
		const float *norm_q,			// l2-norm of query
		MaxK_List *list);				// top-k mip results
	
	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		ret += lsh_->get_memory_usage();
		return ret;
	}

protected:
	int   n_pts_;					// number of data points
	int   dim_;						// dimensionality
	int   m_;						// additional dimension of data
	float U_;						// scale factor
	float M_;						// max norm of data objects
	const float **data_;			// original data objects
	const float **norm_d_;			// l2-norm of data objects
	SRP_LSH *lsh_;					// SRP_LSH
};

} // end namespace mips
