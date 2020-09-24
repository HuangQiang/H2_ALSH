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

class QALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  XBox is used to solve the problem of c-Approximate Maximum Inner Product 
//  (c-AMIP) search.
//
//  the idea was introduced by Yoram Bachrach, Yehuda Finkelstein, Ran 
//  Gilad-Bachrach, Liran Katzir, Noam Koenigstein, Nir Nice, and Ulrich Paquet 
//  in their paper "Speeding up the xbox recommender system using a euclidean 
//  transformation for inner-product spaces", In Proceedings of the 8th ACM 
//  Conference on Recommender systems, pages 257–264, 2014.
//
//  notice that to make a fair comparison with H2-ALSH, we apply QALSH for 
//  ANN search after converting MIP search to NN search by XBox transformation.
// -----------------------------------------------------------------------------
class XBox {
public:
	XBox(							// default constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		float nn_ratio,					// approximation ratio for ANN search
		const float **data, 			// input data
		const float **norm_d);			// l2-norm of data objects

	// -------------------------------------------------------------------------
	~XBox();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kmip(						// c-k-AMIP search
		int   top_k,					// top-k value
		bool  used_new_transform,		// used new transformation
		const float *query,				// input query
		const float *norm_q,			// l2-norm of query
		MaxK_List *list);				// top-k MIP results

	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		ret += lsh_->get_memory_usage();
		return ret;
	}

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	float M_;						// max norm of data objects
	const float **data_;			// original data objects
	const float **norm_d_;			// l2-norm of data objects
	QALSH *lsh_;					// qalsh
};
