#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <cmath>
#include <vector>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"

struct Result;

// -----------------------------------------------------------------------------
//  Query-Aware Locality-Sensitive Hashing (QALSH) is used to solve the problem 
//  of c-Approximate Nearest Neighbor (c-ANN) search.
//
//  the idea was introduced by Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong 
//  Fang, and Wilfred Ng in their paper "Query-aware locality-sensitive hashing 
//  for approximate nearest neighbor search", in Proceedings of the VLDB 
//  Endowment (PVLDB), 9(1), pages 1–12, 2015.
// -----------------------------------------------------------------------------
class QALSH {
public:
	int    n_;						// number of data objects
	int    d_;						// dimensionality
	float  ratio_;					// approximation ratio
	float  w_;						// bucket width
	int    m_;						// number of hash tables
	int    l_;						// collision threshold
	float  **a_;					// lsh functions
	Result **tables_;				// hash tables

	// -------------------------------------------------------------------------
	QALSH(							// constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		float ratio);					// approximation ratio

	// -------------------------------------------------------------------------
	~QALSH();						// destructor

	// -------------------------------------------------------------------------
	float calc_p(					// calc probability
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   tid,						// table id
		const float *data);				// input data

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int knn(						// c-k-ANN search
		int   top_k,					// top-k
		float R,						// limited search range
		const float *query,				// input query
		std::vector<int> &cand);		// NN candidates (return)

	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		ret += SIZEFLOAT * m_ * d_;	// for a_
		ret += sizeof(Result) * m_ * n_; // for tables_
		return ret;
	}
};
