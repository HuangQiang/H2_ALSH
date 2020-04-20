#ifndef __SRP_LSH_H
#define __SRP_LSH_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"

class MaxK_List;

// -----------------------------------------------------------------------------
//  Sign-Random Projection LSH (SRP_LSH) is used to solve the problem of 
//  c-Approximate Maximum Cosine (c-AMC) search
// 
//  the idea was introduced by Moses S. Charikar in his paper "Similarity 
//  estimation techniques from rounding algorithms", In Proceedings of the 
//  thiry-fourth annual ACM symposium on Theory of computing (STOC), pages 
//  380–388, 2002.
// -----------------------------------------------------------------------------
class SRP_LSH {
public:
	SRP_LSH(						// constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		int   K);						// number of hash functions

	// -------------------------------------------------------------------------
	~SRP_LSH();						// destructor

	// -------------------------------------------------------------------------
	bool calc_hash_code(			// calc hash code after random projection
		int   id,						// projection vector id
		const float *data);				// input data

	// -------------------------------------------------------------------------
	void compress_hash_code(		// compress hash code with 64 bits
		const bool *hash_code,			// input hash code
		uint64_t* hash_key);			// hash key (return)

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kmc(						// c-k-AMC search
		int   top_k,					// top-k value
		const float *query,				// input query
		std::vector<int> &cand); 		// MCS candidates  (return)

	// -------------------------------------------------------------------------
	int      n_;					// number of data objects
	int      d_;					// dimensionality
	int      K_;					// number of hash functions
	int      m_;					// number of compressed uint64_t hash code
	float    **proj_;				// random projection vectors
	uint64_t **hash_key_;			// hash code of data objects
	uint32_t *table16_;				// table to record the number of "1" bits

protected:
	// -------------------------------------------------------------------------
	uint32_t bit_count(				// count the number of 1 bits of x
		uint32_t x); 					// input uint32_t value

	// -------------------------------------------------------------------------
	uint32_t table_lookup(			// table lookup the match value
		uint64_t x);					// input uint64_t value
};

#endif // __SRP_LSH_H
