#ifndef __QALSH_H
#define __QALSH_H

class MinK_List;

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
	QALSH(							// constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		float ratio,					// approximation ratio
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~QALSH();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int knn(						// c-k-ANN search
		int   top_k,					// top-k
		float R,						// limited search range
		const float *query,				// input query
		MinK_List *list);				// top-k NN results (return)

protected:
	int    n_pts_;					// number of data objects
	int    dim_;					// dimensionality
	float  appr_ratio_;				// approximation ratio
	const  float **data_;			// data objects

	float  w_;						// bucket width
	float  p1_;						// positive probability
	float  p2_;						// negative probability
	float  alpha_;					// collision threshold percentage
	float  beta_;					// false positive percentage
	float  delta_;					// error probability
	int    m_;						// number of hash tables
	int    l_;						// collision threshold
	float  **a_;					// lsh functions
	Result **tables_;				// hash tables

	int    *freq_;					// frequency		
	int    *lpos_;					// left  position of hash table
	int    *rpos_;					// right position of hash table
	bool   *checked_;				// whether checked
	bool   *bucket_flag_;			// bucket flag
	bool   *range_flag_;			// range flag
	float  *q_val_;					// hash value of query

	// -------------------------------------------------------------------------
	float calc_p(					// calc probability
		float x);						// x = w / (2.0 * r)
};

#endif // __QALSH_H
