#ifndef __QALSH_COL_H
#define __QALSH_COL_H

class MinK_List;

// -----------------------------------------------------------------------------
//  Structure of Hash Table
// -----------------------------------------------------------------------------
struct HashValue {					// structure of HashValue
	int id_;							// data point id
	float proj_;						// projection
};

// -----------------------------------------------------------------------------
//  QALSH_Col: memory version of query-aware locality-sensitive hashing based 
//  on dynamic collision counting for approximate nearest neighbor search.
// -----------------------------------------------------------------------------
class QALSH_Col {
public:
	QALSH_Col(						// constructor
		float** points,					// data points
		int n,							// number of points
		int d,							// dimension of data space
		float ratio);					// approximation ratio

	~QALSH_Col();					// destructor

	int knn(						// knn search based on collision counting
		float R,						// limited search range
		float* query,					// query point
		int top_k,						// top-k
		MinK_List* list);				// top-k nn results (return)

private:
	float** points_;				// data points
	int n_pts_;						// number of points
	int dim_;						// dimension of space
	float appr_ratio_;				// approximation ratio

	float w_;						// bucket width
	float p1_;						// positive probability
	float p2_;						// negative probability

	float alpha_;					// collision threshold percentage
	float beta_;					// false positive percentage
	float delta_;					// error probability

	int m_;							// number of hash tables
	int l_;							// collision threshold
	float* a_array_;				// hash functions
	HashValue** tables_;			// hash tables

	// -------------------------------------------------------------------------
	void calc_paras();				// calculate parameters

	void display_paras();			// display parameters

	void gen_hash_func();			// generate hash functions

	void build_hashtables();		// build hash tables

	// -------------------------------------------------------------------------
	float calc_p(					// calculate probability
		float x);						// parameter

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int table_id,					// table id
		float* point);					// one point
	
	// -------------------------------------------------------------------------
	int binary_search_pos(			// binary search hash tables for position
		int table_id,					// hash table id
		float value);					// hash value
};

// -----------------------------------------------------------------------------
int HashValueQsortComp(				// qsort assistant function
	const void* obj1,					// 1st object
	const void* obj2);					// 2nd object



#endif
