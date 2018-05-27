#include "headers.h"

// -----------------------------------------------------------------------------
SRP_LSH::SRP_LSH(					// constructor
	int n,								// cardinality of dataset
	int d,								// dimensionality of dataset
	int K,								// number of hash tables
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	K_     = K;
	data_  = data;

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	gen_random_vectors();
	bulkload();
}

// -----------------------------------------------------------------------------
SRP_LSH::~SRP_LSH()					// destructor
{
	if (proj_ != NULL) {
		for (int i = 0; i < K_; ++i) {
			delete[] proj_[i];	proj_[i] = NULL;
		}
		delete[] proj_;	proj_ = NULL;
	}

	if (hash_code_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] hash_code_[i];	hash_code_[i] = NULL;
		}
		delete[] hash_code_; hash_code_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::gen_random_vectors()	// generate random projection vectors
{
	proj_ = new float*[K_];
	for (int i = 0; i < K_; ++i) {
		proj_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			proj_[i][j] = gaussian(0.0f, 1.0f);
		}
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::bulkload()			// bulkloading
{
	hash_code_ = new bool*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		hash_code_[i] = new bool[K_];
		get_proj_vector(data_[i], hash_code_[i]);
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::get_proj_vector(		// get vector after random projection
	const float *data,					// input data 
	bool *hash_code)					// hash code of input data (return)
{
	for (int i = 0; i < K_; ++i) {
		float sum = calc_inner_product(dim_, proj_[i], data);

		if (sum >= 0) hash_code[i] = true;
		else hash_code[i] = false;
	}
}

// -----------------------------------------------------------------------------
int SRP_LSH::kmc(					// c-k-AMC search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MC results (return)
{
	bool *mc_query = new bool[K_];
	get_proj_vector(query, mc_query);

	for (int i = 0; i < n_pts_; ++i) {
		int match = 0;
		for (int j = 0; j < K_; ++j) {
			if (hash_code_[i][j] == mc_query[j]) ++match;
		}
		list->insert(match, i);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] mc_query; mc_query = NULL;

	return 0;
}