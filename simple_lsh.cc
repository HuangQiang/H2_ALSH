#include <algorithm>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "srp_lsh.h"
#include "simple_lsh.h"


// -----------------------------------------------------------------------------
Simple_LSH::Simple_LSH(				// constructor
	int   n,							// number of data
	int   d,							// dimension of data
	int   K,							// number of hash tables
	const float **data, 				// input data
	const float **norm_d)				// l2-norm of data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_  = n;
	dim_    = d;
	K_      = K;
	data_   = data;
	norm_d_ = norm_d;

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
Simple_LSH::~Simple_LSH()			// destructor
{
	delete lsh_; lsh_ = NULL;
	for (int i = 0; i < n_pts_; ++i) {
		delete[] simple_lsh_data_[i]; simple_lsh_data_[i] = NULL;
	}
	delete[] simple_lsh_data_; simple_lsh_data_ = NULL;
}

// -----------------------------------------------------------------------------
void Simple_LSH::bulkload()			// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	std::vector<float> norm_sqr(n_pts_, 0.0f);
	float max_norm_sqr = MINREAL;
	for (int i = 0; i < n_pts_; ++i) {
		norm_sqr[i] = norm_d_[i][0] * norm_d_[i][0];
		if (norm_sqr[i] > max_norm_sqr) max_norm_sqr = norm_sqr[i];
	}
	M_ = sqrt(max_norm_sqr);

	// -------------------------------------------------------------------------
	//  construct new format of data
	// -------------------------------------------------------------------------
	simple_lsh_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		simple_lsh_data_[i] = new float[dim_ + 1];
		for (int j = 0; j < dim_; ++j) {
			simple_lsh_data_[i][j] = data_[i][j] / M_;
		}
		simple_lsh_data_[i][dim_] = sqrt(1.0f - norm_sqr[i] / max_norm_sqr);
	}

	// -------------------------------------------------------------------------
	//  indexing the new data using SRP-LSH
	// -------------------------------------------------------------------------
	lsh_ = new SRP_LSH(n_pts_, dim_+1, K_, (const float **) simple_lsh_data_);
}

// -----------------------------------------------------------------------------
void Simple_LSH::display() 			// display parameters
{
	printf("Parameters of Simple_LSH:\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    K = %d\n",   K_);
	printf("    M = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int Simple_LSH::kmip(				// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k MIP results (return) 
{
	float kip   = MINREAL;
	float normq = norm_q[0];

	// -------------------------------------------------------------------------
	//  construct Simple_LSH query
	// -------------------------------------------------------------------------
	float *simple_lsh_query = new float[dim_ + 1];
	for (int i = 0; i < dim_; ++i) {
		simple_lsh_query[i] = query[i] / normq;
	}
	simple_lsh_query[dim_] = 0.0f;

	// -------------------------------------------------------------------------
	//  conduct c-k-AMC search by SRP-LSH
	// -------------------------------------------------------------------------
	std::vector<int> cand;
	lsh_->kmc(top_k, (const float *) simple_lsh_query, cand);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by SRP-LSH
	// -------------------------------------------------------------------------
	int size = (int) cand.size();
	for (int i = 0; i < size; ++i) {
		int id = cand[i];
		if (norm_d_[id][0] * normq <= kip) break;
				
		float ip = calc_inner_product(dim_, kip, data_[id], norm_d_[id], 
			query, norm_q);
		kip = list->insert(ip, id + 1);
	}
	delete[] simple_lsh_query; simple_lsh_query = NULL;

	return 0;
}