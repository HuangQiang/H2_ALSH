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
	data_   = data;
	norm_d_ = norm_d;

	// -------------------------------------------------------------------------
	//  init srp_lsh
	// -------------------------------------------------------------------------
	lsh_ = new SRP_LSH(n, d + 1, K);
	lsh_->display();

	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * n;
	float *norm = new float[n];
	float max_norm = MINREAL;
	for (int i = 0; i < n; ++i) {
		norm[i] = SQR(norm_d[i][0]);
		if (norm[i] > max_norm) max_norm = norm[i];
	}
	M_ = sqrt(max_norm);

	// -------------------------------------------------------------------------
	//  build hash tables for srp_lsh for new format of data
	// -------------------------------------------------------------------------
	g_memory += (SIZEBOOL * K + SIZEFLOAT * (d + 1));
	bool  *hash_code = new bool[K];
	float *simple_lsh_data = new float[d + 1];
	for (int i = 0; i < n; ++i) {
		// construct new format of data by simple-lsh transformation
		for (int j = 0; j < d; ++j) {
			simple_lsh_data[j] = data[i][j] / M_;
		}
		simple_lsh_data[d] = sqrt(1.0f - norm[i] / max_norm);

		// calc hash key for this new format of data
		for (int j = 0; j < K; ++j) {
			hash_code[j] = lsh_->calc_hash_code(j, simple_lsh_data);
		}
		lsh_->compress_hash_code((const bool*) hash_code, lsh_->hash_key_[i]);
	}

	// -------------------------------------------------------------------------
	//  build hash tables for qalsh for new format of data
	// -------------------------------------------------------------------------
	delete[] norm; norm = NULL;
	delete[] hash_code; hash_code = NULL;
	delete[] simple_lsh_data; simple_lsh_data = NULL;

	g_memory -= SIZEFLOAT * n;
	g_memory -= (SIZEBOOL * K + SIZEFLOAT * (d + 1));
}

// -----------------------------------------------------------------------------
Simple_LSH::~Simple_LSH()			// destructor
{
	delete lsh_; lsh_ = NULL;
}

// -----------------------------------------------------------------------------
void Simple_LSH::display() 			// display parameters
{
	printf("Parameters of Simple_LSH:\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    M = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int Simple_LSH::kmip(				// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct Simple_LSH query
	// -------------------------------------------------------------------------
	float normq = norm_q[0];
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
	float kip  = MINREAL;
	int   size = (int) cand.size();
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