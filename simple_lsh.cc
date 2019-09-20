#include "headers.h"

// -----------------------------------------------------------------------------
Simple_LSH::Simple_LSH(				// constructor
	int   n,							// number of data
	int   d,							// dimension of data
	int   K,							// number of hash tables
	FILE  *fp,							// output file pointer
	const float** data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	n_pts_ = n;
	dim_   = d;
	K_     = K;
	data_  = data;

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;	

	// -------------------------------------------------------------------------
	//  display parameters
	// -------------------------------------------------------------------------
	printf("Parameters of Simple_LSH:\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    K = %d\n",   K_);
	printf("    M = %f\n\n", M_);
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	fprintf(fp, "n          = %d\n", n_pts_);
	fprintf(fp, "d          = %d\n", dim_);
	fprintf(fp, "K          = %d\n", K_);
	fprintf(fp, "M          = %f\n", M_);
	fprintf(fp, "index_time = %f Seconds\n\n", indexing_time);
}

// -----------------------------------------------------------------------------
void Simple_LSH::bulkload()			// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	vector<float> norm_sqr(n_pts_, 0.0f);
	float max_norm_sqr = MINREAL;
	for (int i = 0; i < n_pts_; ++i) {
		norm_sqr[i] = calc_inner_product(dim_, data_[i], data_[i]);
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
Simple_LSH::~Simple_LSH()			// destructor
{
	delete lsh_; lsh_ = NULL;
	for (int i = 0; i < n_pts_; ++i) {
		delete[] simple_lsh_data_[i]; simple_lsh_data_[i] = NULL;
	}
	delete[] simple_lsh_data_; simple_lsh_data_ = NULL;
}

// -----------------------------------------------------------------------------
int Simple_LSH::kmip(				// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct Simple_LSH query
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	float *simple_lsh_query = new float[dim_ + 1];

	for (int i = 0; i < dim_; ++i) {
		simple_lsh_query[i] = query[i] / norm_q;
	}
	simple_lsh_query[dim_] = 0.0f;

	// -------------------------------------------------------------------------
	//  conduct c-k-AMC search by SRP-LSH
	// -------------------------------------------------------------------------
	MaxK_List *mcs_list = new MaxK_List(top_k);
	lsh_->kmc(top_k, (const float *) simple_lsh_query, mcs_list);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by SRP-LSH
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int   id = mcs_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}
	delete[] simple_lsh_query; simple_lsh_query = NULL;
	delete mcs_list; mcs_list = NULL;

	return 0;
}