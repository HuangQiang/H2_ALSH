#include "headers.h"

// -----------------------------------------------------------------------------
Simple_LSH::Simple_LSH(				// constructor
	int   n,							// number of data
	int   d,							// dimension of data
	int   K,							// number of hash tables
	float ratio,						// approximation ratio
	FILE  *fp,							// output file pointer
	const float** data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	n_pts_          = n;
	dim_            = d;
	K_              = K;
	appr_ratio_     = ratio;
	data_           = data;
	simple_lsh_dim_ = d + 1;

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
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    K = %d\n", K_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	fprintf(fp, "n          = %d\n", n_pts_);
	fprintf(fp, "d          = %d\n", dim_);
	fprintf(fp, "K          = %d\n", K_);
	fprintf(fp, "c          = %.2f\n", appr_ratio_);
	fprintf(fp, "M          = %.2f\n", M_);
	fprintf(fp, "index_time = %f Seconds\n\n", indexing_time);
}

// -----------------------------------------------------------------------------
void Simple_LSH::bulkload()			// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	M_ = MINREAL;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, data_[i], data_[i]));
		if (norm[i] > M_) M_ = norm[i];
	}

	// -------------------------------------------------------------------------
	//  construct new format of data
	// -------------------------------------------------------------------------
	float scale = 1.0f / M_;
	int exponent = -1;

	simple_lsh_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		simple_lsh_data_[i] = new float[simple_lsh_dim_];

		norm[i] = norm[i] * scale;
		for (int j = 0; j < simple_lsh_dim_; ++j) {
			if (j < dim_) {
				simple_lsh_data_[i][j] = data_[i][j] * scale;
			}
			else {
				simple_lsh_data_[i][j] = sqrt(1.0f - norm[i] * norm[i]);
			}
		}
	}

	// -------------------------------------------------------------------------
	//  indexing the new data using SRP-LSH
	// -------------------------------------------------------------------------
	lsh_ = new SRP_LSH(n_pts_, simple_lsh_dim_, K_, (const float **) simple_lsh_data_);
}

// -----------------------------------------------------------------------------
Simple_LSH::~Simple_LSH()			// destructor
{
	if (simple_lsh_data_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] simple_lsh_data_[i]; simple_lsh_data_[i] = NULL;
		}
		delete[] simple_lsh_data_; simple_lsh_data_ = NULL;
	}

	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
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
	float *simple_lsh_query = new float[simple_lsh_dim_];

	for (int i = 0; i < simple_lsh_dim_; ++i) {
		if (i < dim_) simple_lsh_query[i] = query[i] / norm_q;
		else simple_lsh_query[i] = 0.0f;
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-AMC search by SRP-LSH
	// -------------------------------------------------------------------------
	MaxK_List *mcs_list = new MaxK_List(top_k);
	lsh_->kmc(top_k, (const float *) simple_lsh_query, mcs_list);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by SRP-LSH
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int id = mcs_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] simple_lsh_query; simple_lsh_query = NULL;
	delete mcs_list; mcs_list = NULL;

	return 0;
}