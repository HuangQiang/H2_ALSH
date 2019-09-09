#include "headers.h"

// -----------------------------------------------------------------------------
XBox::XBox(							// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float ratio,						// approximation ratio
	FILE  *fp,							// output file pointer
	const float** data) 				// original data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	n_pts_      = n;
	dim_        = d;
	appr_ratio_ = ratio;
	data_       = data;
	xbox_dim_   = d + 1;

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
	printf("Parameters of XBox:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	fprintf(fp, "n          = %d\n", n_pts_);
	fprintf(fp, "d          = %d\n", dim_);
	fprintf(fp, "c          = %.2f\n", appr_ratio_);
	fprintf(fp, "M          = %.2f\n", M_);
	fprintf(fp, "index_time = %f Seconds\n\n", indexing_time);
}

// -----------------------------------------------------------------------------
void XBox::bulkload()				// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	float max_norm_sqr = MINREAL;
	vector<float> norm_sqr(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm_sqr[i] = calc_inner_product(dim_, data_[i], data_[i]);
		if (norm_sqr[i] > max_norm_sqr) max_norm_sqr = norm_sqr[i];
	}
	M_ = sqrt(max_norm_sqr);

	// -------------------------------------------------------------------------
	//  construct new data and indexing
	// -------------------------------------------------------------------------
	xbox_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		xbox_data_[i] = new float[xbox_dim_];

		for (int j = 0; j < xbox_dim_; ++j) {
			if (j < dim_) xbox_data_[i][j] = data_[i][j];
			else xbox_data_[i][j] = sqrt(max_norm_sqr - norm_sqr[i]);
		}
	}
	
	// -------------------------------------------------------------------------
	//  indexing the new format of data using qalsh
	// -------------------------------------------------------------------------
	lsh_ = new QALSH(n_pts_, xbox_dim_, appr_ratio_, (const float **) xbox_data_);
}

// -----------------------------------------------------------------------------
XBox::~XBox()						// destructor
{
	if (xbox_data_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] xbox_data_[i]; xbox_data_[i] = NULL;
		}
		delete[] xbox_data_; xbox_data_ = NULL;
	}

	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int XBox::kmip(						// c-k-AMIP search
	int   top_k,						// top-k value
	bool  used_new_transform,			// used new transformation
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct XBox query
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	float lambda = M_ / norm_q;
	float *xbox_query = new float[xbox_dim_];

	for (int i = 0; i < xbox_dim_; ++i) {
		if (i < dim_) {
			if (used_new_transform) xbox_query[i] = lambda * query[i];
			else xbox_query[i] = query[i];
		}
		else {
			xbox_query[i] = 0.0f;
		}
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-ANN search by qalsh
	// -------------------------------------------------------------------------
	MinK_List *nn_list = new MinK_List(top_k);
	lsh_->knn(top_k, MAXREAL, (const float *) xbox_query, nn_list);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by qalsh
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int id = nn_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] xbox_query; xbox_query = NULL;
	delete nn_list; nn_list = NULL;

	return 0;
}

