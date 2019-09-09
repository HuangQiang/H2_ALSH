#include "headers.h"

// -----------------------------------------------------------------------------
L2_ALSH::L2_ALSH(					// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	int   m,							// additional dimension of data
	float U,							// scale factor for data
	float ratio,						// approximation ratio
	FILE  *fp,							// output file pointer
	const float** data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	n_pts_       = n;
	dim_         = d;
	m_           = m;
	U_           = U;
	appr_ratio_  = ratio;
	data_        = data;
	l2_alsh_dim_ = d + m;

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
	printf("Parameters of L2_ALSH:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    m = %d\n", m_);
	printf("    U = %.2f\n", U_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	fprintf(fp, "n          = %d\n", n_pts_);
	fprintf(fp, "d          = %d\n", dim_);
	fprintf(fp, "m          = %d\n", m_);
	fprintf(fp, "U          = %.2f\n", U_);
	fprintf(fp, "c          = %.2f\n", appr_ratio_);
	fprintf(fp, "M          = %.2f\n", M_);
	fprintf(fp, "index_time = %f Seconds\n\n", indexing_time);
}

// -----------------------------------------------------------------------------
void L2_ALSH::bulkload()			// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm
	// -------------------------------------------------------------------------
	M_ = MINREAL;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, data_[i], data_[i]));
		if (norm[i] > M_) M_ = norm[i];
	}

	// -------------------------------------------------------------------------
	//  construct new format of data and indexing
	// -------------------------------------------------------------------------
	float scale = U_ / M_;
	int   exponent = -1;

	l2_alsh_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		l2_alsh_data_[i] = new float[l2_alsh_dim_];

		norm[i] = norm[i] * scale;
		for (int j = 0; j < l2_alsh_dim_; ++j) {
			if (j < dim_) {
				l2_alsh_data_[i][j] = data_[i][j] * scale;
			}
			else {
				exponent = (int) pow(2.0f, j - dim_ + 1);
				l2_alsh_data_[i][j] = pow(norm[i], exponent);
			}
		}
	}

	// -------------------------------------------------------------------------
	//  indexing the new format of data using qalsh
	// -------------------------------------------------------------------------
	lsh_ = new QALSH(n_pts_, l2_alsh_dim_, appr_ratio_, (const float **) l2_alsh_data_);
}

// -----------------------------------------------------------------------------
L2_ALSH::~L2_ALSH()					// destructor
{
	if (l2_alsh_data_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] l2_alsh_data_[i]; l2_alsh_data_[i] = NULL;
		}
		delete[] l2_alsh_data_; l2_alsh_data_ = NULL;
	}

	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int L2_ALSH::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct L2_ALSH query
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	float *l2_alsh_query = new float[l2_alsh_dim_];

	for (int i = 0; i < l2_alsh_dim_; ++i) {
		if (i < dim_) l2_alsh_query[i] = query[i] / norm_q;
		else l2_alsh_query[i] = 0.5f;
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-ANN search by qalsh
	// -------------------------------------------------------------------------
	MinK_List *nn_list = new MinK_List(top_k);
	lsh_->knn(top_k, MAXREAL, (const float *) l2_alsh_query, nn_list);

	// -------------------------------------------------------------------------
	//  compute inner product for candidates returned by qalsh
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int id = nn_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] l2_alsh_query; l2_alsh_query = NULL;
	delete nn_list; nn_list = NULL;

	return 0;
}

