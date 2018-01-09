#include "headers.h"

// -----------------------------------------------------------------------------
Simple_LSH::Simple_LSH()			// constructor
{
	n_pts_ = -1;
	dim_ = -1;
	K_ = -1;

	appr_ratio_ = -1.0f;

	simple_lsh_dim_ = -1;
	simple_lsh_data_ = NULL;
	data_ = NULL;
}

// -----------------------------------------------------------------------------
Simple_LSH::~Simple_LSH()			// destructor
{
	if (simple_lsh_data_ != NULL) {
		for (int i = 0; i < n_pts_; i++) {
			delete[] simple_lsh_data_[i]; simple_lsh_data_[i] = NULL;
		}
		delete[] simple_lsh_data_; simple_lsh_data_ = NULL;
	}
	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void Simple_LSH::init(				// init the parameters
	int n,								// number of data
	int d,								// dimension of data
	int K,								// number of hash tables
	float ratio,						// approximation ratio
	float** data)						// input data
{
	n_pts_ = n;
	dim_ = d;
	K_ = K;

	appr_ratio_ = ratio;

	simple_lsh_dim_ = d + 1;
	simple_lsh_data_ = NULL;
	data_ = data;

	pre_processing();
}

// -----------------------------------------------------------------------------
int Simple_LSH::pre_processing()	// pre-processing of data
{
	// -------------------------------------------------------------------------
	//  calculate the norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	float* norm = new float[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		norm[i] = 0.0f;
	}

	M_ = MINREAL;
	for (int i = 0; i < n_pts_; i++) {
		norm[i] = 0.0f;
		for (int j = 0; j < dim_; j++) {
			norm[i] += data_[i][j] * data_[i][j];
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > M_) M_ = norm[i];
	}

	// -------------------------------------------------------------------------
	//  construct new data and indexing
	// -------------------------------------------------------------------------
	float scale = 1.0f / M_;
	int exponent = -1;

	printf("Construct Simple_LSH Data: ");
	simple_lsh_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		simple_lsh_data_[i] = new float[simple_lsh_dim_];

		norm[i] = norm[i] * scale;
		for (int j = 0; j < simple_lsh_dim_; j++) {
			if (j < dim_) {			// construct new data
				simple_lsh_data_[i][j] = data_[i][j] * scale;
			}
			else {
				simple_lsh_data_[i][j] = sqrt(1.0f - norm[i] * norm[i]);
			}
		}
	}
	printf("finish!\n\n");

	// -------------------------------------------------------------------------
	//  indexing the new data using SRP-LSH
	// -------------------------------------------------------------------------
	if (indexing()) return 1;

	display_params();				// display parameters
	delete[] norm; norm = NULL;		// Release space

	return 0;
}

// -----------------------------------------------------------------------------
void Simple_LSH::display_params()	// display parameters
{
	printf("Parameters of Simple_LSH:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    K = %d\n", K_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
}

// -----------------------------------------------------------------------------
int Simple_LSH::indexing()			// indexing the new data
{
	lsh_ = new SRP_LSH();
	lsh_->init(K_, simple_lsh_dim_, n_pts_, simple_lsh_data_);

	return 0;
}

// -----------------------------------------------------------------------------
int Simple_LSH::kmip(				// top-k approximate mip search
	float* query,						// input query
	int top_k,							// top-k value
	MaxK_List* list)					// top-k mip results
{
	int num_of_verf = 0;			// num of verification (NN and MIP calc)

	// -------------------------------------------------------------------------
	//  Construct Simple_LSH query
	// -------------------------------------------------------------------------
	float norm_q = 0.0f;			// calc norm of query
	for (int i = 0; i < dim_; i++) {
		norm_q += query[i] * query[i];
	}
	norm_q = sqrt(norm_q);

	float* simple_lsh_query = new float[simple_lsh_dim_]; // dim + m
	for (int i = 0; i < simple_lsh_dim_; i++) {
		if (i < dim_) {
			simple_lsh_query[i] = query[i] / norm_q;
		}
		else {
			simple_lsh_query[i] = 0.0f;
		}
	}

	// -------------------------------------------------------------------------
	//  Perform kmc search via SRP-LSH
	// -------------------------------------------------------------------------
	int simple_lsh_top_k = top_k;
	MaxK_List *mcs_list = new MaxK_List(simple_lsh_top_k);
	lsh_->kmc(simple_lsh_query, top_k, mcs_list);

	// -------------------------------------------------------------------------
	//  Compute inner product for candidates returned by SRP-LSH
	// -------------------------------------------------------------------------
	for (int i = 0; i < simple_lsh_top_k; i++) {
		int id = mcs_list->ith_largest_id(i);
		float ip = calc_inner_product(data_[id], query, dim_);

		list->insert(ip, id + 1);
	}
	num_of_verf += simple_lsh_top_k;

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] simple_lsh_query; simple_lsh_query = NULL;
	delete mcs_list; mcs_list = NULL;

	return num_of_verf;
}