#include "headers.h"


// -----------------------------------------------------------------------------
Sign_ALSH::Sign_ALSH()				// constructor
{
	n_pts_ = -1;
	dim_ = -1;
	K_ = -1;

	m_ = -1;
	U_ = -1.0f;
	appr_ratio_ = -1.0f;

	sign_alsh_dim_ = -1;
	sign_alsh_data_ = NULL;
	data_ = NULL;
}

// -----------------------------------------------------------------------------
Sign_ALSH::~Sign_ALSH()				// destructor
{
	if (sign_alsh_data_ != NULL) {
		for (int i = 0; i < n_pts_; i++) {
			delete[] sign_alsh_data_[i]; sign_alsh_data_[i] = NULL;
		}
		delete[] sign_alsh_data_; sign_alsh_data_ = NULL;
	}
	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void Sign_ALSH::init(				// init the parameters
	int n,								// number of data
	int d,								// dimension of data
	int K,								// number of hash tables
	int m,								// additional dimension of data
	float U,							// scale factor for data
	float ratio,						// approximation ratio
	float** data)						// input data
{
	n_pts_ = n;
	dim_ = d;
	K_ = K;

	m_ = m;
	U_ = U;
	appr_ratio_ = ratio;

	sign_alsh_dim_ = d + m;
	sign_alsh_data_ = NULL;
	data_ = data;

	pre_processing();
}

// -----------------------------------------------------------------------------
int Sign_ALSH::pre_processing()		// pre-processing of data
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
	float scale = U_ / M_;
	int exponent = -1;

	printf("Construct Sign_ALSH Data: ");
	sign_alsh_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		sign_alsh_data_[i] = new float[sign_alsh_dim_];

		norm[i] = norm[i] * scale;
		for (int j = 0; j < sign_alsh_dim_; j++) {
			if (j < dim_) {
				sign_alsh_data_[i][j] = data_[i][j] * scale;
			}
			else {
				exponent = (int)pow(2.0f, j - dim_ + 1);
				sign_alsh_data_[i][j] = 0.5f - pow(norm[i], exponent);
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
void Sign_ALSH::display_params()	// display parameters
{
	printf("Parameters of Sign_ALSH:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    K = %d\n", K_);
	printf("    m = %d\n", m_);
	printf("    U = %.2f\n", U_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
}

// -----------------------------------------------------------------------------
int Sign_ALSH::indexing()			// indexing the new data
{
	lsh_ = new SRP_LSH();
	lsh_->init(K_, sign_alsh_dim_, n_pts_, sign_alsh_data_);

	return 0;
}

// -----------------------------------------------------------------------------
int Sign_ALSH::kmip(				// top-k approximate mip search
	float* query,						// input query
	int top_k,							// top-k value
	MaxK_List* list)					// top-k mip results
{
	int num_of_verf = 0;			// num of verification (NN and MIP calc)

	// -------------------------------------------------------------------------
	//  Construct Sign_ALSH query
	// -------------------------------------------------------------------------
	float norm_q = 0.0f;			// calc norm of query
	for (int i = 0; i < dim_; i++) {
		norm_q += query[i] * query[i];
	}
	norm_q = sqrt(norm_q);

	float* sign_alsh_query = new float[sign_alsh_dim_]; // dim + m
	for (int i = 0; i < sign_alsh_dim_; i++) {
		if (i < dim_) {
			sign_alsh_query[i] = query[i] / norm_q;
		}
		else {
			sign_alsh_query[i] = 0.0f;
		}
	}

	// -------------------------------------------------------------------------
	//  Perform kmc search via SRP-LSH
	// -------------------------------------------------------------------------
	int alsh_top_k = top_k;
	MaxK_List *mcs_list = new MaxK_List(alsh_top_k);
	lsh_->kmc(sign_alsh_query, top_k, mcs_list);

	// -------------------------------------------------------------------------
	//  Compute inner product for candidates returned by SRP-LSH
	// -------------------------------------------------------------------------
	for (int i = 0; i < alsh_top_k; i++) {
		int id = mcs_list->ith_largest_id(i);
		float ip = calc_inner_product(data_[id], query, dim_);

		list->insert(ip, id + 1);
	}
	num_of_verf += alsh_top_k;

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] sign_alsh_query; sign_alsh_query = NULL;
	delete mcs_list; mcs_list = NULL;

	return num_of_verf;
}