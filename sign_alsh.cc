#include <algorithm>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "srp_lsh.h"
#include "sign_alsh.h"

// -----------------------------------------------------------------------------
Sign_ALSH::Sign_ALSH(				// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	int   K,							// number of hash tables
	int   m,							// additional dimension of data
	float U,							// scale factor for data
	const float **data, 				// input data
	const float **norm_d)				// l2-norm of data objects

{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_         = n;
	dim_           = d;
	K_             = K;
	m_             = m;
	U_             = U;
	data_          = data;
	sign_alsh_dim_ = d + m;

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
Sign_ALSH::~Sign_ALSH()				// destructor
{
	delete lsh_; lsh_ = NULL;
	for (int i = 0; i < n_pts_; ++i) {
		delete[] sign_alsh_data_[i]; sign_alsh_data_[i] = NULL;
	}
	delete[] sign_alsh_data_; sign_alsh_data_ = NULL;
}

// -----------------------------------------------------------------------------
void Sign_ALSH::bulkload()			// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	std::vector<float> norm(n_pts_, 0.0f);
	M_ = MINREAL;
	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, data_[i], data_[i]));
		if (norm[i] > M_) M_ = norm[i];
	}

	// -------------------------------------------------------------------------
	//  construct new format of data
	// -------------------------------------------------------------------------
	float scale = U_ / M_;
	int   exponent = -1;

	sign_alsh_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		sign_alsh_data_[i] = new float[sign_alsh_dim_];

		norm[i] *= scale;
		for (int j = 0; j < sign_alsh_dim_; ++j) {
			if (j < dim_) {
				sign_alsh_data_[i][j] = data_[i][j] * scale;
			}
			else {
				exponent = (int) pow(2.0f, j - dim_ + 1);
				sign_alsh_data_[i][j] = 0.5f - pow(norm[i], exponent);
			}
		}
	}

	// -------------------------------------------------------------------------
	//  indexing the new format of data using srp-lsh
	// -------------------------------------------------------------------------
	lsh_ = new SRP_LSH(n_pts_, sign_alsh_dim_, K_, (const float **) sign_alsh_data_);
}

// -----------------------------------------------------------------------------
void Sign_ALSH::display()			// display parameters
{
	printf("Parameters of Sign_ALSH:\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    K = %d\n",   K_);
	printf("    m = %d\n",   m_);
	printf("    U = %.2f\n", U_);
	printf("    M = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int Sign_ALSH::kmip(				// c-k-AMIP search
	int   top_k,						// top-k value
	const float* query,					// input query
	MaxK_List* list)					// top-k mip results
{
	// -------------------------------------------------------------------------
	//  construct Sign_ALSH query
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	float *sign_alsh_query = new float[sign_alsh_dim_]; // dim + m
	
	for (int i = 0; i < sign_alsh_dim_; ++i) {
		if (i < dim_) sign_alsh_query[i] = query[i] / norm_q;
		else sign_alsh_query[i] = 0.0f;
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-AMC search by SRP-LSH
	// -------------------------------------------------------------------------
	MaxK_List *mcs_list = new MaxK_List(top_k);
	lsh_->kmc(top_k, (const float *) sign_alsh_query,  mcs_list);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by SRP-LSH
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int   id = mcs_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}
	delete[] sign_alsh_query; sign_alsh_query = NULL;
	delete mcs_list; mcs_list = NULL;

	return 0;
}