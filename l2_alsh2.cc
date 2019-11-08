#include <algorithm>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"
#include "l2_alsh2.h"

// -----------------------------------------------------------------------------
L2_ALSH2::L2_ALSH2(					// constructor
	int   n,							// number of data objects
	int   qn,							// number of queries
	int   d,							// dimension of data objects
	int   m,							// additional dimension of data
	float U,							// scale factor for data
	float nn_ratio,						// approximation ratio for ANN search
	const float **data, 				// input data
	const float **norm_d,				// l2-norm of data objects
	const float **norm_q)				// queries
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_        = n;
	dim_          = d;
	m_            = m;
	U_            = U;
	nn_ratio_     = nn_ratio;
	data_         = data;
	norm_d_       = norm_d;
	norm_q_       = norm_q;
	l2_alsh2_dim_ = d + 2 * m;	

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload(qn, norm_q);
}

// -----------------------------------------------------------------------------
L2_ALSH2::~L2_ALSH2()				// destructor
{
	delete lsh_; lsh_ = NULL;
	for (int i = 0; i < n_pts_; ++i) {
		delete[] l2_alsh2_data_[i]; l2_alsh2_data_[i] = NULL;
	}
	delete[] l2_alsh2_data_; l2_alsh2_data_ = NULL;
}

// -----------------------------------------------------------------------------
void L2_ALSH2::bulkload(			// pre-processing of data
	int   qn,							// number of queries
	const float** norm_q)				// queries
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of 
	//  data objects and queries
	// -------------------------------------------------------------------------
	std::vector<float> norm(n_pts_, 0.0f);
	M_ = MINREAL;
	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = norm_d_[i][0];
		if (norm[i] > M_) M_ = norm[i];
	}

	for (int i = 0; i < qn; ++i) {
		if (norm_q[i][0] > M_) M_ = norm_q[i][0];
	}

	// -------------------------------------------------------------------------
	//  construct new format of data
	// -------------------------------------------------------------------------
	float scale = U_ / M_;
	int   exponent = -1;

	l2_alsh2_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		l2_alsh2_data_[i] = new float[l2_alsh2_dim_];

		norm[i] *= scale;
		for (int j = 0; j < l2_alsh2_dim_; ++j) {
			if (j < dim_) {
				l2_alsh2_data_[i][j] = data_[i][j] * scale;
			}
			else if (j < dim_ + m_) {
				exponent = (int) pow(2.0f, j - dim_ + 1);
				l2_alsh2_data_[i][j] = pow(norm[i], exponent);
			}
			else {
				l2_alsh2_data_[i][j] = 0.5f;
			}
		}
	}

	// -------------------------------------------------------------------------
	//  indexing the new format of data using qalsh
	// -------------------------------------------------------------------------
	lsh_ = new QALSH(n_pts_, l2_alsh2_dim_, nn_ratio_, (const float **) l2_alsh2_data_);
}

// -----------------------------------------------------------------------------
void L2_ALSH2::display()			// display parameters
{
	printf("Parameters of L2_ALSH2:\n");
	printf("    n  = %d\n",   n_pts_);
	printf("    d  = %d\n",   dim_);
	printf("    m  = %d\n",   m_);
	printf("    U  = %.2f\n", U_);
	printf("    c0 = %.1f\n", nn_ratio_);
	printf("    M  = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int L2_ALSH2::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct L2_ALSH2 query
	// -------------------------------------------------------------------------
	int   exponent = -1;
	float scale = U_ / M_;
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	float *l2_alsh2_query = new float[l2_alsh2_dim_];

	norm_q *= scale;
	for (int i = 0; i < l2_alsh2_dim_; ++i) {
		if (i < dim_) {
			l2_alsh2_query[i] = query[i] * scale;
		}
		else if (i < dim_ + m_) {
			l2_alsh2_query[i] = 0.5f;
		}
		else {
			exponent = (int) pow(2.0f, i - dim_ - m_ + 1);
			l2_alsh2_query[i] = pow(norm_q, exponent);
		}
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-ANN search by qalsh
	// -------------------------------------------------------------------------
	MinK_List *nn_list = new MinK_List(top_k);
	lsh_->knn(top_k, MAXREAL, l2_alsh2_query, nn_list);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by qalsh
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int   id = nn_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}
	delete[] l2_alsh2_query; l2_alsh2_query = NULL;
	delete nn_list; nn_list = NULL;

	return 0;
}
