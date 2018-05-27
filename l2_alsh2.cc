#include "headers.h"

// -----------------------------------------------------------------------------
L2_ALSH2::L2_ALSH2()				// default constructor
{
	n_pts_         = -1;
	dim_           = -1;
	m_             = -1;
	U_             = -1.0f;
	appr_ratio_    = -1.0f;
	data_          = NULL;

	M_             = -1.0f;
	l2_alsh2_dim_  = -1;
	l2_alsh2_data_ = NULL;
	lsh_           = NULL;
}

// -----------------------------------------------------------------------------
L2_ALSH2::~L2_ALSH2()				// destructor
{
	if (l2_alsh2_data_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] l2_alsh2_data_[i]; l2_alsh2_data_[i] = NULL;
		}
		delete[] l2_alsh2_data_; l2_alsh2_data_ = NULL;
	}

	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void L2_ALSH2::build(				// build index
	int   n,							// number of data objects
	int   qn,							// number of queries
	int   d,							// dimension of data objects
	int   m,							// additional dimension of data
	float U,							// scale factor for data
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const float **query)				// queries
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_        = n;
	dim_          = d;
	m_            = m;
	U_            = U;
	appr_ratio_   = ratio;
	data_         = data;
	l2_alsh2_dim_ = d + 2 * m;	

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload(qn, query);
	display();
}

// -----------------------------------------------------------------------------
int L2_ALSH2::bulkload(				// pre-processing of data
	int   qn,							// number of queries
	const float** query)				// queries
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of 
	//  data objects and queries
	// -------------------------------------------------------------------------
	M_ = MINREAL;
	vector<float> norm(n_pts_, 0.0f);
	
	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, data_[i], data_[i]));
		if (norm[i] > M_) M_ = norm[i];
	}

	float tmp_norm = -1.0f;
	for (int i = 0; i < qn; ++i) {
		tmp_norm = sqrt(calc_inner_product(dim_, query[i], query[i]));
		if (tmp_norm > M_) M_ = tmp_norm;
	}

	// -------------------------------------------------------------------------
	//  construct new format of data
	// -------------------------------------------------------------------------
	float scale = U_ / M_;
	int   exponent = -1;

	printf("Construct L2_ALSH2 Data\n\n");
	l2_alsh2_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		l2_alsh2_data_[i] = new float[l2_alsh2_dim_];

		norm[i] = norm[i] * scale;
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
	lsh_ = new QALSH(n_pts_, l2_alsh2_dim_, appr_ratio_, 
		(const float **) l2_alsh2_data_);

	return 0;
}

// -----------------------------------------------------------------------------
void L2_ALSH2::display()			// display parameters
{
	printf("Parameters of L2_ALSH2:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    m = %d\n", m_);
	printf("    U = %.2f\n", U_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n", M_);
	printf("\n");
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

	norm_q = norm_q * scale;
	for (int i = 0; i < l2_alsh2_dim_; ++i) {
		if (i < dim_) {
			l2_alsh2_query[i] = query[i] * scale;
		}
		else if (i < dim_ + m_) {
			l2_alsh2_query[i] = 0.5f;
		}
		else {
			exponent = (int)pow(2.0f, i - dim_ - m_ + 1);
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
		int id = nn_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] l2_alsh2_query; l2_alsh2_query = NULL;
	delete nn_list; nn_list = NULL;

	return 0;
}
