#include "l2_alsh.h"

// -----------------------------------------------------------------------------
L2_ALSH::L2_ALSH(					// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	int   m,							// additional dimension of data
	float U,							// scale factor for data
	float nn_ratio,						// approximation ratio for ANN search
	const float **data, 				// input data
	const float **norm_d)				// l2-norm of data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_       = n;
	dim_         = d;
	m_           = m;
	U_           = U;
	data_        = data;
	norm_d_      = norm_d;

	// -------------------------------------------------------------------------
	//  init qalsh
	// -------------------------------------------------------------------------
	int l2_alsh_dim = d + m;
	lsh_ = new QALSH(n, l2_alsh_dim, nn_ratio);
	lsh_->display();
	
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * n;
	float *norm = new float[n];
	M_ = MINREAL;
	for (int i = 0; i < n; ++i) {
		norm[i] = norm_d_[i][0];
		if (norm[i] > M_) M_ = norm[i];
	}

	// -------------------------------------------------------------------------
	//  build hash tables for qalsh for new format of data
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * l2_alsh_dim;
	float *l2_alsh_data = new float[l2_alsh_dim];
	float scale    = U / M_;
	int   exponent = -1;
	int   lsh_m = lsh_->m_;

	for (int i = 0; i < n; ++i) {
		// construct new format of data by l2_alsh transformation
		norm[i] *= scale;
		for (int j = 0; j < l2_alsh_dim; ++j) {
			if (j < d) {
				l2_alsh_data[j] = data[i][j] * scale;
			}
			else {
				exponent = (int) pow(2.0f, j-d+1);
				l2_alsh_data[j] = pow(norm[i], exponent);
			}
		}
		// calc hash value for new format of data
		for (int j = 0; j < lsh_m; ++j) {
			lsh_->tables_[j][i].id_  = i;
			lsh_->tables_[j][i].key_ = lsh_->calc_hash_value(j, l2_alsh_data);
		}
	}
	for (int i = 0; i < lsh_m; ++i) {
		qsort(lsh_->tables_[i], n, sizeof(Result), ResultComp);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm; norm = NULL;
	delete[] l2_alsh_data; l2_alsh_data = NULL;

	g_memory -= SIZEFLOAT * n;
	g_memory -= SIZEFLOAT * l2_alsh_dim;
}

// -----------------------------------------------------------------------------
L2_ALSH::~L2_ALSH()					// destructor
{
	delete lsh_; lsh_ = NULL;
}

// -----------------------------------------------------------------------------
void L2_ALSH::display()				// display parameters
{
	printf("Parameters of L2_ALSH:\n");
	printf("    n  = %d\n",   n_pts_);
	printf("    d  = %d\n",   dim_);
	printf("    m  = %d\n",   m_);
	printf("    U  = %.2f\n", U_);
	printf("    M  = %f\n",   M_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int L2_ALSH::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct L2_ALSH query
	// -------------------------------------------------------------------------
	int   l2_alsh_dim = dim_ + m_;
	float normq = norm_q[0];
	float *l2_alsh_query = new float[l2_alsh_dim];

	for (int i = 0; i < l2_alsh_dim; ++i) {
		if (i < dim_) l2_alsh_query[i] = query[i] / normq;
		else l2_alsh_query[i] = 0.5f;
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-ANN search by qalsh
	// -------------------------------------------------------------------------
	std::vector<int> cand;
	lsh_->knn(top_k, MAXREAL, (const float *) l2_alsh_query, cand);

	// -------------------------------------------------------------------------
	//  compute inner product for candidates returned by qalsh
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
	delete[] l2_alsh_query; l2_alsh_query = NULL;

	return 0;
}

