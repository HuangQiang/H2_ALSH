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
	n_pts_  = n;
	dim_    = d;
	m_      = m;
	U_      = U;
	data_   = data;
	norm_d_ = norm_d;
	
	// -------------------------------------------------------------------------
	//  init srp_lsh
	// -------------------------------------------------------------------------
	int sign_alsh_dim = d + m;
	lsh_ = new SRP_LSH(n, sign_alsh_dim, K);
	lsh_->display();

	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * n;
	float *norm = new float[n];
	M_ = MINREAL;
	for (int i = 0; i < n; ++i) {
		norm[i] = norm_d[i][0];
		if (norm[i] > M_) M_ = norm[i];
	}

	// -------------------------------------------------------------------------
	//  build hash tables for srp_lsh for new format of data
	// -------------------------------------------------------------------------
	g_memory += (SIZEBOOL * K + SIZEFLOAT * sign_alsh_dim);
	bool  *hash_code = new bool[K];
	float *sign_alsh_data = new float[sign_alsh_dim];
	float scale = U / M_;
	int   exponent = -1;

	for (int i = 0; i < n; ++i) {
		// construct new format of data by sign-alsh transformation
		norm[i] *= scale;
		for (int j = 0; j < sign_alsh_dim; ++j) {
			if (j < d) {
				sign_alsh_data[j] = data[i][j] * scale;
			}
			else {
				exponent = (int) pow(2.0f, j - d + 1);
				sign_alsh_data[j] = 0.5f - pow(norm[i], exponent);
			}
		}

		// calc hash key for this new format of data
		for (int j = 0; j < K; ++j) {
			hash_code[j] = lsh_->calc_hash_code(j, sign_alsh_data);
		}
		lsh_->compress_hash_code((const bool*) hash_code, lsh_->hash_key_[i]);
	}

	// -------------------------------------------------------------------------
	//  build hash tables for qalsh for new format of data
	// -------------------------------------------------------------------------
	delete[] norm; norm = NULL;
	delete[] hash_code; hash_code = NULL;
	delete[] sign_alsh_data; sign_alsh_data = NULL;

	g_memory -= SIZEFLOAT * n;
	g_memory -= (SIZEBOOL * K + SIZEFLOAT * sign_alsh_dim);
}

// -----------------------------------------------------------------------------
Sign_ALSH::~Sign_ALSH()				// destructor
{
	delete lsh_; lsh_ = NULL;
}

// -----------------------------------------------------------------------------
void Sign_ALSH::display()			// display parameters
{
	printf("Parameters of Sign_ALSH:\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    m = %d\n",   m_);
	printf("    U = %.2f\n", U_);
	printf("    M = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int Sign_ALSH::kmip(				// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k mip results
{
	// -------------------------------------------------------------------------
	//  construct Sign_ALSH query
	// -------------------------------------------------------------------------
	int   sign_alsh_dim = dim_ + m_;
	float normq = norm_q[0];
	float *sign_alsh_query = new float[sign_alsh_dim];

	for (int i = 0; i < sign_alsh_dim; ++i) {
		if (i < dim_) sign_alsh_query[i] = query[i] / normq;
		else sign_alsh_query[i] = 0.0f;
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-AMC search by SRP-LSH
	// -------------------------------------------------------------------------
	std::vector<int> cand;
	lsh_->kmc(top_k, (const float *) sign_alsh_query, cand);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by SRP-LSH
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
	delete[] sign_alsh_query; sign_alsh_query = NULL;

	return 0;
}