#include "l2_alsh2.h"

namespace mips {

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
	: n_pts_(n), dim_(d), m_(m), U_(U), data_(data), norm_d_(norm_d), norm_q_(norm_q)
{
	// -------------------------------------------------------------------------
	//  indexing the new format of data using qalsh
	// -------------------------------------------------------------------------
	int l2_alsh2_dim = d + 2 * m;
	lsh_ = new QALSH(n, l2_alsh2_dim, nn_ratio); 
	lsh_->display();
	
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of 
	//  data objects and queries
	// -------------------------------------------------------------------------
	float *norm = new float[n];
	M_ = MINREAL;
	for (int i = 0; i < n; ++i) {
		norm[i] = norm_d[i][0];
		if (norm[i] > M_) M_ = norm[i];
	}

	for (int i = 0; i < qn; ++i) {
		if (norm_q[i][0] > M_) M_ = norm_q[i][0];
	}

	// -------------------------------------------------------------------------
	//  construct new format of data
	// -------------------------------------------------------------------------
	float *l2_alsh2_data = new float[l2_alsh2_dim];
	float scale = U / M_;
	int   exponent = -1;
	int   lsh_m = lsh_->m_;
	
	for (int i = 0; i < n; ++i) {
		// construct new format of data by l2_alsh2 transformation
		norm[i] *= scale;
		for (int j = 0; j < l2_alsh2_dim; ++j) {
			if (j < d) {
				l2_alsh2_data[j] = data_[i][j] * scale;
			}
			else if (j < d + m) {
				exponent = (int) pow(2.0f, j - d + 1);
				l2_alsh2_data[j] = pow(norm[i], exponent);
			}
			else {
				l2_alsh2_data[j] = 0.5f;
			}
		}
		// calc hash value for new format of data
		for (int j = 0; j < lsh_m; ++j) {
			lsh_->tables_[j][i].id_  = i;
			lsh_->tables_[j][i].key_ = lsh_->calc_hash_value(j, l2_alsh2_data);
		}
	}
	for (int i = 0; i < lsh_m; ++i) {
		qsort(lsh_->tables_[i], n, sizeof(Result), ResultComp);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm;
	delete[] l2_alsh2_data;
}

// -----------------------------------------------------------------------------
L2_ALSH2::~L2_ALSH2()				// destructor
{
	if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
}

// -----------------------------------------------------------------------------
void L2_ALSH2::display()			// display parameters
{
	printf("Parameters of L2_ALSH2:\n");
	printf("    n  = %d\n",   n_pts_);
	printf("    d  = %d\n",   dim_);
	printf("    m  = %d\n",   m_);
	printf("    c0 = %.1f\n", lsh_->ratio_);
	printf("    U  = %.2f\n", U_);
	printf("    M  = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int L2_ALSH2::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct L2_ALSH2 query
	// -------------------------------------------------------------------------
	int   l2_alsh2_dim = dim_ + 2 * m_;
	int   exponent = -1;
	float scale = U_ / M_;
	float normq = norm_q[0];
	float *l2_alsh2_query = new float[l2_alsh2_dim];

	normq *= scale;
	for (int i = 0; i < l2_alsh2_dim; ++i) {
		if (i < dim_) {
			l2_alsh2_query[i] = query[i] * scale;
		}
		else if (i < dim_ + m_) {
			l2_alsh2_query[i] = 0.5f;
		}
		else {
			exponent = (int) pow(2.0f, i - dim_ - m_ + 1);
			l2_alsh2_query[i] = pow(normq, exponent);
		}
	}

	// -------------------------------------------------------------------------
	//  conduct c-k-ANN search by qalsh
	// -------------------------------------------------------------------------
	std::vector<int> cand;
	lsh_->knn(top_k, MAXREAL, (const float *) l2_alsh2_query, cand);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by qalsh
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
	delete[] l2_alsh2_query;

	return 0;
}

} // end namespace mips
