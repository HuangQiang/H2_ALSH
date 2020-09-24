#include "xbox.h"

// -----------------------------------------------------------------------------
XBox::XBox(							// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float nn_ratio,						// approximation ratio for ANN search
	const float **data, 				// input data
	const float **norm_d)				// l2-norm of data objects
	: n_pts_(n), dim_(d), data_(data), norm_d_(norm_d)
{
	// -------------------------------------------------------------------------
	//  init qalsh
	// -------------------------------------------------------------------------
	lsh_ = new QALSH(n, d + 1, nn_ratio);
	lsh_->display();

	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	float *norm = new float[n];
	float max_norm = MINREAL;
	for (int i = 0; i < n; ++i) {
		norm[i] = SQR(norm_d_[i][0]);
		if (norm[i] > max_norm) max_norm = norm[i];
	}
	M_ = sqrt(max_norm);

	// -------------------------------------------------------------------------
	//  build hash tables for qalsh for new format of data
	// -------------------------------------------------------------------------
	int   m = lsh_->m_;	
	float *xbox_data = new float[d + 1];
	for (int i = 0; i < n; ++i) {
		// construct new format of data by xbox transformation
		for (int j = 0; j < d; ++j) {
			xbox_data[j] = data[i][j];
		}
		xbox_data[d] = sqrt(max_norm - norm[i]);

		// calc hash value for new format of data
		for (int j = 0; j < m; ++j) {
			lsh_->tables_[j][i].id_  = i;
			lsh_->tables_[j][i].key_ = lsh_->calc_hash_value(j, xbox_data);
		}
	}
	for (int i = 0; i < m; ++i) {
		qsort(lsh_->tables_[i], n, sizeof(Result), ResultComp);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm;
	delete[] xbox_data;
}

// -----------------------------------------------------------------------------
XBox::~XBox()						// destructor
{
	if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
}

// -----------------------------------------------------------------------------
void XBox::display()				// display parameters
{
	printf("Parameters of XBox:\n");
	printf("    n  = %d\n",   n_pts_);
	printf("    d  = %d\n",   dim_);
	printf("    c0 = %.1f\n", lsh_->ratio_);
	printf("    M  = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int XBox::kmip(						// c-k-AMIP search
	int   top_k,						// top-k value
	bool  used_new_transform,			// used new transformation
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct XBox query
	// -------------------------------------------------------------------------
	float normq  = norm_q[0];
	float lambda = used_new_transform ? M_ / normq : 1.0f;

	float *xbox_query = new float[dim_ + 1];
	for (int i = 0; i < dim_; ++i) {
		xbox_query[i] = lambda * query[i];
	}
	xbox_query[dim_] = 0.0f;

	// -------------------------------------------------------------------------
	//  find candidates by qalsh
	// -------------------------------------------------------------------------
	std::vector<int> cand;
	lsh_->knn(top_k, MAXREAL, (const float *) xbox_query, cand);

	// -------------------------------------------------------------------------
	//  check candidates by calculating actual inner product value with query
	// -------------------------------------------------------------------------
	int   size = (int) cand.size();
	float kip  = MINREAL;	

	for (int i = 0; i < size; ++i) {
		int id = cand[i];
		if (norm_d_[id][0] * normq <= kip) break;
			
		float ip = calc_inner_product(dim_, kip, data_[id], norm_d_[id], 
			query, norm_q);
		kip = list->insert(ip, id + 1);
	}
	delete[] xbox_query;

	return 0;
}

