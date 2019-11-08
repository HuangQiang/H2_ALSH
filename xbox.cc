#include <algorithm>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"
#include "xbox.h"

// -----------------------------------------------------------------------------
XBox::XBox(							// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float nn_ratio,						// approximation ratio for ANN search
	const float **data, 				// input data
	const float **norm_d)				// l2-norm of data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	n_pts_      = n;
	dim_        = d;
	nn_ratio_   = nn_ratio;
	data_       = data;
	norm_d_     = norm_d;

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
XBox::~XBox()						// destructor
{
	delete lsh_; lsh_ = NULL;
	for (int i = 0; i < n_pts_; ++i) {
		delete[] xbox_data_[i]; xbox_data_[i] = NULL;
	}
	delete[] xbox_data_; xbox_data_ = NULL;
}

// -----------------------------------------------------------------------------
void XBox::bulkload()				// bulkloading
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	std::vector<float> norm_sqr(n_pts_, 0.0f);
	float max_norm_sqr = MINREAL;
	for (int i = 0; i < n_pts_; ++i) {
		norm_sqr[i] = calc_inner_product(dim_, data_[i], data_[i]);
		if (norm_sqr[i] > max_norm_sqr) max_norm_sqr = norm_sqr[i];
	}
	M_ = sqrt(max_norm_sqr);

	// -------------------------------------------------------------------------
	//  construct new data and indexing
	// -------------------------------------------------------------------------
	xbox_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		xbox_data_[i] = new float[dim_ + 1];
		for (int j = 0; j < dim_; ++j) {
			xbox_data_[i][j] = data_[i][j];
		}
		xbox_data_[i][dim_] = sqrt(max_norm_sqr - norm_sqr[i]);
	}
	
	// -------------------------------------------------------------------------
	//  indexing the new format of data using qalsh
	// -------------------------------------------------------------------------
	lsh_ = new QALSH(n_pts_, dim_+1, nn_ratio_, (const float **) xbox_data_);
}

// -----------------------------------------------------------------------------
void XBox::display()				// display parameters
{
	printf("Parameters of XBox:\n");
	printf("    n  = %d\n",   n_pts_);
	printf("    d  = %d\n",   dim_);
	printf("    c0 = %.1f\n", nn_ratio_);
	printf("    M  = %f\n\n", M_);
}

// -----------------------------------------------------------------------------
int XBox::kmip(						// c-k-AMIP search
	int   top_k,						// top-k value
	bool  used_new_transform,			// used new transformation
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  construct XBox query
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	float lambda = used_new_transform ? M_ / norm_q : 1.0f;
	float *xbox_query = new float[dim_ + 1];

	for (int i = 0; i < dim_; ++i) {
		xbox_query[i] = lambda * query[i];
	}
	xbox_query[dim_] = 0.0f;

	// -------------------------------------------------------------------------
	//  conduct c-k-ANN search by qalsh
	// -------------------------------------------------------------------------
	MinK_List *nn_list = new MinK_List(top_k);
	lsh_->knn(top_k, MAXREAL, (const float *) xbox_query, nn_list);

	// -------------------------------------------------------------------------
	//  calc inner product for candidates returned by qalsh
	// -------------------------------------------------------------------------
	for (int i = 0; i < top_k; ++i) {
		int   id = nn_list->ith_id(i);
		float ip = calc_inner_product(dim_, data_[id], query);

		list->insert(ip, id + 1);
	}
	delete[] xbox_query; xbox_query = NULL;
	delete nn_list; nn_list = NULL;

	return 0;
}

