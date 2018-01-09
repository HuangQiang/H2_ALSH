#include "headers.h"

// -----------------------------------------------------------------------------
XBox::XBox()						// constructor
{
	n_pts_ = -1;
	dim_ = -1;
	appr_ratio_ = -1.0f;

	M_ = -1.0f;
	xbox_dim_ = -1;
	
	data_ = NULL;
	xbox_data_ = NULL;
	lsh_ = NULL;
}

// -----------------------------------------------------------------------------
XBox::~XBox()						// destructor
{
	if (xbox_data_ != NULL) {
		for (int i = 0; i < n_pts_; i++) {
			delete[] xbox_data_[i]; xbox_data_[i] = NULL;
		}
		delete[] xbox_data_; xbox_data_ = NULL;
	}

	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void XBox::init(					// init the parameters
	int n,								// number of data objects
	int d,								// dimension of data objects
	float ratio,						// approximation ratio
	float** data)						// original data objects
{
	n_pts_ = n;
	dim_ = d;
	appr_ratio_ = ratio;
	data_ = data;

	M_ = -1.0f;
	xbox_dim_ = d + 1;

	xbox_data_ = NULL;
	lsh_ = NULL;

	pre_processing();
}

// -----------------------------------------------------------------------------
int XBox::pre_processing()			// pre-processing of data
{
	// -------------------------------------------------------------------------
	//  calculate the norm of data and find the maximum norm of data
	// -------------------------------------------------------------------------
	float max_norm_sqr = MINREAL;
	float* norm_sqr = new float[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		norm_sqr[i] = 0.0f;
	}

	for (int i = 0; i < n_pts_; i++) {
		norm_sqr[i] = 0.0f;
		for (int j = 0; j < dim_; j++) {
			norm_sqr[i] += data_[i][j] * data_[i][j];
		}

		if (norm_sqr[i] > max_norm_sqr) max_norm_sqr = norm_sqr[i];
	}
	M_ = sqrt(max_norm_sqr);		// init <M_>

	// -------------------------------------------------------------------------
	//  construct new data and indexing
	// -------------------------------------------------------------------------
	printf("Construct XBox Data: ");
	xbox_data_ = new float*[n_pts_];// init <alsh_data_>
	for (int i = 0; i < n_pts_; i++) {
		xbox_data_[i] = new float[xbox_dim_];

		for (int j = 0; j < xbox_dim_; j++) {
			if (j < dim_) {			// construct new data
				xbox_data_[i][j] = data_[i][j];
			}
			else {
				xbox_data_[i][j] = sqrt(max_norm_sqr - norm_sqr[i]);
			}
		}
	}
	printf("finish!\n\n");

	// -------------------------------------------------------------------------
	//  indexing the new data using qalsh
	// -------------------------------------------------------------------------
	if (indexing()) return 1;

	display_params();
	delete[] norm_sqr; norm_sqr = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void XBox::display_params()			// display parameters
{
	printf("Parameters of XBox:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
}

// -----------------------------------------------------------------------------
int XBox::indexing()				// indexing the new data
{
	lsh_ = new QALSH_Col(xbox_data_, n_pts_, xbox_dim_, appr_ratio_);

	if (lsh_ != NULL) return 0;
	else return 1;
}

// -----------------------------------------------------------------------------
int XBox::kmip(						// top-k approximate mip search
	float* query,						// input query
	int top_k,							// top-k value
	bool used_new_transform,			// used new transformation
	MaxK_List* list)					// top-k mip results
{
	int num_of_verf = 0;			// num of verification (NN and MIP calc)

	// -------------------------------------------------------------------------
	//  Construct XBox query
	// -------------------------------------------------------------------------
	float norm_q = 0.0f;
	for (int i = 0; i < dim_; i++) {
		norm_q += query[i] * query[i];
	}
	norm_q = sqrt(norm_q);

	float lambda = M_ / norm_q;
	float* xbox_query = new float[xbox_dim_];
	for (int i = 0; i < xbox_dim_; i++) {
		if (i < dim_) {
			if (used_new_transform) xbox_query[i] = lambda * query[i];
			else xbox_query[i] = query[i];
		}
		else {
			xbox_query[i] = 0.0f;
		}
	}

	// -------------------------------------------------------------------------
	//  Perform knn search via qalsh
	// -------------------------------------------------------------------------
	int xbox_top_k = top_k;
	MinK_List *nn_list = new MinK_List(xbox_top_k);

	num_of_verf += lsh_->knn(MAXREAL, xbox_query, xbox_top_k, nn_list);

	// -------------------------------------------------------------------------
	//  Compute inner product for candidates returned by qalsh
	// -------------------------------------------------------------------------
	for (int i = 0; i < xbox_top_k; i++) {
		int id = nn_list->ith_smallest_id(i);
		float ip = calc_inner_product(data_[id], query, dim_);

		list->insert(ip, id + 1);
	}
	num_of_verf += xbox_top_k;

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] xbox_query; xbox_query = NULL;
	delete nn_list; nn_list = NULL;

	return num_of_verf;
}

