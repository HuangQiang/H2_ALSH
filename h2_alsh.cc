#include "headers.h"

// -----------------------------------------------------------------------------
H2_ALSH::H2_ALSH(					// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	FILE  *fp,							// output file pointer
	const float **data)					// input data
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	n_pts_     = n;	
	dim_       = d;
	nn_ratio_  = nn_ratio;
	mip_ratio_ = mip_ratio;
	data_      = data;

	b_ = sqrt((pow(nn_ratio_,4.0f) - 1) / (pow(nn_ratio_,4.0f) - mip_ratio_));

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;	

	// -------------------------------------------------------------------------
	//  display parameters
	// -------------------------------------------------------------------------
	printf("Parameters of H2_ALSH:\n");
	printf("    n          = %d\n",   n_pts_);
	printf("    d          = %d\n",   dim_);
	printf("    c          = %.2f\n", nn_ratio_);
	printf("    c0         = %.2f\n", mip_ratio_);
	printf("    M          = %.2f\n", M_);
	printf("    num_blocks = %d\n\n", num_blocks_);
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	fprintf(fp, "n          = %d\n", n_pts_);
	fprintf(fp, "d          = %d\n", dim_);
	fprintf(fp, "c          = %f\n", nn_ratio_);
	fprintf(fp, "c0         = %f\n", mip_ratio_);
	fprintf(fp, "M          = %f\n", M_);
	fprintf(fp, "num_blocks = %d\n", num_blocks_);
	fprintf(fp, "index_time = %f Seconds\n\n", indexing_time);
}

// -----------------------------------------------------------------------------
void H2_ALSH::bulkload()			// bulkloading
{
	// -------------------------------------------------------------------------
	//  sort data objects by their Euclidean norms under the ascending order
	// -------------------------------------------------------------------------
	vector<pair<float, int> > norm(n_pts_);
	for (int i = 0; i < n_pts_; ++i) {
		norm[i].second = i;			// data object id
		norm[i].first  = sqrt(calc_inner_product(dim_, data_[i], data_[i]));
	}
	sort(norm.begin(), norm.end());
	M_ = norm[n_pts_ - 1].first;

	// -------------------------------------------------------------------------
	//  construct new data
	// -------------------------------------------------------------------------
	h2_alsh_data_ = new float*[n_pts_];
	num_blocks_ = 0;

	int node_index = n_pts_ - 1;	// id for <norm> (descreasing)
	int data_index = 0;				// id for <h2_alsh_data> (increasing)

	while (node_index >= 0) {
		int   data_count = 0;
		float M = norm[node_index].first;
		float m = M * b_;

		while (node_index >= 0) {
			if (data_count >= MAX_BLOCK_NUM) break;

			int norm_d = norm[node_index].first;
			int id     = norm[node_index].second;
			if (data_count >= CANDIDATES && norm_d < m) {
				break;
			}

			float *data = new float[dim_ + 1];
			for (int i = 0; i < dim_; ++i) {
				data[i] = data_[id][i];
			}
			data[dim_]= sqrt(M * M - norm_d * norm_d);
			h2_alsh_data_[data_index] = data; 

			--node_index;
			++data_index;
			++data_count;
		}

		Block *block = new Block();
		block->n_pts_ = data_count;
		block->M_     = M;
		block->index_ = new int[data_count];

		int base = node_index + data_count;
		for (int i = 0; i < data_count; ++i) {
			block->index_[i] = norm[base - i].second;
		}

		if (data_count > CANDIDATES) {
			int start_index = data_index - data_count;
			block->lsh_ = new QALSH(data_count, dim_ + 1, nn_ratio_, 
				(const float **) h2_alsh_data_ + start_index);
		}
		blocks_.push_back(block);
		num_blocks_++;
	}
}

// -----------------------------------------------------------------------------
H2_ALSH::~H2_ALSH()					// destructor
{
	if (h2_alsh_data_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] h2_alsh_data_[i]; h2_alsh_data_[i] = NULL;
		}
		delete[] h2_alsh_data_; h2_alsh_data_ = NULL;
	}

	if (!blocks_.empty()) {
		for (int i = 0; i < num_blocks_; ++i) {
			delete blocks_[i]; blocks_[i] = NULL;
		}
		blocks_.clear(); blocks_.shrink_to_fit();
	}
}

// -----------------------------------------------------------------------------
int H2_ALSH::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  initialize parameters
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));
	int   *index = NULL;
	Block *block = NULL;

	float *h2_alsh_query = new float[dim_ + 1];
	MinK_List *nn_list = new MinK_List(top_k);

	// -------------------------------------------------------------------------
	//  c-k-AMIP search
	// -------------------------------------------------------------------------
	for (int i = 0; i < num_blocks_; ++i) {
		block = blocks_[i];
		index = block->index_;
		if (block->M_ * norm_q <= list->min_key()) break;

		int n = block->n_pts_;
		if (n <= CANDIDATES) {
			// -----------------------------------------------------------------
			//  MIP search by linear scan
			// -----------------------------------------------------------------
			for (int j = 0; j < n; ++j) {
				int   id = index[j];
				float ip = calc_inner_product(dim_, data_[id], query);

				list->insert(ip, id + 1);
			}
		}
		else {
			float lambda = block->M_ / norm_q;
			for (int j = 0; j < dim_; ++j) {
				h2_alsh_query[j] = lambda * query[j];
			}
			h2_alsh_query[dim_] = 0.0f;

			// -----------------------------------------------------------------
			//  conduct c-k-ANN search by qalsh
			// -----------------------------------------------------------------
			nn_list->reset();
			float R = sqrt(2.0f * (block->M_*block->M_ - lambda*list->min_key()));
			block->lsh_->knn(top_k, R, (const float *) h2_alsh_query, nn_list);

			// -----------------------------------------------------------------
			//  compute inner product for the candidates returned by qalsh
			// -----------------------------------------------------------------
			int size = (int) nn_list->size();
			for (int j = 0; j < size; ++j) {
				int   id = index[nn_list->ith_id(j)];
				float ip = calc_inner_product(dim_, data_[id], query);
				
				list->insert(ip, id + 1);
			}
		}
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] h2_alsh_query; h2_alsh_query = NULL;
	delete nn_list; nn_list = NULL;

	return 0;
}
