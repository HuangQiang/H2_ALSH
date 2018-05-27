#include "headers.h"

// -----------------------------------------------------------------------------
H2_ALSH::H2_ALSH()					// default constructor
{
	n_pts_        = -1;
	dim_          = -1;
	nn_ratio_     = -1.0f;
	mip_ratio_    = -1.0f;
	data_         = NULL;

	b_            = -1.0f;
	M_            = -1.0f;
	h2_alsh_dim_  = -1;
	h2_alsh_data_ = NULL;
	num_blocks_   = 0;
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
void H2_ALSH::build(				// build index
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const float **data)					// input data
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_     = n;	
	dim_       = d;
	nn_ratio_  = nn_ratio;
	mip_ratio_ = mip_ratio;
	data_      = data;

	h2_alsh_dim_ = d + 1;
	b_ = sqrt((pow(nn_ratio_,4.0f) - 1) / (pow(nn_ratio_,4.0f) - mip_ratio_));

	// -------------------------------------------------------------------------
	//  build index
	// -------------------------------------------------------------------------
	bulkload();
	display();
}

// -----------------------------------------------------------------------------
int H2_ALSH::bulkload()				// bulkloading
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
	printf("Construct H2_ALSH Data\n\n");
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
			if (data_count >= CANDIDATES && norm[node_index].first < m) {
				break;
			}
			int id = norm[node_index].second;
			h2_alsh_data_[data_index] = new float[h2_alsh_dim_];

			for (int j = 0; j < h2_alsh_dim_; ++j) {
				if (j < dim_) {
					h2_alsh_data_[data_index][j] = data_[id][j];
				}
				else {
					h2_alsh_data_[data_index][j] = sqrt(M * M -
						norm[node_index].first * norm[node_index].first);
				}
			}
			--node_index;
			++data_index;
			++data_count;
		}

		Block *block = new Block();
		block->n_pts_ = data_count;
		block->M_     = M;
		block->index_ = new int[block->n_pts_];
		for (int i = 0; i < block->n_pts_; ++i) {
			block->index_[i] = norm[node_index + data_count - i].second;
		}

		if (data_count > CANDIDATES) {
			int start_index = data_index - data_count;
			block->lsh_ = new QALSH(data_count, h2_alsh_dim_, nn_ratio_, 
				(const float **) h2_alsh_data_ + start_index);
		}
		num_blocks_++;
		blocks_.push_back(block);
	}

	return 0;
}

// -----------------------------------------------------------------------------
void H2_ALSH::display()				// display parameters
{
	printf("Parameters of H2_ALSH:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    c          = %.2f\n", nn_ratio_);
	printf("    c0         = %.2f\n", mip_ratio_);
	printf("    M          = %.2f\n", M_);
	printf("    num_blocks = %d\n", num_blocks_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int H2_ALSH::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MIP results (return) 
{
	// -------------------------------------------------------------------------
	//  calc the Euclidean norm of input query
	// -------------------------------------------------------------------------
	float norm_q = sqrt(calc_inner_product(dim_, query, query));

	// -------------------------------------------------------------------------
	//  c-k-AMIP search
	// -------------------------------------------------------------------------
	float *h2_alsh_query = new float[h2_alsh_dim_];
	MinK_List *nn_list = new MinK_List(top_k);

	for (int block_id = 0; block_id < num_blocks_; ++block_id) {
		float ub = blocks_[block_id]->M_ * norm_q;
		if (ub <= list->min_key()) break;

		if (blocks_[block_id]->n_pts_ <= CANDIDATES) {
			// -----------------------------------------------------------------
			//  MIP search by linear scan
			// -----------------------------------------------------------------
			for (int j = 0; j < blocks_[block_id]->n_pts_; ++j) {
				int id = blocks_[block_id]->index_[j];
				float ip = calc_inner_product(dim_, data_[id], query);

				list->insert(ip, id + 1);
			}
		}
		else {
			float lambda = blocks_[block_id]->M_ / norm_q;
			for (int i = 0; i < h2_alsh_dim_; ++i) {
				if (i < dim_) h2_alsh_query[i] = lambda * query[i];
				else h2_alsh_query[i] = 0.0f;
			}

			// -----------------------------------------------------------------
			//  conduct c-k-ANN search by qalsh
			// -----------------------------------------------------------------
			nn_list->reset();

			float R = 2.0f * blocks_[block_id]->M_ * blocks_[block_id]->M_;
			R -= 2.0f * lambda * list->min_key();
			R = sqrt(R);
			
			blocks_[block_id]->lsh_->knn(top_k, R, (const float *) h2_alsh_query, 
				nn_list);

			// -----------------------------------------------------------------
			//  compute inner product for the candidates returned by qalsh
			// -----------------------------------------------------------------
			int size = (int) nn_list->size();
			for (int i = 0; i < size; ++i) {
				int nn_id = nn_list->ith_id(i);

				int id = blocks_[block_id]->index_[nn_id];
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
