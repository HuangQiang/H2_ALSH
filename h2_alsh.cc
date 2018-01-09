#include "headers.h"


// -----------------------------------------------------------------------------
H2_ALSH::H2_ALSH()					// constructor
{
	n_pts_ = -1;
	dim_ = -1;
	nn_ratio_ = -1.0f;
	mip_ratio_ = -1.0f;
	data_ = NULL;

	N_ = -1;
	cmpr_ratio_ = -1.0f;
	M_ = -1.0f;

	h2_alsh_dim_ = -1;
	h2_alsh_data_ = NULL;
	block_count_ = 0;
}

// -----------------------------------------------------------------------------
H2_ALSH::~H2_ALSH()					// destructor
{
	if (h2_alsh_data_ != NULL) {
		for (int i = 0; i < n_pts_; i++) {
			delete[] h2_alsh_data_[i]; h2_alsh_data_[i] = NULL;
		}
		delete[] h2_alsh_data_; h2_alsh_data_ = NULL;
	}

	for (int i = 0; i < block_count_; i++) {
		delete blocks_[i]; blocks_[i] = NULL;
	}
	blocks_.clear();
}

// -----------------------------------------------------------------------------
void H2_ALSH::init(					// init the parameters
	int n,								// number of data
	int d,								// dimension of data
	float nn_ratio,						// approximation ratio for NN
	float mip_ratio,					// approximation ratio for MIP
	float** data)						// input data
{
	n_pts_ = n;	
	dim_ = d;
	nn_ratio_ = nn_ratio;
	mip_ratio_ = mip_ratio;
	data_ = data;

	h2_alsh_dim_ = d + 1;
	N_ = CANDIDATES;

	cmpr_ratio_ = sqrt((pow(nn_ratio_, 4.0f) - 1) /
		(pow(nn_ratio_, 4.0f) - mip_ratio_));
	
	M_ = -1.0f;
	block_count_ = -1;
	h2_alsh_data_ = NULL;

	pre_processing();
}

// -----------------------------------------------------------------------------
int H2_ALSH::pre_processing()		// pre-processing of data
{
	// -------------------------------------------------------------------------
	//  sort the data by norm
	// -------------------------------------------------------------------------
	pair<float, int>* pq_nodes = new pair<float, int>[n_pts_];

	for (int i = 0; i < n_pts_; i++) {
		pq_nodes[i].second = i;

		pq_nodes[i].first = 0.0f;
		for (int j = 0; j < dim_; j++) {
			pq_nodes[i].first += data_[i][j] * data_[i][j];
		}
		pq_nodes[i].first = sqrt(pq_nodes[i].first);
	}
	sort(pq_nodes, pq_nodes + n_pts_);
	M_ = pq_nodes[n_pts_ - 1].first;

	// -------------------------------------------------------------------------
	//  construct new data and indexing
	// -------------------------------------------------------------------------
	printf("Construct H2_ALSH Data: ");
	h2_alsh_data_ = new float*[n_pts_];
	block_count_ = 0;

	//int block_index = 0;
	int node_index = n_pts_ - 1;	// id for pq_node (descreasing)
	int data_index = 0;				// id for h2_alsh_data (increasing)

	while (node_index >= 0) {
		int data_count = 0;
		float M = pq_nodes[node_index].first;
		float m = M * cmpr_ratio_;

		while (node_index >= 0) {
			if (data_count >= 5000) break;
		
			if (data_count >= N_ && pq_nodes[node_index].first < m) {
				break;
			}

			int id = pq_nodes[node_index].second;
			h2_alsh_data_[data_index] = new float[h2_alsh_dim_];

			for (int j = 0; j < h2_alsh_dim_; j++) {
				if (j < dim_) {
					h2_alsh_data_[data_index][j] = data_[id][j];
				}
				else {
					h2_alsh_data_[data_index][j] = sqrt(M * M -
						pq_nodes[node_index].first * pq_nodes[node_index].first);
				}
			}

			node_index--;
			data_index++;
			data_count++;
		}

		Block *block = new Block();
		block->n_pts_ = data_count;
		block->M_ = M;
		block->index_ = new int[block->n_pts_];
		for (int i = 0; i < block->n_pts_; i++) {
			block->index_[i] = pq_nodes[node_index + data_count - i].second;
		}
		//display_block_params(block_count_, block);

		if (data_count > N_) {
			int start_index = data_index - data_count;
			//cout << start_index << " " << data_count << endl;
			if (indexing(start_index, data_count, block_count_, block)) {
				return 1;
			}
		}
		block_count_++;
		blocks_.push_back(block);
	}
	display_params();
	
	delete[] pq_nodes; pq_nodes = NULL;
	return 0;
}

// -----------------------------------------------------------------------------
int H2_ALSH::indexing(				// indexing the new data
	int start_index,
	int n_pts,
	int block_index,
	Block* block)
{
	block->lsh_ = new QALSH_Col(h2_alsh_data_ + start_index, n_pts, 
		h2_alsh_dim_, nn_ratio_);
	
	if (block->lsh_ != NULL) return 0;
	else return 1;
}

// -----------------------------------------------------------------------------
void H2_ALSH::display_params()		// display parameters
{
	printf("Parameters of H2_ALSH:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    c = %.2f\n", nn_ratio_);
	printf("    c' = %.2f\n", mip_ratio_);
	printf("    M = %.2f\n", M_);
	printf("    O = %d\n\n", block_count_);
}

// -----------------------------------------------------------------------------
void H2_ALSH::display_block_params(	// display parameters for block
	int block_index,					// block id
	Block* block)						// input block
{
	printf("Parameters of block %d:\n", block_index);
	printf("    n = %d\n", block->n_pts_);
	printf("    M = %.2f\n\n", block->M_);
}

// -----------------------------------------------------------------------------
int H2_ALSH::kmip(					// top-k approximate mip search
	float* query,						// input query
	int top_k,							// top-k value
	MaxK_List* list)					// top-k mip results
{
	int num_of_verf = 0;			// num of verification (NN and MIP calc)

	// -------------------------------------------------------------------------
	//  Construct H2_ALSH query
	// -------------------------------------------------------------------------
	float* h2_alsh_query = new float[h2_alsh_dim_];

	float norm_q = 0.0f;
	for (int i = 0; i < dim_; i++) {
		norm_q += query[i] * query[i];
	}
	norm_q = sqrt(norm_q);

	// -------------------------------------------------------------------------
	//  kMIP search
	// -------------------------------------------------------------------------
	int h2_alsh_top_k = top_k;
	MinK_List *nn_list = new MinK_List(h2_alsh_top_k);

	int block_id = 0;
	for (block_id = 0; block_id < block_count_; block_id++) {
		//cout << " block id = " << block_id << " " << blocks_[block_id]->M_ << endl;

		float mip = blocks_[block_id]->M_ * norm_q;
		if (mip <= list->min_key()) { // * mip_ratio_
			break;
		}

		if (blocks_[block_id]->n_pts_ <= N_) {
			for (int j = 0; j < blocks_[block_id]->n_pts_; j++) {
				int id = blocks_[block_id]->index_[j];
				float ip = calc_inner_product(data_[id], query, dim_);

				list->insert(ip, id + 1);
			}
			num_of_verf += blocks_[block_id]->n_pts_;
		}
		else {
			float lambda = blocks_[block_id]->M_ / norm_q;
			for (int i = 0; i < h2_alsh_dim_; i++) {
				if (i < dim_) {
					h2_alsh_query[i] = lambda * query[i];
				}
				else {
					h2_alsh_query[i] = 0.0f;
				}
			}

			// -----------------------------------------------------------------
			//  Perform knn search via qalsh
			// -----------------------------------------------------------------
			nn_list->reset();

			float R = 2.0f * blocks_[block_id]->M_ * blocks_[block_id]->M_;
			R -= 2.0f * lambda * list->min_key();
			R = sqrt(R);
			
			num_of_verf += blocks_[block_id]->lsh_->knn(R, h2_alsh_query, 
				h2_alsh_top_k, nn_list);

			// -----------------------------------------------------------------
			//  Compute inner product for candidates returned by qalsh
			// -----------------------------------------------------------------
			int size = (int)nn_list->size();
			for (int i = 0; i < size; i++) {
				int nn_id = nn_list->ith_smallest_id(i);

				int id = blocks_[block_id]->index_[nn_id];
				float ip = calc_inner_product(data_[id], query, dim_);
				list->insert(ip, id + 1);
			}
			num_of_verf += size;
		}
	}

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete[] h2_alsh_query; h2_alsh_query = NULL;
	delete nn_list; nn_list = NULL;

	return num_of_verf;
}
