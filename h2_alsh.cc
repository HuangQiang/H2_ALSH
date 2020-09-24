
#include "h2_alsh.h"

// -----------------------------------------------------------------------------
H2_ALSH::H2_ALSH(					// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const float **data, 				// input data
	const float **norm_d)				// l2-norm of data objects
	: n_pts_(n), dim_(d), ratio_(mip_ratio), data_(data), norm_d_(norm_d)
{
	// -------------------------------------------------------------------------
	//  sort data objects by their Euclidean norms under the ascending order
	// -------------------------------------------------------------------------
	Result *order = new Result[n];
	for (int i = 0; i < n; ++i) {
		order[i].key_ = norm_d_[i][0];
		order[i].id_  = i;			// data object id		
	}
	qsort(order, n, sizeof(Result), ResultCompDesc);

	h2_alsh_id_ = new int[n];
	for (int i = 0; i < n; ++i) h2_alsh_id_[i] = order[i].id_;

	M_ = order[0].key_;
	b_ = sqrt((pow(nn_ratio,4.0f) - 1) / (pow(nn_ratio,4.0f) - mip_ratio));

	// -------------------------------------------------------------------------
	//  divide datasets into blocks and build qalsh for each block
	// -------------------------------------------------------------------------
	float *h2_alsh_data = new float[d + 1];
	int   start = 0;

	while (start < n) {
		// divide one block
		float M = order[start].key_;
		float min_radius = M * b_, M_sqr = SQR(M);
		int   idx = start, cnt = 0;

		while (idx < n && order[idx].key_ >= min_radius) {
			++idx;
			if (++cnt >= MAX_BLOCK_NUM) break;
		}

		// build qalsh for this block
		Block *block  = new Block();
		block->n_pts_ = cnt;
		block->M_     = M;
		block->index_ = h2_alsh_id_ + start;

		if (cnt > N_THRESHOLD) {
			block->lsh_ = new QALSH(cnt, d + 1, nn_ratio);
			
			// build hash tables for qalsh
			int m = block->lsh_->m_;
			for (int i = 0; i < cnt; ++i) {
				// construct new format of data by h2_alsh transformation
				int id = block->index_[i];
				for (int j = 0; j < d; ++j) {
					h2_alsh_data[j] = data[id][j];
				}
				h2_alsh_data[d] = sqrt(M_sqr - SQR(norm_d[id][0]));

				// calc hash value for new format of data
				for (int j = 0; j < m; ++j) {
					float val = block->lsh_->calc_hash_value(j, h2_alsh_data);
					block->lsh_->tables_[j][i].id_  = i;
					block->lsh_->tables_[j][i].key_ = val;
				}
			}
			for (int i = 0; i < m; ++i) {
				qsort(block->lsh_->tables_[i], cnt, sizeof(Result), ResultComp);
			}
		}
		blocks_.push_back(block);
		start += cnt;
	}
	assert(start == n);
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] order;
	delete[] h2_alsh_data;
}

// -----------------------------------------------------------------------------
H2_ALSH::~H2_ALSH()					// destructor
{
	delete[] h2_alsh_id_; h2_alsh_id_ = NULL; 	
	for (auto block : blocks_) {
		delete block; block = NULL;
	}
	blocks_.clear(); blocks_.shrink_to_fit();
}

// -------------------------------------------------------------------------
void H2_ALSH::display()				// display parameters
{
	printf("Parameters of H2_ALSH:\n");
	printf("    n          = %d\n",   n_pts_);
	printf("    d          = %d\n",   dim_);
	printf("    c          = %.1f\n", ratio_);
	printf("    M          = %f\n",   M_);
	printf("    num_blocks = %d\n\n", (int) blocks_.size());
}

// -----------------------------------------------------------------------------
int H2_ALSH::kmip(					// c-k-AMIP search
	int   top_k,						// top-k value
	const float *query,					// input query
	const float *norm_q,				// l2-norm of query
	MaxK_List *list)					// top-k MIP results (return)  
{
	// -------------------------------------------------------------------------
	//  initialize parameters
	// -------------------------------------------------------------------------
	float kip   = MINREAL;
	float normq = norm_q[0];
	float *h2_alsh_query = new float[dim_ + 1];
	std::vector<int> cand;

	// -------------------------------------------------------------------------
	//  c-k-AMIP search
	// -------------------------------------------------------------------------
	for (auto block : blocks_) {
		int   *index = block->index_;
		int   n      = block->n_pts_;
		float M      = block->M_;
		if (M * normq <= kip) break;

		if (n <= N_THRESHOLD) {
			// -----------------------------------------------------------------
			//  MIP search by linear scan
			// -----------------------------------------------------------------
			for (int j = 0; j < n; ++j) {
				int id = index[j];
				if (norm_d_[id][0] * normq <= kip) break;
				
				float ip = calc_inner_product(dim_, kip, data_[id], 
					norm_d_[id], query, norm_q);
				kip = list->insert(ip, id + 1);
			}
		}
		else {
			// -----------------------------------------------------------------
			//  conduct c-k-ANN search by qalsh
			// -----------------------------------------------------------------
			float lambda = M / normq;
			float R = sqrt(2.0f * (M * M - lambda * kip));
			for (int j = 0; j < dim_; ++j) {
				h2_alsh_query[j] = lambda * query[j];
			}
			h2_alsh_query[dim_] = 0.0f;

			cand.clear();
			block->lsh_->knn(top_k, R, (const float *) h2_alsh_query, cand);

			// -----------------------------------------------------------------
			//  compute inner product for the candidates returned by qalsh
			// -----------------------------------------------------------------
			int size = (int) cand.size();
			for (int j = 0; j < size; ++j) {
				int id = index[cand[j]];

				if (norm_d_[id][0] * normq > kip) {
					float ip = calc_inner_product(dim_, kip, data_[id], 
						norm_d_[id], query, norm_q);
					kip = list->insert(ip, id + 1);
				}
			}
		}
	}
	delete[] h2_alsh_query; h2_alsh_query = NULL;

	return 0;
}
