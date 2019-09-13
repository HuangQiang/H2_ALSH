#include "headers.h"

// -----------------------------------------------------------------------------
QALSH::QALSH(						// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	appr_ratio_ = ratio;
	data_       = data;
	beta_       = (float) CANDIDATES / (float) n_pts_;
	delta_      = 1.0f / E;

	w_  = sqrt((8.0f * appr_ratio_ * appr_ratio_ * log(appr_ratio_))
			/ (appr_ratio_ * appr_ratio_ - 1.0f));
	p1_ = calc_p(w_ / 2.0f);
	p2_ = calc_p(w_ / (2.0f * appr_ratio_));

	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);
	float eta   = para1 / para2;

	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);
	m_     = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_     = (int) ceil((p1_ * para1 + p2_ * para2) * (para1 + para2) / para3);

	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
	int size = m_ * dim_;
	a_array_ = new float[size];
	for (int i = 0; i < size; ++i) { // chosen from N(0.0, 1.0)
		a_array_[i] = gaussian(0.0F, 1.0F);
	}
	
	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	freq_        = new int[n_pts_];
	lpos_        = new int[m_];
	rpos_        = new int[m_];
	checked_     = new bool[n_pts_];
	bucket_flag_ = new bool[m_];
	range_flag_  = new bool[m_];
	q_val_       = new float[m_];

	tables_      = new Result*[m_];
	for (int i = 0; i < m_; ++i) {
		tables_[i] = new Result[n_pts_];
		for (int j = 0; j < n_pts_; ++j) {
			tables_[i][j].id_ = j;
			tables_[i][j].key_ = calc_hash_value(i, data_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(Result), ResultComp);
	}
	// display();
}

// -----------------------------------------------------------------------------
QALSH::~QALSH()						// destructor
{
	delete[] a_array_;     a_array_ = NULL;
	delete[] freq_;        freq_ = NULL;
	delete[] lpos_;        lpos_ = NULL;
	delete[] rpos_;        rpos_ = NULL;
	delete[] checked_;     checked_ = NULL;
	delete[] bucket_flag_; bucket_flag_ = NULL;
	delete[] range_flag_;  range_flag_ = NULL;
	delete[] q_val_;       q_val_ = NULL;
	
	for (int i = 0; i < m_; ++i) {
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] tables_; tables_ = NULL;
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_p(			// calc probability
	float x)							// x = w / (2.0 * r)
{
	return new_cdf(x, 0.001f);		// cdf of [-x, x]
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_hash_value( // calc hash value
	int   table_id,						// table id
	const float *data)					// input data
{
	float result = 0.0f;
	int   base   = table_id * dim_;
	for (int i = 0; i < dim_; ++i) {
		result += a_array_[base + i] * data[i];
	}
	return result;
}

// -----------------------------------------------------------------------------
void QALSH::display()				// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("    n     = %d\n",    n_pts_);
	printf("    d     = %d\n",    dim_);
	printf("    ratio = %.2f\n",  appr_ratio_);
	printf("    w     = %0.4f\n", w_);
	printf("    p1    = %0.4f\n", p1_);
	printf("    p2    = %0.4f\n", p2_);
	printf("    alpha = %0.6f\n", alpha_);
	printf("    beta  = %0.6f\n", beta_);
	printf("    delta = %0.6f\n", delta_);
	printf("    m     = %d\n",    m_);
	printf("    l     = %d\n\n",  l_);
}

// -----------------------------------------------------------------------------
int QALSH::knn(						// c-k-ANN search
	int   top_k,						// top-k
	float R,							// limited search range
	const float *query,					// input query
	MinK_List *list)					// top-k NN results (return)
{
	int   candidates = CANDIDATES + top_k - 1; // candidate size
	float knn_dist   = MAXREAL;		// k-th ANN distance

	// -------------------------------------------------------------------------
	//  linear scan if the number of data is small enough
	// -------------------------------------------------------------------------
	if (candidates >= n_pts_) {
		float dist = -1.0f;
		for (int i = 0; i < n_pts_; ++i) {
			dist = calc_l2_sqr(dim_, data_[i], query);
			list->insert(dist, i);
		}
		return n_pts_;
	}

	// -------------------------------------------------------------------------
	//  initialize parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);
	memset(bucket_flag_, true, m_ * SIZEBOOL);
	memset(range_flag_, true, m_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		int pos = binary_search_pos(i, q_val_[i]);
		if (pos == 0) {
			lpos_[i] = -1;  rpos_[i] = pos;
		} else {
			lpos_[i] = pos; rpos_[i] = pos + 1;
		}
	}

	// -------------------------------------------------------------------------
	//  k-nn search via dynamic collision counting
	// -------------------------------------------------------------------------
	int   dist_cnt     = 0;			// number of candidates computation
	int   num_range    = 0;			// number of search range flag
	float radius       = 1.0f;		// search radius
	float bucket_width = radius * w_ / 2.0f;  // bucket width
	
	float range = -1.0f;			// limited search range
	if (R > MAXREAL - 1.0f) range = MAXREAL;
	else range = R * w_ / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-NN search
		// ---------------------------------------------------------------------
		int count = 0;
		while (num_bucket < m_ && num_range < m_) {
			float ldist = -1.0f;	// left  proj dist to query
			float rdist = -1.0f;	// right proj dist to query

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				// -------------------------------------------------------------
				//  step 2.1: scan the left part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos_[j] >= 0) {
						ldist = abs(q_val_[j] - tables_[j][lpos_[j]].key_);
					}
					if (ldist > bucket_width) break;
					if (ldist > range) break;

					int id = tables_[j][lpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						float dist = calc_l2_sqr(dim_, data_[id], query);
						knn_dist = list->insert(dist, id);

						if (++dist_cnt > candidates) break;
					}
					--lpos_[j];
					++count;
				}
				if (dist_cnt > candidates) break;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					rdist = MAXREAL;
					if (rpos_[j] < n_pts_) {
						rdist = abs(q_val_[j] - tables_[j][rpos_[j]].key_);
					}
					if (rdist > bucket_width) break;
					if (rdist > range) break;					

					int id = tables_[j][rpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						float dist = calc_l2_sqr(dim_, data_[id], query);
						knn_dist = list->insert(dist, id);

						if (++dist_cnt > candidates) break;
					}
					++rpos_[j];
					++count;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (ldist > bucket_width && rdist > bucket_width) {
					bucket_flag_[j] = false;
					++num_bucket;
				}
				if (ldist > range && rdist > range) {
					bucket_flag_[j] = false;
					++num_bucket;
					if (range_flag_[j]) {
						range_flag_[j] = false;
						++num_range;
					}
				}
			}
			if (dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop condition T1 and T2
		// ---------------------------------------------------------------------
		if (sqrt(knn_dist) < appr_ratio_ * radius) break;
		if (num_range >= m_ || dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = appr_ratio_ * radius;
		bucket_width = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH::binary_search_pos(		// binary search position
	int   table_id,						// hash table is
	float value)						// hash value
{
	int left = 0;
	int right = n_pts_ - 1;
	int mid = 0;

	while (left < right) {
		mid = (left + right + 1) / 2;
		if (fabs(tables_[table_id][mid].key_ - value) < FLOATZERO) {
			return mid;
		}
		if (tables_[table_id][mid].key_ < value) left = mid;
		else right = mid - 1;
	}
	assert(left >= 0 && left < n_pts_);

	return left;
}
