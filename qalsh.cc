#include <algorithm>
#include <cstring>

#include "def.h"
#include "random.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"

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
	beta_       = (float) CANDIDATES / n;
	delta_      = 1.0f / E;

	w_  = sqrt((8.0f * ratio * ratio * log(ratio)) / (ratio * ratio - 1.0f));
	p1_ = calc_p(w_ / 2.0f);
	p2_ = calc_p(w_ / (2.0f * ratio));

	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);
	float eta   = para1 / para2;

	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);
	m_     = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_     = (int) ceil(alpha_ * m_);

	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) { // chosen from N(0.0, 1.0)
		a_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			a_[i][j] = gaussian(0.0F, 1.0F);
		}
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

	tables_ = new Result*[m_];
	for (int i = 0; i < m_; ++i) {
		tables_[i] = new Result[n_pts_];

		const float *a = a_[i];
		for (int j = 0; j < n_pts_; ++j) {
			tables_[i][j].id_  = j;
			tables_[i][j].key_ = calc_inner_product(dim_, a, data_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(Result), ResultComp);
	}
}

// -----------------------------------------------------------------------------
QALSH::~QALSH()						// destructor
{
	delete[] freq_;        freq_        = NULL;
	delete[] lpos_;        lpos_        = NULL;
	delete[] rpos_;        rpos_        = NULL;
	delete[] checked_;     checked_     = NULL;
	delete[] bucket_flag_; bucket_flag_ = NULL;
	delete[] range_flag_;  range_flag_  = NULL;
	delete[] q_val_;       q_val_       = NULL;
	
	for (int i = 0; i < m_; ++i) {
		delete[] a_[i];      a_[i]      = NULL;
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] a_;      a_      = NULL;
	delete[] tables_; tables_ = NULL;
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_p(			// calc probability
	float x)							// x = w / (2.0 * r)
{
	return new_cdf(x, 0.001f);		// cdf of [-x, x]
}

// -----------------------------------------------------------------------------
void QALSH::display()				// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    c0    = %.1f\n", appr_ratio_);
	printf("    w     = %f\n",   w_);
	printf("    p1    = %f\n",   p1_);
	printf("    p2    = %f\n",   p2_);
	printf("    alpha = %f\n",   alpha_);
	printf("    beta  = %f\n",   beta_);
	printf("    delta = %f\n",   delta_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n\n", l_);
}

// -----------------------------------------------------------------------------
int QALSH::knn(						// c-k-ANN search
	int   top_k,						// top-k
	float R,							// limited search range
	const float *query,					// input query
	std::vector<int> &cand)				// NN candidates (return)
{
	int candidates = CANDIDATES + top_k - 1; // candidate size
	// float kdist = MAXREAL;			// k-th ANN distance
	
	// -------------------------------------------------------------------------
	//  initialize parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);
	memset(bucket_flag_, true, m_ * SIZEBOOL);
	memset(range_flag_, true, m_ * SIZEBOOL);
	
	Result tmp;
	Result *table = NULL;
	for (int i = 0; i < m_; ++i) {
		tmp.key_= calc_inner_product(dim_, (const float *) a_[i], query);
		q_val_[i] = tmp.key_;

		table = tables_[i];
		int pos = std::lower_bound(table, table+n_pts_, tmp, cmp) - table;
		if (pos <= 0) {
			lpos_[i] = -1; rpos_[i] = pos;
		}
		else {
			lpos_[i] = pos - 1; rpos_[i] = pos;
		}
	}

	// -------------------------------------------------------------------------
	//  k-nn search via dynamic collision counting
	// -------------------------------------------------------------------------
	int dist_cnt  = 0;				// number of candidates computation
	int num_range = 0;				// number of search range flag

	float radius = 1.0f;			// search radius
	float bucket = radius * w_ / 2.0f; // bucket width
	float range  = -1.0f;			// limited search range
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
		int cnt = -1, pos = -1, id = -1;
		while (num_bucket < m_ && num_range < m_) {
			float ldist = -1.0f;	// left  proj dist to query
			float rdist = -1.0f;	// right proj dist to query
			float q_val = -1.0f;	// hash value of 
			float dist  = -1.0f;	// l2-sqr dist

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				table = tables_[j];
				q_val = q_val_[j];
				// -------------------------------------------------------------
				//  step 2.1: scan the left part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				pos = lpos_[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (pos >= 0) {
						ldist = fabs(q_val - table[pos].key_);
					}
					if (ldist > bucket || ldist > range) break;

					id = table[pos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						cand.push_back(id);
						// dist = calc_l2_sqr(dim_, kdist, data_[id], query);
						// kdist = list->insert(dist, id);

						if (++dist_cnt >= candidates) break;
					}
					--pos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				lpos_[j] = pos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				pos = rpos_[j];
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (pos < n_pts_) {
						rdist = fabs(q_val - table[pos].key_);
					}
					if (rdist > bucket || rdist > range) break;

					id = table[pos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						cand.push_back(id);
						// dist = calc_l2_sqr(dim_, kdist, data_[id], query);
						// kdist = list->insert(dist, id);

						if (++dist_cnt >= candidates) break;
					}
					++pos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				rpos_[j] = pos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (ldist > bucket && rdist > bucket) {
					bucket_flag_[j] = false;
					if (++num_bucket > m_) break;
				}
				if (ldist > range && rdist > range) {
					if (bucket_flag_[j]) {
						bucket_flag_[j] = false;
						if (++num_bucket > m_) break;
					}
					if (range_flag_[j]) {
						range_flag_[j] = false;
						if (++num_range > m_) break;
					}
				}
			}
			if (dist_cnt >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop condition T1 and T2
		// ---------------------------------------------------------------------
		// if (sqrt(kdist) < appr_ratio_ * radius) break;
		if (num_range >= m_ || dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = appr_ratio_ * radius;
		bucket = radius * w_ / 2.0f;
	}
	return 0;
}
