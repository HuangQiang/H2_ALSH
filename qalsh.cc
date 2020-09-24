#include "qalsh.h"

// -----------------------------------------------------------------------------
QALSH::QALSH(						// constructor
	int   n,							// number of data objects
	int   d,							// dimension of data objects
	float ratio)						// approximation ratio
	: n_(n), d_(d), ratio_(ratio)
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	w_ = sqrt((8.0f * ratio * ratio * log(ratio)) / (ratio * ratio - 1.0f));
	
	float p1    = calc_p(w_ / 2.0f);
	float p2    = calc_p(w_ / (2.0f * ratio));
	float beta  = (float) CANDIDATES / n;
	float delta = 1.0f / E;

	float para1 = sqrt(log(2.0f / beta));
	float para2 = sqrt(log(1.0f / delta));
	float para3 = 2.0f * (p1 - p2) * (p1 - p2);
	float eta   = para1 / para2;
	float alpha = (eta * p1 + p2) / (1.0f + eta);
	
	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha * m_);

	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) { // chosen from N(0.0, 1.0)
		a_[i] = new float[d];
		for (int j = 0; j < d; ++j) {
			a_[i][j] = gaussian(0.0F, 1.0F);
		}
	}

	// -------------------------------------------------------------------------
	//  allocate space for hash tables
	// -------------------------------------------------------------------------
	tables_ = new Result*[m_];
	for (int i = 0; i < m_; ++i) tables_[i] = new Result[n];
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_p(			// calc probability
	float x)							// x = w / (2.0 * r)
{
	return new_cdf(x, 0.001f);		// cdf of [-x, x]
}

// -----------------------------------------------------------------------------
QALSH::~QALSH()						// destructor
{
	for (int i = 0; i < m_; ++i) {
		delete[] a_[i]; a_[i] = NULL;
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] a_; a_ = NULL;
	delete[] tables_; tables_ = NULL;
}

// -----------------------------------------------------------------------------
float QALSH::calc_hash_value(		// calc hash value
	int   tid,							// table id
	const float *data)					// input data
{
	return calc_inner_product(d_, a_[tid], data);
}

// -----------------------------------------------------------------------------
void QALSH::display()				// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("    n     = %d\n",   n_);
	printf("    d     = %d\n",   d_);
	printf("    c0    = %.1f\n", ratio_);
	printf("    w     = %f\n",   w_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n",   l_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int QALSH::knn(						// c-k-ANN search
	int   top_k,						// top-k
	float R,							// limited search range
	const float *query,					// input query
	std::vector<int> &cand)				// NN candidates (return)
{
	int   *freq        = new int[n_];
	int   *lpos        = new int[m_];
	int   *rpos        = new int[m_];
	bool  *checked     = new bool[n_];
	bool  *bucket_flag = new bool[m_];
	bool  *range_flag  = new bool[m_];
	float *q_val       = new float[m_];
	
	// -------------------------------------------------------------------------
	//  initialize parameters
	// -------------------------------------------------------------------------
	memset(freq,        0,     n_ * SIZEFLOAT);
	memset(checked,     false, n_ * SIZEBOOL);
	memset(range_flag,  true,  m_ * SIZEBOOL);
	
	Result tmp;
	Result *table = NULL;
	for (int i = 0; i < m_; ++i) {
		tmp.key_= calc_inner_product(d_, (const float *) a_[i], query);
		q_val[i] = tmp.key_;

		table = tables_[i];
		int pos = std::lower_bound(table, table+n_, tmp, cmp) - table;
		if (pos <= 0) {
			lpos[i] = -1; rpos[i] = pos;
		}
		else {
			lpos[i] = pos - 1; rpos[i] = pos;
		}
	}

	// -------------------------------------------------------------------------
	//  k-nn search via dynamic collision counting
	// -------------------------------------------------------------------------
	int   candidates = CANDIDATES + top_k - 1; // candidate size
	int   cand_cnt   = 0;			// candidate counter
	int   num_range  = 0;			// number of search range flag

	float radius = 1.0f;			// search radius
	float width  = radius * w_ / 2.0f;	// bucket width
	float range  = R > MAXREAL-1.0f ? MAXREAL : R * w_ / 2.0f; // search range

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-NN search
		// ---------------------------------------------------------------------
		while (num_bucket < m_ && num_range < m_) {
			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag[j]) continue;

				table = tables_[j];
				float q_v = q_val[j], ldist = -1.0f, rdist = -1.0f;
				// -------------------------------------------------------------
				//  step 2.1: scan the left part of hash table
				// -------------------------------------------------------------
				int cnt = 0;
				int pos = lpos[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (pos >= 0) {
						ldist = fabs(q_v - table[pos].key_);
					}
					else break;
					if (ldist > width || ldist > range) break;

					int id = table[pos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						cand.push_back(id);

						if (++cand_cnt >= candidates) break;
					}
					--pos; ++cnt;
				}
				if (cand_cnt >= candidates) break;
				lpos[j] = pos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				pos = rpos[j];
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (pos < n_) {
						rdist = fabs(q_v - table[pos].key_);
					}
					else break;
					if (rdist > width || rdist > range) break;

					int id = table[pos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						cand.push_back(id);

						if (++cand_cnt >= candidates) break;
					}
					++pos; ++cnt;
				}
				if (cand_cnt >= candidates) break;
				rpos[j] = pos;

				// -------------------------------------------------------------
				//  step 2.3: whether this bucket width is finished scanned
				// -------------------------------------------------------------
				if (ldist > width && rdist > width) {
					bucket_flag[j] = false;
					if (++num_bucket > m_) break;
				}
				if (ldist > range && rdist > range) {
					if (bucket_flag[j]) {
						bucket_flag[j] = false;
						if (++num_bucket > m_) break;
					}
					if (range_flag[j]) {
						range_flag[j] = false;
						if (++num_range > m_) break;
					}
				}
			}
			if (num_bucket > m_ || num_range > m_ || cand_cnt >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop condition
		// ---------------------------------------------------------------------
		if (num_range >= m_ || cand_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = ratio_ * radius;
		width  = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] freq;        freq        = NULL;
	delete[] lpos;        lpos        = NULL;
	delete[] rpos;        rpos        = NULL;
	
	delete[] checked;     checked     = NULL;
	delete[] bucket_flag; bucket_flag = NULL;
	delete[] range_flag;  range_flag  = NULL;
	delete[] q_val;       q_val       = NULL;

	return 0;
}
