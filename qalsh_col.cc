#include "headers.h"


// -----------------------------------------------------------------------------
//  QALSH_Col: memory version of query-aware locality-sensitive hashing based 
//  on dynamic collision counting for approximate nearest neighbor search.
// -----------------------------------------------------------------------------
QALSH_Col::QALSH_Col(				// constructor
	float** points,						// data points
	int n,								// cardinality
	int d,								// dimensionality
	float ratio)						// approximation ratio
{
	points_ = points;
	n_pts_ = n;
	dim_ = d;
	appr_ratio_ = ratio;

	calc_paras();

	gen_hash_func();
	build_hashtables();
	//display_paras();
}

// -----------------------------------------------------------------------------
QALSH_Col::~QALSH_Col()				// destructor
{
	for (int i = 0; i < m_; i++) {	// release <tables_>
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] tables_; tables_ = NULL;

	if (a_array_) {					// release <a_array_>
		delete[] a_array_; a_array_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void QALSH_Col::calc_paras()		// init parameters
{
	beta_ = 100.0f / n_pts_;		// init <beta_>
	delta_ = 1.0f / E;				// init <delta_>

									// calc <w_> to minimize <m_>
	w_ = sqrt((8.0f * appr_ratio_ * appr_ratio_ * log(appr_ratio_))
			/ (appr_ratio_ * appr_ratio_ - 1.0f));

	p1_ = calc_p(w_ / 2.0f);		// calc <p1_> and <p2_>
	p2_ = calc_p(w_ / (2.0f * appr_ratio_));

	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);

	float eta = para1 / para2;		// calc <alpha_>
	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);

									// calc <m_> and <l_>
	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil((p1_ * para1 + p2_ * para2) * (para1 + para2) / para3);
}

// -----------------------------------------------------------------------------
float QALSH_Col::calc_p(			// calc probability
	float x)							// integral border
{
	return new_cdf(x, 0.001f);		// cdf of [-x, x]
}

// -----------------------------------------------------------------------------
void QALSH_Col::gen_hash_func()		// generate hash functions
{
	int size = m_ * dim_;
	a_array_ = new float[size];
	for (int i = 0; i < size; i++) {// <a_array_> chosen from N(0.0, 1.0)
		a_array_[i] = gaussian(0.0F, 1.0F);
	}
}

// -----------------------------------------------------------------------------
void QALSH_Col::build_hashtables() 	// build hash tables
{
	tables_ = new HashValue*[m_];	// init <tables_>

	for (int i = 0; i < m_; i++) {
		tables_[i] = new HashValue[n_pts_];
		for (int j = 0; j < n_pts_; j++) {
			tables_[i][j].id_ = j;
			tables_[i][j].proj_ = calc_hash_value(i, points_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(HashValue), HashValueQsortComp);
	}
}

// -----------------------------------------------------------------------------
float QALSH_Col::calc_hash_value(	// calc hash value
	int table_id,						// hash table id
	float* point)						// one point
{
	float result = 0.0f;
	for (int i = 0; i < dim_; i++) {
		result += (a_array_[table_id * dim_ + i] * point[i]);
	}
	return result;
}

// -----------------------------------------------------------------------------
void QALSH_Col::display_paras()		// display parameters
{
	printf("Parameters of QALSH_Col:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    ratio      = %.2f\n", appr_ratio_);
	//printf("    w          = %0.4f\n", w_);
	//printf("    p1         = %0.4f\n", p1_);
	//printf("    p2         = %0.4f\n", p2_);
	//printf("    alpha      = %0.6f\n", alpha_);
	//printf("    beta       = %0.6f\n", beta_);
	//printf("    delta      = %0.6f\n", delta_);
	printf("    m          = %d\n", m_);
	printf("    l          = %d\n", l_);
	//printf("    candidates = %d\n", 100);
	printf("\n");
}


// -----------------------------------------------------------------------------
int QALSH_Col::knn(					// knn search with collision counting
	float R,							// limited search range
	float* query,						// query point
	int top_k,							// top-k value
	MinK_List* list)					// k-nn results
{
	// TODO: restrict the returned objects <= n_pts_
	if (top_k + 99 >= n_pts_) {
		float dist = -1.0f;
		for (int i = 0; i < n_pts_; i++) {
			dist = calc_l2_dist(points_[i], query, dim_);
			list->insert(dist, i);
		}
		
		return n_pts_;
	}

	int* frequency = new int[n_pts_];
	bool* checked = new bool[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		frequency[i] = 0;			// collision time
		checked[i] = false;			// check whether calc l2 dist
	}
	
	float* q_val = new float[m_];	// hash value of query
	bool* flag = new bool[m_];		// flag of each bucket
	int* lpos  = new int[m_];		// left position of hash table
	int* rpos  = new int[m_];		// right position of hash table

	for (int i = 0; i < m_; i++) {	// calc hash value
		q_val[i] = calc_hash_value(i, query);
		flag[i]  = true;
									// find position of queries
		int pos = binary_search_pos(i, q_val[i]);
		if (pos == 0) {
			lpos[i] = -1;  rpos[i] = pos;
		} else {
			lpos[i] = pos; rpos[i] = pos + 1;
		}
	}

	// -------------------------------------------------------------------------
	//  k-nn search via dynamic collision counting
	// -------------------------------------------------------------------------
	int cand_size = 99 + top_k;		// candidate size
	int num_cand = 0;				// number of candidates
	int flag_num = 0;				// flag number

	float radius = 1.0f;			// current radius
	float bucket_width = radius * w_ / 2.0f;
	float knn_dist = MAXREAL;		// k-th nn distance

	float left_dist = MAXREAL;
	float right_dist = MAXREAL;

	// -------------------------------------------------------------------------
	int range_num = 0;
	bool *range_flag = new bool[m_];
	for (int i = 0; i < m_; i++) {
		range_flag[i] = true;
	}

	float range_width = -1.0f;
	if (R > MAXREAL - 1.0f) range_width = MAXREAL;
	else range_width = R * w_ / 2.0f;
	// -------------------------------------------------------------------------
	
	while (true) {
		// ---------------------------------------------------------------------
		//  Step 1: initialization
		// ---------------------------------------------------------------------
		flag_num = 0;
		for (int j = 0; j < m_; j++) {
			flag[j] = true;
		}

		// ---------------------------------------------------------------------
		//  Step 2: (r,c)-nn search
		// ---------------------------------------------------------------------
		int count = 0;
		while (flag_num < m_ && range_num < m_) {
			for (int j = 0; j < m_; j++) {
				if (!flag[j]) continue;

				// -------------------------------------------------------------
				//  Step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					left_dist = MAXREAL;
					if (lpos[j] >= 0) {
						left_dist = abs(q_val[j] - tables_[j][lpos[j]].proj_);
					}
					if (left_dist > bucket_width) break;
					if (left_dist > range_width) break;

					int id = tables_[j][lpos[j]].id_;
					frequency[id]++;
					if (frequency[id] >= l_ && !checked[id]) {
						float dist = calc_l2_dist(points_[id], query, dim_);
						list->insert(dist, id);
						knn_dist = list->max_key();

						checked[id] = true;
						num_cand++;
						if (num_cand > cand_size) break;
					}
					lpos[j]--;		// move one step to left
					count++;
				}
				if (num_cand > cand_size) break;

				// -------------------------------------------------------------
				//  Step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					right_dist = MAXREAL;
					if (rpos[j] < n_pts_) {
						right_dist = abs(q_val[j] - tables_[j][rpos[j]].proj_);
					}
					if (right_dist > bucket_width) break;
					if (right_dist > range_width) break;

					int id = tables_[j][rpos[j]].id_;
					frequency[id]++;
					if (frequency[id] >= l_ && !checked[id]) {
						float dist = calc_l2_dist(points_[id], query, dim_);
						list->insert(dist, id);
						knn_dist = list->max_key();

						checked[id] = true;
						num_cand++;
						if (num_cand > cand_size) break;
					}
					rpos[j]++;		// move one step to right
					count++;
				}
				if (num_cand >= cand_size) break;

				// -------------------------------------------------------------
				//  Step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (left_dist > bucket_width && right_dist > bucket_width) {
					flag[j] = false;
					flag_num++;
				}
				if (left_dist > range_width && right_dist > range_width) { // TODO Qiang added
					flag[j] = false;
					flag_num++;
					if (range_flag[j]) {
						range_flag[j] = false;
						range_num++;
						//cout << range_num << endl;
					}
				}
			}
			if (num_cand >= cand_size) break;
		}

		// ---------------------------------------------------------------------
		//  Step 3: stop condition T1 and T2
		// ---------------------------------------------------------------------
		if (knn_dist < appr_ratio_ * radius || num_cand >= cand_size) {
			break;
		}
		if (range_num >= m_) break;

		// ---------------------------------------------------------------------
		//  Step 4: update radius
		// ---------------------------------------------------------------------
		radius = appr_ratio_ * radius;
		bucket_width = radius * w_ / 2.0f;
	}

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (frequency  != NULL || checked != NULL) {
		delete[] frequency; frequency = NULL;
		delete[] checked; checked = NULL;
	}
	if (q_val != NULL || lpos != NULL || rpos != NULL) {
		delete[] q_val; q_val = NULL;
		delete[] lpos;  lpos = NULL;
		delete[] rpos;  rpos = NULL;
	}

	return num_cand;
}


// -----------------------------------------------------------------------------
int QALSH_Col::binary_search_pos(	// binary search position
	int table_id,						// hash table is
	float value)						// hash value
{
	int left = 0;
	int right = n_pts_ - 1;
	int mid = 0;

	while (left < right) {
		mid = (left + right + 1) / 2;
		if (tables_[table_id][mid].proj_ == value) return mid;

		if (tables_[table_id][mid].proj_ <= value) left = mid;
		else right = mid - 1;
	}
									// check whether out of range
	if (left < 0 || left >= n_pts_) {
		error("QALSH_Col::binary_search_location error!", false);
	}
	return left;
}


// -----------------------------------------------------------------------------
int HashValueQsortComp(				// qsort assistant function
	const void* obj1,					// 1st object
	const void* obj2)					// 2nd object
{
	int ret = 0;
	HashValue* e1 = (HashValue *) obj1;
	HashValue* e2 = (HashValue *) obj2;

	if (e1->proj_ < e2->proj_) {
		ret = -1;
	} else if (e1->proj_ > e2->proj_) {
		ret = 1;
	} else {
		if (e1->id_ < e2->id_) ret = -1;
		else if (e1->id_ > e2->id_) ret = 1;
		else ret = 0;
	}
	return ret;
}
