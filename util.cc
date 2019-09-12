#include "headers.h"


timeval g_start_time;
timeval g_end_time;

float   g_runtime  = -1.0f;
float   g_ratio    = -1.0f;
float   g_recall   = -1.0f;

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	Result *item1 = (Result*) e1;
	Result *item2 = (Result*) e2;

	if (item1->key_ < item2->key_) {
		ret = -1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = 1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	Result *item1 = (Result*) e1;
	Result *item2 = (Result*) e2;

	if (item1->key_ < item2->key_) {
		ret = 1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = -1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}

// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path)							// input path
{
	int len = (int) strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';
									// check whether the directory exists
			int ret = access(path, F_OK);
			if (ret != 0) {			// create the directory
				ret = mkdir(path, 0755);
				if (ret != 0) {
					printf("Could not create directory %s\n", path);
				}
			}
			path[i + 1] = ch;
		}
	}
}

// -----------------------------------------------------------------------------
int read_data(						// read data/query set from disk
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data)						// data/query objects (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int i   = 0;
	int tmp = -1;
	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &tmp);
		for (int j = 0; j < d; ++j) {
			fscanf(fp, " %f", &data[i][j]);
		}
		fscanf(fp, "\n");

		++i;
	}
	assert(feof(fp) && i == n);
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R)							// ground truth results (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int tmp1 = -1;
	int tmp2 = -1;
	fscanf(fp, "%d %d\n", &tmp1, &tmp2);
	assert(tmp1 == qn && tmp2 == MAXK);

	for (int i = 0; i < qn; ++i) {
		for (int j = 0; j < MAXK; ++j) {
			fscanf(fp, "%d %f ", &R[i][j].id_, &R[i][j].key_);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float ret = 0.0f;
	for (int i = 0; i < dim; ++i) {
		ret += (p1[i] * p2[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc L2 square distance
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float diff = 0.0f;
	float ret  = 0.0f;
	for (int i = 0; i < dim; ++i) {
		diff = p1[i] - p2[i];
		ret += diff * diff;
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L2 distance
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	return sqrt(calc_l2_sqr(dim, p1, p2));
}

// -----------------------------------------------------------------------------
float calc_l1_dist(					// calc L1 distance
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float ret = 0.0f;
	for (int i = 0; i < dim; ++i) {
		ret += fabs(p1[i] - p2[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list)					// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && R[last].key_ - list->ith_key(i) > FLOATZERO) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
int get_hits(						// get the number of hits between two ID list
	int   k,							// top-k value
	int   t,							// top-t value
	const Result *R,					// ground truth results 
	MaxK_List *list)					// results returned by algorithms
{
	int i = k - 1;
	int last = t - 1;
	while (i >= 0 && R[last].key_ - list->ith_key(i) > FLOATZERO) {
		i--;
	}
	return min(t, i + 1);
}

// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth MIP results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	MaxK_List *list = new MaxK_List(MAXK);
	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {	
			float ip = calc_inner_product(d, data[j], query[i]);
			list->insert(ip, j + 1);
		}

		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float truth_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int norm_distribution(				// analyse norm distribution of data
	int   n,							// number of data points
	int   d,							// dimension of space
	const float **data,					// data set
	const char  *out_path) 				// output path
{
	// -------------------------------------------------------------------------
	//  calc norm for all data objects
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	vector<float> norm(n, 0.0f);
	float max_norm = MINREAL;

	for (int i = 0; i < n; ++i) {
		norm[i] = sqrt(calc_inner_product(d, data[i], data[i]));
		if (norm[i] > max_norm) max_norm = norm[i];
	}

	// -------------------------------------------------------------------------
	//  get the percentage of frequency of norm
	// -------------------------------------------------------------------------
	int m = 25;
	float interval = max_norm / m;
	printf("m = %d, max_norm = %f, interval = %f\n", m, max_norm, interval);

	vector<int> freq(m, 0);
	for (int i = 0; i < n; ++i) {
		int id = (int) ceil(norm[i] / interval) - 1;
		if (id < 0) id = 0;
		if (id >= m) id = m - 1;
		freq[id]++;
	}

	// -------------------------------------------------------------------------
	//  write norm distribution
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%snorm_distribution.out", out_path);

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	float num = 0.5f / m; 
	float step = 1.0f / m;
	for (int i = 0; i < m; ++i) {
		fprintf(fp, "%.1f\t%f\n", (num + step * i) * 100.0, freq[i] * 100.0 / n);
	}
	fprintf(fp, "\n");
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
		g_start_time.tv_usec) / 1000000.0f;
	printf("Norm distribution: %.6f Seconds\n\n", runtime);

	return 0;
}