#include "headers.h"

// -----------------------------------------------------------------------------
//  Global functions
// -----------------------------------------------------------------------------
void error(							// an error message
	char* msg,							// an message
	bool is_exit)						// whether exit the program
{
	printf(msg);
	if (is_exit) exit(1);
}

// -----------------------------------------------------------------------------
int read_data(						// read data from disk
	int n,								// number of data points
	int d,								// dimensionality
	char* fname,						// address of data
	float** data)						// data (return)
{
	int i = 0;
	int j = 0;
	FILE* fp = NULL;

	fp = fopen(fname, "r");			// open the file
	if (!fp) {
		error("Could not open the file.\n", true);
		return 1;
	}

	i = 0;
	while (!feof(fp) && i < n) {	// read data
		fscanf(fp, "%d", &j);
		for (j = 0; j < d; j++) {
			fscanf(fp, " %f", &data[i][j]);
		}
		fscanf(fp, "\n");

		i++;
	}
	if (!feof(fp) && i == n) {		// check the size of dataset
		error("The size of dataset is larger than you set\n", true);
		return 1;
	}
	else if (feof(fp) && i < n) {
		printf("Set the size of dataset to be %d. ", i);
		error("And try again\n", true);
		return 1;
	}
	fclose(fp);						// close data file
	return 0;
}

// -----------------------------------------------------------------------------
int check_path(						// check input path whether exist
	char* path)							// input path
{
	int len = (int) strlen(path);	// ensure the path is a folder
	if (path[len - 1] != '/') {
		path[len] = '/';
		path[len + 1] = '\0';
	}

	// -------------------------------------------------------------------------
	//  Check whether the directory exists. If the directory does not exist, we 
	//  create the directory for each folder.
	// -------------------------------------------------------------------------
#ifdef LINUX_
	len = (int) strlen(path);		// create directory under Linux
	for (int i = 0; i < len; i++) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';
									// check whether the directory exists
			int ret = access(path, F_OK);
			if (ret != 0) {			// create the directory
				ret = mkdir(path, 0755);
				if (ret != 0) {
					printf("Could not create directory %s\n", path);
					error("write_data_new_form error\n", true);
					return 1;		// fail to return
				}
			}
			path[i + 1] = ch;
		}
	}
#else
	len = (int) strlen(path);		// create directory under Windows
	for (int i = 0; i < len; i++) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';
									// check whether the directory exists
			int ret = _access(path, 0);
			if (ret != 0) {			// create the directory
				ret = _mkdir(path);
				if (ret != 0) {
					printf("Could not create directory %s\n", path);
					error("write_data_new_form error\n", true);
					return 1;		// fail to return
				}
			}
			path[i + 1] = ch;
		}
	}
#endif
	return 0;
}

// -----------------------------------------------------------------------------
float uniform(						// r.v. from Uniform(min, max)
	float min,							// min value
	float max)							// max value
{
	int num = rand();
	float base = (float) RAND_MAX - 1.0F;
	float frac  = ((float) num) / base;

	return (max - min) * frac + min;
}

// -----------------------------------------------------------------------------
//	Given a mean and a standard deviation, gaussian generates a normally 
//		distributed random number.
//
//	Algorithm:  Polar Method, p.104, Knuth, vol. 2
// -----------------------------------------------------------------------------
float gaussian(						// r.v. from Gaussian(mean, sigma)
	float mean,							// mean value
	float sigma)						// std value
{
	float v1, v2;
	float s;
	float x;

	do {
		v1 = 2.0F * uniform(0.0F, 1.0F) - 1.0F;
		v2 = 2.0F * uniform(0.0F, 1.0F) - 1.0F;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.0F);
	x = v1 * sqrt (-2.0F * log (s) / s);

	x = x * sigma + mean; 			// x is distributed from N(0, 1)
	return x;
}

// -----------------------------------------------------------------------------
float normal_pdf(					// pdf of Guassian(mean, std)
	float x,							// variable
	float u,							// mean
	float sigma)						// standard error
{
	float ret = exp(-(x - u) * (x - u) / (2.0f * sigma * sigma));
	ret /= sigma * sqrt(2.0f * PI);
	return ret;
}

// -----------------------------------------------------------------------------
float normal_cdf(					// cdf of N(0, 1) in range (-inf, x]
	float _x,							// integral border
	float _step)						// step increment
{
	float ret = 0.0;
	for (float i = -10.0; i < _x; i += _step) {
		ret += _step * normal_pdf(i, 0.0f, 1.0f);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float new_cdf(						// cdf of N(0, 1) in range [-x, x]
	float x,							// integral border
	float step)							// step increment
{
	float result = 0.0f;
	for (float i = -x; i <= x; i += step) {
		result += step * normal_pdf(i, 0.0f, 1.0f);
	}
	return result;
}

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim)							// dimension
{
	float ret = 0.0F;
	for (int i = 0; i < dim; i++) {
		ret += (p1[i] * p2[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc L2 square distance
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim)							// dimension
{
	float diff = 0.0F;
	float ret  = 0.0F;
	for (int i = 0; i < dim; i++) {
		diff = p1[i] - p2[i];
		ret += diff * diff;
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L2 distance
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim)							// dimension
{
	float ret = calc_l2_sqr(p1, p2, dim);
	return sqrt(ret);
}

// -----------------------------------------------------------------------------
float calc_l1_dist(					// calc L1 distance
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim)							// dimension
{
	float ret = 0.0F;
	for (int i = 0; i < dim; i++) {
		ret += fabs(p1[i] - p2[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall between two ID list
	int n,								// length of ID list
	MaxK_List* list,					// ID list of k-nn method
	map<int, int>& id_rank)				// ground truth ID list
{
	int count = 0;
	for (int i = 0; i < n; i++) {
		int obj_id = list->ith_largest_id(i);

		if (id_rank.find(obj_id) != id_rank.end() && id_rank[obj_id] <= n) {
			count++;
		}
	}
	float recall = (float)count / (float)n;
	return recall;
}

// -----------------------------------------------------------------------------
int get_hits(						// get the number of hits between two ID list
	int k,								// length of ID list
	int t,								// top-t
	MaxK_List* list,					// ID list of k-nn method
	map<int, int>& id_rank)				// ground truth ID list
{
	int count = 0;
	for (int i = 0; i < k; i++) {
		int obj_id = list->ith_largest_id(i);

		if (id_rank.find(obj_id) != id_rank.end() && id_rank[obj_id] <= t) {
			count++;
		}
	}
	return min(t, count);
}
