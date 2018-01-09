#include "headers.h"


// -----------------------------------------------------------------------------
SRP_LSH::SRP_LSH()					// constructor
{
	K_ = -1;
	dim_ = -1;
	n_pts_ = -1;

	proj_ = NULL;
	hash_code_ = NULL;
}

// -----------------------------------------------------------------------------
SRP_LSH::~SRP_LSH()					// destructor
{
	if (proj_ != NULL) {
		for (int i = 0; i < K_; i++) {
			delete[] proj_[i];	proj_[i] = NULL;
		}
		delete[] proj_;	proj_ = NULL;
	}

	if (hash_code_ != NULL) {
		for (int i = 0; i < n_pts_; i++) {
			delete[] hash_code_[i];	hash_code_[i] = NULL;
		}
		delete[] hash_code_; hash_code_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::init(					// initialize the parameters
	int K,								// number of hash tables
	int d,								// dimensionality of dataset
	int n,								// cardinality of dataset
	float** data)						// input data objects
{
	n_pts_ = n;
	K_ = K;
	dim_ = d;
	data_ = data;

	pre_processing();
}

// -----------------------------------------------------------------------------
void SRP_LSH::pre_processing()		// pre-processing of data
{
	generate();

	hash_code_ = new bool*[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		hash_code_[i] = new bool[K_];
		get_projection_vector(data_[i], hash_code_[i]);
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::generate()			// generate random projection vectors
{
	proj_ = new float*[K_];
	for (int k = 0; k < K_; k++) {
		proj_[k] = new float[dim_];
		for (int d = 0; d < dim_; d++) {
			proj_[k][d] = gaussian(0.0f, 1.0f);
		}
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::get_projection_vector(// get vector after random projection
	float* data,						// input data 
	bool* hash_code)					// hash code of input data (return)
{
	for (int k = 0; k < K_; k++) {
		float sum = 0.0f;
		for (int d = 0; d < dim_; d++) {
			sum += proj_[k][d] * data[d];
		}
		if (sum >= 0) {
			hash_code[k] = true;
		}
		else {
			hash_code[k] = false;
		}
	}
}

// -----------------------------------------------------------------------------
int SRP_LSH::kmc(					// top-k approximate maximum cosine search
	float* query,						// input query
	int top_k,							// top-k value
	MaxK_List* list)					// top-k mip results
{
	bool* mc_query = new bool[K_];
	get_projection_vector(query, mc_query);

	for (int i = 0; i < n_pts_; i++) {
		int match = 0;
		for (int k = 0; k < K_; k++) {
			if (hash_code_[i][k] == mc_query[k]) {
				match++;
			}
		}
		//cout << i << " " << dist << endl;
		list->insert(match, i);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] mc_query; mc_query = NULL;

	return 0;
}