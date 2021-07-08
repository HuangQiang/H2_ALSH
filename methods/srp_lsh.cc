#include "srp_lsh.h"

namespace mips {

// -----------------------------------------------------------------------------
SRP_LSH::SRP_LSH(					// constructor
	int   n,							// cardinality of dataset
	int   d,							// dimensionality of dataset
	int   K)							// number of hash tables
	: n_(n), d_(d), K_(K)
{
	m_ = (int) ceil(K / 64.0f);

	// -------------------------------------------------------------------------
	//  generate random projection vectors
	// -------------------------------------------------------------------------
	proj_ = new float*[K];
	for (int i = 0; i < K; ++i) {
		proj_[i] = new float[d];
		for (int j = 0; j < d; ++j) {
			proj_[i][j] = gaussian(0.0f, 1.0f);
		}
	}

	// -------------------------------------------------------------------------
	//  initialize lookup table for all uint16_t values
	// -------------------------------------------------------------------------
	int size = 1 << 16;
	table16_ = new uint32_t[size];
	for (int i = 0; i < size; ++i) {
		table16_[i] = bit_count(i);
	}

	// -------------------------------------------------------------------------
	//  allocate space for hash_key
	// -------------------------------------------------------------------------
	hash_key_ = new uint64_t*[n];
	for (int i = 0; i < n; ++i) hash_key_[i] = new uint64_t[m_];
}

// -----------------------------------------------------------------------------
SRP_LSH::~SRP_LSH()					// destructor
{
	for (int i = 0; i < K_; ++i) {
		delete[] proj_[i]; proj_[i] = NULL;
	}
	delete[] proj_;	proj_ = NULL; 

	delete[] table16_; table16_ = NULL; 
	
	for (int i = 0; i < n_; ++i) {
		delete[] hash_key_[i]; hash_key_[i] = NULL;
	}
	delete[] hash_key_; hash_key_ = NULL;
}

// -----------------------------------------------------------------------------
inline uint32_t SRP_LSH::bit_count(	// count the number of 1 bits of x
	uint32_t x) 						// input uint32_t value
{
    uint32_t num = x - ((x >> 1) & 033333333333) - ((x >> 2) & 011111111111);
    return ((num + (num >> 3)) & 030707070707) % 63;
}

// -----------------------------------------------------------------------------
bool SRP_LSH::calc_hash_code( 		// calc hash code after random projection
	int   id,							// projection vector id
	const float *data)					// input data
{
	return calc_inner_product(d_, proj_[id], data) >= 0 ? true : false;
}

// -----------------------------------------------------------------------------
void SRP_LSH::compress_hash_code( 	// compress hash code with 64 bits
	const bool *hash_code,				// input hash code
	uint64_t* hash_key)					// hash key (return)
{
	int shift = 0;
	for (int i = 0; i < m_; ++i) {
		int size = (i == m_-1 && K_%64 != 0) ? (K_ % 64) : 64;
		uint64_t val = 0;
		for (int j = 0; j < size; ++j) {
			int idx = j + shift;
			if (hash_code[idx]) val |= ((uint64_t) 1 << (63-j));
		}
		hash_key[i] = val;
		shift += size;
	}
}

// -----------------------------------------------------------------------------
void SRP_LSH::display()				// display parameters
{
	printf("Parameters of SRP_LSH:\n");
	printf("    n = %d\n", n_);
	printf("    d = %d\n", d_);
	printf("    K = %d\n", K_);
	printf("    m = %d\n", m_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int SRP_LSH::kmc(					// c-k-AMC search
	int   top_k,						// top-k value
	const float *query,					// input query
	std::vector<int> &cand) 			// MCS candidates  (return)
{
	// -------------------------------------------------------------------------
	//  calculate the hash key (compressed hash code) of query
	// -------------------------------------------------------------------------
	bool *hash_code_q = new bool[K_];
	for (int i = 0; i < K_; ++i) {
		hash_code_q[i] = calc_hash_code(i, query);
	}
	uint64_t *hash_key_q = new uint64_t[m_];
	compress_hash_code((const bool*) hash_code_q, hash_key_q);

	// -------------------------------------------------------------------------
	//  find the candidates with largest matched values
	// -------------------------------------------------------------------------
	MaxK_List *list = new MaxK_List(CANDIDATES + top_k - 1);
	int total_bits = 64 * m_;
	for (int i = 0; i < n_; ++i) {
		uint32_t match = 0;
		for (int j = 0; j < m_; ++j) {
			match += table_lookup(hash_key_[i][j] ^ hash_key_q[j]);
		}
		list->insert((float) (total_bits - match), i);
	}

	int size = list->size();
	for (int i = 0; i < size; ++i) {
		cand.push_back(list->ith_id(i));
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] hash_code_q; hash_code_q = NULL;
	delete[] hash_key_q;  hash_key_q  = NULL;
	delete list; list = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
inline uint32_t SRP_LSH::table_lookup( // table lookup the match value
	uint64_t x)							// input uint64_t value
{
	return table16_[x & 0xffff] + table16_[(x>>16) & 0xffff] + 
		table16_[(x>>32) & 0xffff] + table16_[(x>>48) & 0xffff];
}

} // end namespace mips
