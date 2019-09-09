#ifndef __L2_ALSH_H
#define __L2_ALSH_H

class QALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  L2_ALSH is used to solve the problem of c-Approximate Maximum Inner Product 
//  (c-AMIP) search.
//
//  the idea was introduced by Anshumali Shrivastava and Ping Li in their paper 
//  "Asymmetric LSH (ALSH) for sublinear time Maximum Inner Product Search 
//  (MIPS)", In Advances in Neural Information Processing Systems (NIPS), pages
//  2321–2329, 2014.
//  
//  notice that in order to make a fair comparison with H2-ALSH, we apply 
//  QALSH for ANN search after converting MIP search to NN search by the 
//  L2_ALSH transformation.
// -----------------------------------------------------------------------------
class L2_ALSH {
public:
	L2_ALSH(						// constructor
		int   n,						// number of data objects
		int   d,						// dimension of data objects
		int   m,						// additional dimension of data
		float U,						// scale factor for data
		float ratio,					// approximation ratio
		FILE  *fp,						// output file pointer
		const float** data);			// data objects

	// -------------------------------------------------------------------------
	~L2_ALSH();						// destructor

	// -------------------------------------------------------------------------
	int kmip(						// c-k-AMIP search
		int   top_k,					// top-k value
		const float *query,				// input query
		MaxK_List *list);				// top-k MIP results (return) 

protected:
	int   n_pts_;					// number of data points
	int   dim_;						// dimension of data
	int   m_;						// additional dimension of data
	float U_;						// scale factor
	float appr_ratio_;				// approximation ratio for ANN search
	const float **data_;			// data objects 

	float M_;						// max norm of data
	int   l2_alsh_dim_;				// dimension of l2_alsh data (dim_ + m_)
	float **l2_alsh_data_;			// l2_alsh data
	QALSH *lsh_;					// qalsh

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading
};

#endif // __L2_ALSH_H
