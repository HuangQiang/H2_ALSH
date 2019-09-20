#ifndef __XBOX_H
#define __XBOX_H

class QALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  XBox is used to solve the problem of c-Approximate Maximum Inner Product 
//  (c-AMIP) search.
//
//  the idea was introduced by Yoram Bachrach, Yehuda Finkelstein, Ran 
//  Gilad-Bachrach, Liran Katzir, Noam Koenigstein, Nir Nice, and Ulrich Paquet 
//  in their paper "Speeding up the xbox recommender system using a euclidean 
//  transformation for inner-product spaces", In Proceedings of the 8th ACM 
//  Conference on Recommender systems, pages 257–264, 2014.
//
//  notice that in order to make a fair comparison with H2-ALSH, we apply 
//  QALSH for ANN search after converting MIP search to NN search by the 
//  XBox transformation.
// -----------------------------------------------------------------------------
class XBox {
public:
	XBox(							// default constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		float ratio,					// approximation ratio
		FILE  *fp,						// output file pointer
		const float** data);			// original data objects

	// -------------------------------------------------------------------------
	~XBox();						// destructor

	// -------------------------------------------------------------------------
	int kmip(						// c-k-AMIP search
		int   top_k,					// top-k value
		bool  used_new_transform,		// used new transformation
		const float *query,				// input query
		MaxK_List *list);				// top-k MIP results

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	float nn_ratio_;				// approximation ratio for ANN search
	const float **data_;			// data objects

	float M_;						// max norm of data objects
	float **xbox_data_;				// xbox data
	QALSH *lsh_;					// qalsh

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading
};

#endif // __XBOX_H
