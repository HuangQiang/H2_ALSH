#ifndef __L2_ALSH2_H
#define __L2_ALSH2_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"

class QALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  L2_ALSH2 is used to solve the problem of c-Approximate Maximum Inner 
//  Product (c-AMIP) search.
//
//  the idea was introduced by Anshumali Shrivastava and Ping Li in their paper 
//  "Asymmetric LSH (ALSH) for sublinear time Maximum Inner Product Search 
//  (MIPS)", In Advances in Neural Information Processing Systems (NIPS), pages
//  2321–2329, 2014.
//
//  notice that in order to make a fair comparison with H2-ALSH, we apply 
//  QALSH for ANN search after converting MIP search to NN search by the 
//  L2_ALSH2 transformation. 
// 
//  in the problem definition section of our KDD 2018 paper, we assume we do
//  NOT know the Euclidean norm of queries before c-AMIP search. However, this  
//  transformation requires to know this information before c-AMIP search. Thus, 
//  we did NOT compare H2-ALSH with it in our KDD paper. Nevertheless, based on 
//  the results over five real datasets (Mnist, Sift, Gist, Netflix, and Yahoo) 
//  we use, H2-ALSH significantly outperforms L2-ALSH2.  
// -----------------------------------------------------------------------------
class L2_ALSH2 {
public:
	L2_ALSH2(						// constructor
		int   n,						// number of data objects
		int   qn,						// number of query objects
		int   d,						// dimension of data objects
		int   m,						// additional dimension of data
		float U,						// scale factor for data
		float nn_ratio,					// approximation ratio for ANN search
		const float **data, 			// input data
		const float **norm_d,			// l2-norm of data objects
		const float **norm_q);			// l2-norm of query objects

	// -------------------------------------------------------------------------
	~L2_ALSH2();					// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kmip(						// c-k-AMIP search
		int   top_k,					// top-k value
		const float *query,				// input query
		const float *norm_q,			// l2-norm of query
		MaxK_List *list);				// top-k MIP results (return) 

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	int   m_;						// additional dimension of data
	float U_;						// scale factor
	float M_;						// max norm of data and query
	const float **data_;			// original data objects
	const float **norm_d_;			// l2-norm of data objects
	const float **norm_q_;			// l2-norm of query objects
	// int   l2_alsh2_dim_;			// dim of l2_alsh2 data (dim_ + 2 * m_)
	// float **l2_alsh2_data_;			// l2_alsh2 data
	QALSH *lsh_;					// qalsh
};

#endif // __L2_ALSH2_H
