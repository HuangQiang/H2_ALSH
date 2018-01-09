#ifndef __PRI_QUEUE_H
#define __PRI_QUEUE_H


// -----------------------------------------------------------------------------
//  MinK_List: a structure which maintains the smallest k values (of type float)
//  and associated object id (of type int).
//
//  This structure is used for ANN search
//
//  It is currently implemented using an array with k items. Items are stored 
//  in descending order by key, and its insertion is made through standard 
//  insertion mode. (It is quite inefficient if k is large, but current 
//  applications are able to handle since the value of k is small.)
//
//  Note that the priority queue contains k + 1 entries, while the last entry
//  is used as a simple place holder and is otherwise ignored.
// -----------------------------------------------------------------------------
class MinK_List {
public:
	MinK_List(int max);				// constructor (given max size)
	~MinK_List();					// destructor

	void reset();					// make exist list empty

	float min_key();				// return minimum key
	float max_key();				// return maximum key

	float ith_smallest_key(int i);	// return ith key
	int   ith_smallest_id(int i);	// return ith id

	bool isFull();					// is full?

	int size();						// number of elements

	void insert(					// insert item (inline for speed)
		float key,						// key of item
		int id);						// id of item

private:
	struct MinK_Node {				// node in MinK_List
		float key_;						// key value
		int id_;						// object id
	};

	int k_;							// max number of keys
	int num_;						// number of key current active
	MinK_Node* mk_;					// the list itself
};


// -----------------------------------------------------------------------------
//  MaxK_List: An MaxK_List structure is one which maintains the largest k 
//  values (of type float) and associated object id (of type int).
//
//  This structure is used for MIP search
//
//  It is currently implemented using an array with k items. Items are stored
//  in descending order by key, and its insertion is made through standard 
//  insertion mode. (It is quite inefficient if k is large, but current 
//  applications are able to handle since the value of k is small.)
//
//  Note that the priority queue contains k + 1 entries, while the last entry
//  is used as a simple place holder and is otherwise ignored.
// -----------------------------------------------------------------------------
class MaxK_List {
public:
	MaxK_List(int max);				// constructor (given max size)
	~MaxK_List();					// destructor

	void reset();					// make exist list empty

	float min_key();				// return minimum key
	float max_key();				// return maximum key

	float ith_largest_key(int i);	// return ith key
	int ith_largest_id(int i);		// return ith id

	int size();						// number of elements

	void insert(					// insert item (inline for speed)
		float key,						// key of item
		int id);						// id of item

private:
	struct MaxK_Node {				// node in MinK_List
		float key_;						// key value
		int id_;						// object id
	};

	int k_;							// max numner of keys
	int num_;						// number of key current active
	MaxK_Node* mk_;					// the list itself
};

#endif
