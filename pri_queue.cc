#include "headers.h"


// -----------------------------------------------------------------------------
//  MinK_List: a structure which maintains the smallest k values (of type float)
//  and associated object id (of type int).
// -----------------------------------------------------------------------------
MinK_List::MinK_List(				// constructor (given max size)
	int max)							// max size
{
	num_ = 0;
	k_ = max;
	mk_ = new MinK_Node[max + 1];
}

// -----------------------------------------------------------------------------
MinK_List::~MinK_List() 			// destructor
{
	if (mk_ != NULL) {
		delete[] mk_; mk_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void MinK_List::reset()				// make exist queue empty
{
	num_ = 0;						// setting current active num = 0
}

// -----------------------------------------------------------------------------
float MinK_List::min_key()			// return minimum key
{
	return (num_ > 0 ? mk_[0].key_ : MAXREAL);
}

// -----------------------------------------------------------------------------
float MinK_List::max_key()			// return maximum key
{
	return (num_ >= k_ ? mk_[k_-1].key_ : MAXREAL);
}

// -----------------------------------------------------------------------------
float MinK_List::ith_smallest_key(	// return ith key
	int i)								// ith position
{
	return (i < num_ ? mk_[i].key_ : MAXREAL);
}

// -----------------------------------------------------------------------------
int MinK_List::ith_smallest_id(		// return ith info
	int i)								// ith position
{
	return (i < num_ ? mk_[i].id_ : MININT);
}

// -----------------------------------------------------------------------------
bool MinK_List::isFull()			// is full?
{
	if (num_ >= k_) return true;
	else return false;
}

// -----------------------------------------------------------------------------
int MinK_List::size()				// number of elements
{
	return num_;
}

// -----------------------------------------------------------------------------
void MinK_List::insert(				// insert item (inline for speed)
	float key,							// key of item
	int id)								// id of item
{
	int i = 0;
	for (i = num_; i > 0; i--) {
		if (mk_[i-1].key_ > key) mk_[i] = mk_[i-1];
		else break;
	}
	mk_[i].key_ = key;				// store new item here
	mk_[i].id_ = id;
	if (num_ < k_) num_++;			// increament number of items
}

// -----------------------------------------------------------------------------
//  MaxK_List: An MaxK_List structure is one which maintains the largest k 
//  values (of type float) and associated object id (of type int).
// -----------------------------------------------------------------------------
MaxK_List::MaxK_List(				// constructor (given max size)
	int max)							// max size
{
	num_ = 0;
	k_ = max;
	mk_ = new MaxK_Node[max + 1];
}

// -----------------------------------------------------------------------------
MaxK_List::~MaxK_List() 			// destructor
{
	if (mk_ != NULL) {
		delete[] mk_;
		mk_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void MaxK_List::reset()				// make exist queue empty
{
	num_ = 0;						// setting current active num = 0
}

// -----------------------------------------------------------------------------
float MaxK_List::max_key()			// return minimum key
{
	return (num_ > 0 ? mk_[0].key_ : MINREAL);
}

// -----------------------------------------------------------------------------
float MaxK_List::min_key()			// return maximum key
{
	return (num_ == k_ ? mk_[k_ - 1].key_ : MINREAL);
}

// -----------------------------------------------------------------------------
float MaxK_List::ith_largest_key(	// return ith key
	int i)								// ith position
{
	return (i < num_ ? mk_[i].key_ : MINREAL);
}

// -----------------------------------------------------------------------------
int MaxK_List::ith_largest_id(		// return ith info
	int i)								// ith position
{
	return (i < num_ ? mk_[i].id_ : MININT);
}

// -----------------------------------------------------------------------------
int MaxK_List::size()				// number of elements
{
	return num_;
}

// -----------------------------------------------------------------------------
void MaxK_List::insert(				// insert item (inline for speed)
	float key,							// key of item
	int id)								// id of item
{
	int i = 0;
	for (i = num_; i > 0; i--) {
		if (mk_[i - 1].key_ < key) {
			mk_[i].key_ = mk_[i - 1].key_;
			mk_[i].id_ = mk_[i - 1].id_;
		}
		else {
			break;
		}
	}
	mk_[i].key_ = key;				// store new item here
	mk_[i].id_ = id;
	if (num_ < k_) num_++;			// increament number of items
}


