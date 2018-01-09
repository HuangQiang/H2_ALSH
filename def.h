#ifndef __DEF_H
#define __DEF_H

// -----------------------------------------------------------------------------
//  Macros
// -----------------------------------------------------------------------------
#define MIN(a, b)	(((a) < (b)) ? (a) : (b))
#define MAX(a, b)	(((a) > (b)) ? (a) : (b))
#define SQR(x)		((x) * (x))
#define SUM(x, y)	((x) + (y))
#define DIFF(x, y)	((y) - (x))
#define SWAP(x, y)	{int tmp=x; x=y; y=tmp;}

// -----------------------------------------------------------------------------
//  Constants
// -----------------------------------------------------------------------------
const float E  = 2.71828182845904F;	// math constants
const float PI = 3.14159265358979F;
									// max real value
const float MAXREAL = 3.402823466e+38F;
const float MINREAL = -MAXREAL;		// min real value
const float FLOATZERO = 1e-6F;		// accuracy

const int MAXINT = 2147483647;		// max integer value
const int MININT = -MAXINT;			// min integer value

const int MAXK = 10;				// max top-k value
const int STRING_LEN = 200;			// max length of file name
const int SCAN_SIZE = 512;			// half size of one scan
const int CANDIDATES = 100;			// candidate for qalsh
									// max memory, 8 GB
const long long MAXMEMORY = 8589934591;

const int SIZEINT   = (int) sizeof(int);
const int SIZECHAR  = (int) sizeof(char);
const int SIZEFLOAT = (int) sizeof(float);
const int SIZEBOOL  = (int) sizeof(bool);

#endif
