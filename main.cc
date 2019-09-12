#include "headers.h"

// -----------------------------------------------------------------------------
void usage() 						// display the usage of this package
{
	printf("\n"
		"-------------------------------------------------------------------\n"
		" Usage of the package for c-Approximate MIP (c-AMIP) search\n"
		"-------------------------------------------------------------------\n"
		"    -alg  {integer}  options of algorithms (0 - 11)\n"
		"    -n    {integer}  cardinality of the dataset\n"
		"    -d    {integer}  dimensionality of the dataset\n"
		"    -qn   {integer}  number of queries\n"
		"    -K    {integer}  #hash tables for Sign_ALSH and Simple_LSH\n"
		"    -m    {integer}  extra dim for L2_ALSH, L2_ALSH2, Sign_ALSH\n"
		"    -U    {real}     range (0,1] for L2_ALSH, L2_ALSH2, Sign_ALSH\n"
		"    -c0   {real}     approximation ratio of ANN search (c0 > 1)\n"
		"    -c    {real}     approximation ratio of AMIP search (0 < c < 1)\n"
		"    -ds   {string}   address of the data  set\n"
		"    -qs   {string}   address of the query set\n"
		"    -ts   {string}   address of the truth set\n"
		"    -op   {string}   output path\n"
		"\n"
		"-------------------------------------------------------------------\n"
		" The options of algorithms are:\n"
		"-------------------------------------------------------------------\n"
		"    0  - Ground-Truth\n"
		"         Parameters: -alg 0 -n -qn -d -ds -qs -ts\n"
		"\n"
		"    1  - AMIP search by H2_ALSH\n"
		"         Parameters: -alg 1 -n -qn -d -c0 -c -ds -qs -ts -op\n"
		"\n"
		"    2  - AMIP search by L2_ALSH\n"
		"         Parameters: -alg 2 -n -qn -d -m -U -c0 -ds -qs -ts -op\n"
		"\n"
		"    3  - AMIP search by L2_ALSH2\n"
		"         Parameters: -alg 3 -n -qn -d -m -U -c0 -ds -qs -ts -op\n"
		"\n"
		"    4  - AMIP search by XBOX and H2-ALSH-\n"
		"         Parameters: -alg 4 -n -qn -d -c0 -ds -qs -ts -op\n"
		"\n"
		"    5  - AMIP search by Sign_ALSH\n"
		"         Parameters: -alg 5 -n -qn -d -K -m -U -c0 -ds -qs -ts -op\n"
		"\n"
		"    6  - AMIP search by Simple_LSH\n"
		"         Parameters: -alg 6 -n -qn -d -K -c0 -ds -qs -ts -op\n"
		"\n"
		"    7  - MIP search by Linear_Scan\n"
		"         Parameters: -alg 7 -n -qn -d -B -qs -ts -df -op\n"
		"\n"
		"    8  - Precision-Recall Curve of AMIP search by H2_ALSH\n"
		"         Parameters: -alg 8 -n -qn -d -c0 -c -ds -qs -ts -op\n"
		"\n"
		"    9  - Precision-Recall Curve of AMIP search by Sign_ALSH\n"
		"         Parameters: -alg 9 -n -qn -d -K -m -U -c0 -ds -qs -ts -op\n"
		"\n"
		"    10 - Precision-Recall Curve of AMIP search by Simple_LSH\n"
		"         Parameters: -alg 10 -n -qn -d -K -c0 -ds -qs -ts -op\n"
		"\n"
		"    11 - Norm Distributiuon\n"
		"         Parameters: -alg 11 -n -d -ds -op\n"
		"\n"
		"-------------------------------------------------------------------\n"
		" Authors: Qiang Huang (huangq2011@gmail.com)                       \n"
		"          Guihong Ma  (maguihong@vip.qq.com)                       \n"
		"-------------------------------------------------------------------\n"
		"\n\n\n");
}

// -----------------------------------------------------------------------------
int main(int nargs, char **args)
{
	srand((unsigned) time(NULL));	// set the random seed
	//usage();

	int    alg       = -1;			// which algorithm?
	int    n         = -1;			// cardinality
	int    qn        = -1;			// query number
	int    d         = -1;			// dimensionality
	int    K         = -1;			// #tables for sign-alsh and simple-lsh
	int    m         = -1;			// param for l2-alsh, l2-alsh2, sign-alsh
	float  U         = -1.0f;		// param for l2-alsh, l2-alsh2, sign-alsh
	float  nn_ratio  = -1.0f;		// approximation ratio of ANN search
	float  mip_ratio = -1.0f;		// approximation ratio of AMIP search
	
	char   data_set[200];			// address of data set
	char   query_set[200];			// address of query set
	char   truth_set[200];			// address of ground truth file
	char   out_path[200];			// output path

	float  **pre     = NULL;		// precision array
	float  **recall  = NULL;		// recall array
	float  **data    = NULL;		// data set
	float  **query   = NULL;		// query set
	Result **R       = NULL;		// truth set
	bool   failed    = false;
	int    cnt       = 1;
	
	while (cnt < nargs && !failed) {
		if (strcmp(args[cnt], "-alg") == 0) {
			alg = atoi(args[++cnt]);
			printf("alg       = %d\n", alg);
			if (alg < 0 || alg > 11) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-n") == 0) {
			n = atoi(args[++cnt]);
			printf("n         = %d\n", n);
			if (n <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-d") == 0) {
			d = atoi(args[++cnt]);
			printf("d         = %d\n", d);
			if (d <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-qn") == 0) {
			qn = atoi(args[++cnt]);
			printf("qn        = %d\n", qn);
			if (qn <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-K") == 0) {
			K = atoi(args[++cnt]);
			printf("K         = %d\n", K);
			if (K <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-m") == 0) {
			m = atoi(args[++cnt]);
			printf("m         = %d\n", m);
			if (m <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-U") == 0) {
			U = (float) atof(args[++cnt]);
			printf("U         = %.2f\n", U);
			if (U <= 0.0f || U > 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-c0") == 0) {
			nn_ratio = (float) atof(args[++cnt]);
			printf("c0        = %.2f\n", nn_ratio);
			if (nn_ratio <= 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-c") == 0) {
			mip_ratio = (float) atof(args[++cnt]);
			printf("c         = %.2f\n", mip_ratio);
			if (mip_ratio <= 0.0f || mip_ratio >= 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-ds") == 0) {
			strncpy(data_set, args[++cnt], sizeof(data_set));
			printf("data_set  = %s\n", data_set);
		}
		else if (strcmp(args[cnt], "-qs") == 0) {
			strncpy(query_set, args[++cnt], sizeof(query_set));
			printf("query_set = %s\n", query_set);
		}
		else if (strcmp(args[cnt], "-ts") == 0) {
			strncpy(truth_set, args[++cnt], sizeof(truth_set));
			printf("truth_set = %s\n", truth_set);
		}
		else if (strcmp(args[cnt], "-op") == 0) {
			strncpy(out_path, args[++cnt], sizeof(out_path));
			printf("out_path  = %s\n", out_path);

			int len = (int) strlen(out_path);
			if (out_path[len - 1] != '/') {
				out_path[len] = '/';
				out_path[len + 1] = '\0';
			}
			create_dir(out_path);
		}
		else {
			failed = true;
			usage();
			break;
		}
		cnt++;
	}
	printf("\n");

	// -------------------------------------------------------------------------
	//  read data set, query set, and ground truth file
	// -------------------------------------------------------------------------
	data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		return 1;
	}

	if (alg >= 0 && alg <= 10) {
        query = new float*[qn];
        for (int i = 0; i < qn; ++i) query[i] = new float[d];
        if (read_data(qn, d, query_set, query) == 1) {
            return 1;
        }
    }

	if (alg >= 1 && alg <= 10) {
		R = new Result*[qn];
		for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
		if (read_ground_truth(qn, truth_set, R) == 1) {
			return 1;
		}
	}

	if (alg >= 8 && alg <= 10) {
		pre    = new float*[MAX_ROUND];
		recall = new float*[MAX_ROUND];

		for (int round = 0; round < MAX_ROUND; ++round) {
			pre[round]    = new float[MAX_T];
			recall[round] = new float[MAX_T];

			for (int t = 0; t < MAX_T; ++t) {
				pre[round][t]    = 0.0f;
				recall[round][t] = 0.0f;
			}
		}
	}

	// -------------------------------------------------------------------------
	//  methods
	// -------------------------------------------------------------------------
	switch (alg) {
	case 0:
		ground_truth(n, qn, d, (const float **) data, 
            (const float **) query, truth_set);
		break;
	case 1:
		h2_alsh(n, qn, d, nn_ratio, mip_ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 2:
		l2_alsh(n, qn, d, m, U, nn_ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 3:
		l2_alsh2(n, qn, d, m, U, nn_ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 4:
		xbox(n, qn, d, nn_ratio, (const float **) data, (const float **) query, 
			(const Result **) R, out_path);
		break;
	case 5:
		sign_alsh(n, qn, d, K, m, U, nn_ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 6:
		simple_lsh(n, qn, d, K, nn_ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 7:
		linear_scan(n, qn, d, (const float **) data, (const float **) query,
			(const Result **) R, out_path);
		break;
	case 8:
		h2_alsh_precision_recall(n, qn, d, nn_ratio, mip_ratio, pre, recall, 
			(const float **) data, (const float **) query, 
			(const Result **) R, out_path);
		break;
	case 9:
		sign_alsh_precision_recall(n, qn, d, K, m, U, nn_ratio, pre, recall, 
			(const float **) data, (const float **) query, 
			(const Result **) R, out_path);
		break;
	case 10:
		simple_lsh_precision_recall(n, qn, d, K, nn_ratio, pre, recall, 
			(const float **) data, (const float **) query, 
			(const Result **) R, out_path);
		break;
	case 11:
		norm_distribution(n, d, (const float **) data, out_path);
		break;
	default:
		printf("Parameters error!\n");
		usage();
		break;
	}
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data  = NULL;

    if (alg >= 0 && alg <= 10) {
        for (int i = 0; i < qn; ++i) {
            delete[] query[i]; query[i] = NULL;
        }
        delete[] query; query = NULL;
    }

	if (alg >= 1 && alg <= 10) {
		delete[] R; R = NULL;
	}

	if (alg >= 8 && alg <= 10) {
		for (int i = 0; i < MAX_ROUND; ++i) {
			delete[] pre[i];	pre[i] = NULL;
			delete[] recall[i];	recall[i] = NULL;
		}
		delete[] pre;	 pre = NULL;
		delete[] recall; recall = NULL;
	}
	
	return 0;
}