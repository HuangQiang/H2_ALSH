#include "headers.h"

// -----------------------------------------------------------------------------
void usage() 						// display the usage of this package
{
	printf("\n"
		"-------------------------------------------------------------------\n"
		" Usage of the package for internal approximate MIP search\n"
		"-------------------------------------------------------------------\n"
		"    -alg  {integer}  options of algorithms (0 - 13)\n"
		"    -n    {integer}  cardinality of the dataset\n"
		"    -d    {integer}  dimensionality of the dataset\n"
		"    -qn   {integer}  number of queries\n"
		"    -K    {integer}  number of hash tables for Sign_ALSH and Simple_LSH\n"
		"    -m    {integer}  additional dimension for L2_ALSH, L2_ALSH2, and Sign_ALSH\n"
		"    -U    {real}     a real value (0, 1] for L2_ALSH, L2_ALSH2, and Sign_ALSH\n"
		"    -c    {real}     approximation ratio of nn (c > 1)\n"
		"    -C    {real}     approximation ratio of mip (0 < C < 1) for H2_ALSH\n"
		"    -ds   {string}   address of the dataset\n"
		"    -qs   {string}   address of the query set\n"
		"    -ts   {string}   address of the ground truth set\n"
		"    -of   {string}   output folder to store output results\n"
		"\n"
		"-------------------------------------------------------------------\n"
		" The options of algorithms are:\n"
		"-------------------------------------------------------------------\n"
		"    0 - Generating ground-truth results\n"
		"        Parameters: -alg 0 -n -qn -d -ds -qs -ts\n"
		"\n"
		"    1 - Running AMIP search by L2_ALSH\n"
		"        Parameters: -alg 1 -n -qn -d -m -U -c -ds -qs -ts -of\n"
		"\n"
		"    2 - Running AMIP search by L2_ALSH2\n"
		"        Parameters: -alg 2 -n -qn -d -m -U -c -ds -qs -ts -of\n"
		"\n"
		"    3 - Running AMIP search by XBox and XBox2\n"
		"        Parameters: -alg 3 -n -qn -d -c -ds -qs -ts -of\n"
		"\n"
		"    4 - Running AMIP search by H2_ALSH\n"
		"        Parameters: -alg 4 -n -qn -d -c -C -ds -qs -ts -of\n"
		"\n"
		"    5 - Running AMIP search by Sign_ALSH\n"
		"        Parameters: -alg 5 -n -qn -d -K -m -U -c -ds -qs -ts -of\n"
		"\n"
		"    6 - Running AMIP search by Simple_LSH\n"
		"        Parameters: -alg 6 -n -qn -d -K -c -ds -qs -ts -of\n"
		"\n"
		"    7 - Running MIP search by Linear_Scan\n"
		"        Parameters: -alg 6 -n -qn -d -B -qs -ts -df -of\n"
		"\n"
		"    8 - Running Precision-Recall Curve of AMIP search by L2_ALSH\n"
		"        Parameters: -alg 1 -n -qn -d -m -U -c -ds -qs -ts -of\n"
		"\n"
		"    9 - Running Precision-Recall Curve of AMIP search by L2_ALSH2\n"
		"        Parameters: -alg 2 -n -qn -d -m -U -c -ds -qs -ts -of\n"
		"\n"
		"    10 - Running Precision-Recall Curve of AMIP search by XBox and XBox2\n"
		"        Parameters: -alg 3 -n -qn -d -c -ds -qs -ts -of\n"
		"\n"
		"    11 - Running Precision-Recall Curve of AMIP search by H2_ALSH\n"
		"        Parameters: -alg 4 -n -qn -d -c -C -ds -qs -ts -of\n"
		"\n"
		"    12 - Running Precision-Recall Curve of AMIP search by Sign_ALSH\n"
		"        Parameters: -alg 5 -n -qn -d -K -m -U -c -ds -qs -ts -of\n"
		"\n"
		"    13 - Running Precision-Recall Curve of AMIP search by Simple_LSH\n"
		"        Parameters: -alg 6 -n -qn -d -K -c -ds -qs -ts -of\n"
		"\n"
		"-------------------------------------------------------------------\n"
		" Authors:\n"
		"-------------------------------------------------------------------\n"
		"    Qiang Huang (huangq2011@gmail.com)\n"
		"    Guihong Ma  (maguihong@vip.qq.com)\n\n\n");
}

// -----------------------------------------------------------------------------
int main(int nargs, char **args)
{
	srand((unsigned) time(NULL));	// set the random seed
	usage();

	int alg = -1;
	int n  = -1;					// cardinality
	int qn = -1;					// query number
	int d  = -1;					// dimensionality
	int K = -1;						// number of hash tables

	int   m = -1;
	float U = -1.0f;
	float nn_ratio = -1.0f;
	float mip_ratio = -1.0f;

	char data_set[200];				// address of data set
	char query_set[200];			// address of query set
	char truth_set[200];			// address of ground truth file
	char output_folder[200];		// output folder

	bool failed = false;
	int  cnt = 1;
	
	while (cnt < nargs && !failed) {
		if (strcmp(args[cnt], "-alg") == 0) {
			alg = atoi(args[++cnt]);
			printf("alg = %d\n", alg);
			if (alg < 0 || alg > 13) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-n") == 0) {
			n = atoi(args[++cnt]);
			printf("n = %d\n", n);
			if (n <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-d") == 0) {
			d = atoi(args[++cnt]);
			printf("d = %d\n", d);
			if (d <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-qn") == 0) {
			qn = atoi(args[++cnt]);
			printf("qn = %d\n", qn);
			if (qn <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-K") == 0) {
			K = atoi(args[++cnt]);
			printf("K = %d\n", K);
			if (K <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-m") == 0) {
			m = atoi(args[++cnt]);
			printf("m = %d\n", m);
			if (m <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-U") == 0) {
			U = (float) atof(args[++cnt]);
			printf("U = %.2f\n", U);
			if (U <= 0.0f || U > 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-c") == 0) {
			nn_ratio = (float) atof(args[++cnt]);
			printf("c = %.2f\n", nn_ratio);
			if (nn_ratio <= 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-C") == 0) {
			mip_ratio = (float) atof(args[++cnt]);
			printf("C = %.2f\n", mip_ratio);
			if (mip_ratio <= 0.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-ds") == 0) {
			strncpy(data_set, args[++cnt], sizeof(data_set));
			printf("dataset = %s\n", data_set);
		}
		else if (strcmp(args[cnt], "-qs") == 0) {
			strncpy(query_set, args[++cnt], sizeof(query_set));
			printf("query set = %s\n", query_set);
		}
		else if (strcmp(args[cnt], "-ts") == 0) {
			strncpy(truth_set, args[++cnt], sizeof(truth_set));
			printf("truth set = %s\n", truth_set);
		}
		else if (strcmp(args[cnt], "-of") == 0) {
			strncpy(output_folder, args[++cnt], sizeof(output_folder));
			printf("output folder = %s\n", output_folder);
			check_path(output_folder);
		}
		else {
			failed = true;
			error("Parameters error!\n", false);
			usage();
			break;
		}
		cnt++;
	}

	// -------------------------------------------------------------------------
	//  Read data set and query set
	// -------------------------------------------------------------------------
	clock_t startTime = (clock_t) -1;
	clock_t endTime   = (clock_t) -1;

	startTime = clock();
	float** data = new float*[n];
	for (int i = 0; i < n; i++) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		error("Reading dataset error!\n", true);
		return 1;
	}

	float** query = new float*[qn];
	for (int i = 0; i < qn; i++) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		error("Reading query set error!\n", true);
		return 1;
	}
	endTime = clock();
	printf("Read Dataset and Query Set: %.6f seconds\n\n",
		((float) endTime - startTime) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  Methods
	// -------------------------------------------------------------------------
	switch (alg) {
	case 0:
		ground_truth(n, qn, d, data, query, truth_set);
		break;
	case 1:
		l2_alsh(n, qn, d, m, U, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 2:
		l2_alsh2(n, qn, d, m, U, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 3:
		xbox(n, qn, d, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 4:
		h2_alsh(n, qn, d, nn_ratio, mip_ratio, data, query, truth_set, output_folder);
		break;
	case 5:
		sign_alsh(n, qn, d, K, m, U, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 6:
		simple_lsh(n, qn, d, K, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 7:
		linear_scan(n, qn, d, data, query, truth_set, output_folder);
		break;
	case 8:
		l2_alsh_precision_recall(n, qn, d, m, U, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 9:
		l2_alsh2_precision_recall(n, qn, d, m, U, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 10:
		xbox_precision_recall(n, qn, d, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 11:
		h2_alsh_precision_recall(n, qn, d, nn_ratio, mip_ratio, data, query, truth_set, output_folder);
		break;
	case 12:
		sign_alsh_precision_recall(n, qn, d, K, m, U, nn_ratio, data, query, truth_set, output_folder);
		break;
	case 13:
		simple_lsh_precision_recall(n, qn, d, K, nn_ratio, data, query, truth_set, output_folder);
		break;
	default:
		error("Parameters error!\n", true);
		usage();
		break;
	}
	
	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n; i++) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data  = NULL;

	for (int i = 0; i < qn; i++) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;

	return 0;
}


