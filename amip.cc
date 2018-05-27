#include "headers.h"

// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	MaxK_List *list = new MaxK_List(MAXK);
	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {	
			float ip = calc_inner_product(d, data[j], query[i]);
			list->insert(ip, j + 1);
		}

		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&end_time, NULL);
	float truth_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int h2_alsh(						// mip search via h2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	H2_ALSH *lsh = new H2_ALSH();
	lsh->build(n, d, nn_ratio, mip_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via H2_ALSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%sh2_alsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k c-AMIP of H2_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int l2_alsh(						// mip search via l2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	L2_ALSH *lsh = new L2_ALSH();
	lsh->build(n, d, m, U, nn_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via L2_ALSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%sl2_alsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k c-AMIP of L2_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int l2_alsh2(						// mip search via l2_alsh2
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	L2_ALSH2 *lsh = new L2_ALSH2();
	lsh->build(n, qn, d, m, U, nn_ratio, data, query);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via L2_ALSH2
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%sl2_alsh2.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k c-AMIP of L2_ALSH2: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int xbox(							// mip search via xbox
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	XBox *lsh = new XBox();
	lsh->build(n, d, nn_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via XBox
	// -------------------------------------------------------------------------
	char output_set[200];
	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	sprintf(output_set, "%sxbox.out", output_folder);
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	printf("Top-k c-AMIP of XBox: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, false, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  c-AMIP search via H2-ALSH-
	// -------------------------------------------------------------------------	
	sprintf(output_set, "%sh2_alsh_minus.out", output_folder);
	fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	printf("Top-k c-AMIP of H2-ALSH-: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, true, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int sign_alsh(						// mip search via sign_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	Sign_ALSH *lsh = new Sign_ALSH();
	lsh->build(n, d, K, m, U, nn_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via Sign_ALSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%ssign_alsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k c-AMIP of Sign_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh(						// mip search via simple_lsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	Simple_LSH *lsh = new Simple_LSH();
	lsh->build(n, d, K, nn_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via Simple_LSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%ssimple_lsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k c-AMIP of Simple_LSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear_scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  c-AMIP search via linear scan
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%slinear.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMIPs[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k MIP of Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			for (int j = 0; j < n; ++j) {
				float ip = calc_inner_product(d, data[j], query[i]);
				list->insert(ip, j + 1);
			}
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision recall curve of h2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	H2_ALSH *lsh = new H2_ALSH();
	lsh->build(n, d, nn_ratio, mip_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of H2_ALSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%sh2_alsh_precision_recall.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int tMIPs[] = { 1, 2, 5, 10 };
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxT_round = 4;
	int maxK_round = 16;

	float **pre    = new float*[maxT_round];
	float **recall = new float*[maxT_round];
	
	for (int t_round = 0; t_round < maxT_round; ++t_round) {
		pre[t_round]    = new float[maxK_round];
		recall[t_round] = new float[maxK_round];

		for (int k_round = 0; k_round < maxK_round; ++k_round) {
			pre[t_round][k_round]    = 0;
			recall[t_round][k_round] = 0;
		}
	}

	printf("Top-t c-AMIP of H2_ALSH: \n");
	for (int k_round = 0; k_round < maxK_round; ++k_round) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);

			for (int t_round = 0; t_round < maxT_round; ++t_round) {
				int top_t = tMIPs[t_round];
				int hits = get_hits(top_k, top_t, R[i], list);

				pre[t_round][k_round]    += hits / (float) top_k;
				recall[t_round][k_round] += hits / (float) top_t;
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxT_round; ++t_round) {
		int top_t = tMIPs[t_round];
		printf("Top-%d\t\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxK_round; ++k_round) {
			int top_k = kMIPs[k_round];
			pre[t_round][k_round]    = pre[t_round][k_round]    * 100.0f / qn;
			recall[t_round][k_round] = recall[t_round][k_round] * 100.0f / qn;

			printf("%4d\t\t%.2f\t\t%.2f\n", top_k,
				recall[t_round][k_round], pre[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", top_k,
				recall[t_round][k_round], pre[t_round][k_round]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;
	for (int i = 0; i < maxT_round; ++i) {
		delete[] pre[i];	pre[i] = NULL;
		delete[] recall[i];	recall[i] = NULL;
	}
	delete[] pre;	 pre = NULL;
	delete[] recall; recall = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int sign_alsh_precision_recall(		// precision recall curve of sign_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	Sign_ALSH *lsh = new Sign_ALSH();
	lsh->build(n, d, K, m, U, nn_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Sign-ALSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%ssign_alsh_precision_recall.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int tMIPs[] = { 1, 2, 5, 10 };
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxT_round = 4;
	int maxK_round = 16;

	float **pre    = new float*[maxT_round];
	float **recall = new float*[maxT_round];
	
	for (int t_round = 0; t_round < maxT_round; ++t_round) {
		pre[t_round]    = new float[maxK_round];
		recall[t_round] = new float[maxK_round];

		for (int k_round = 0; k_round < maxK_round; ++k_round) {
			pre[t_round][k_round]    = 0;
			recall[t_round][k_round] = 0;
		}
	}

	printf("Top-t c-AMIP of Sign_ALSH: \n");
	for (int k_round = 0; k_round < maxK_round; ++k_round) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);

			for (int t_round = 0; t_round < maxT_round; ++t_round) {
				int top_t = tMIPs[t_round];
				int hits = get_hits(top_k, top_t, R[i], list);

				pre[t_round][k_round]    += hits / (float) top_k;
				recall[t_round][k_round] += hits / (float) top_t;
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxT_round; ++t_round) {
		int top_t = tMIPs[t_round];
		printf("Top-%d\t\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxK_round; ++k_round) {
			int top_k = kMIPs[k_round];
			pre[t_round][k_round]    = pre[t_round][k_round]    * 100.0f / qn;
			recall[t_round][k_round] = recall[t_round][k_round] * 100.0f / qn;

			printf("%4d\t\t%.2f\t\t%.2f\n", top_k,
				recall[t_round][k_round], pre[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", top_k,
				recall[t_round][k_round], pre[t_round][k_round]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;

	for (int i = 0; i < maxT_round; ++i) {
		delete[] pre[i];	pre[i] = NULL;
		delete[] recall[i];	recall[i] = NULL;
	}
	delete[] pre;	pre = NULL;
	delete[] recall;	recall = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision recall curve of simple_lsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	Simple_LSH *lsh = new Simple_LSH();
	lsh->build(n, d, K, nn_ratio, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Simple_LSH
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%ssimple_lsh_precision_recall.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int tMIPs[] = { 1, 2, 5, 10 };
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxT_round = 4;
	int maxK_round = 16;

	float **pre    = new float*[maxT_round];
	float **recall = new float*[maxT_round];
	
	for (int t_round = 0; t_round < maxT_round; ++t_round) {
		pre[t_round]    = new float[maxK_round];
		recall[t_round] = new float[maxK_round];

		for (int k_round = 0; k_round < maxK_round; ++k_round) {
			pre[t_round][k_round]    = 0;
			recall[t_round][k_round] = 0;
		}
	}

	printf("Top-t c-AMIP of Simple_LSH: \n");
	for (int k_round = 0; k_round < maxK_round; ++k_round) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);

			for (int t_round = 0; t_round < maxT_round; ++t_round) {
				int top_t = tMIPs[t_round];
				int hits = get_hits(top_k, top_t, R[i], list);

				pre[t_round][k_round]    += hits / (float) top_k;
				recall[t_round][k_round] += hits / (float) top_t;
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxT_round; ++t_round) {
		int top_t = tMIPs[t_round];
		printf("Top-%d\t\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxK_round; ++k_round) {
			int top_k = kMIPs[k_round];
			pre[t_round][k_round]    = pre[t_round][k_round]    * 100.0f / qn;
			recall[t_round][k_round] = recall[t_round][k_round] * 100.0f / qn;

			printf("%4d\t\t%.2f\t\t%.2f\n", top_k,
				recall[t_round][k_round], pre[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", top_k,
				recall[t_round][k_round], pre[t_round][k_round]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	delete[] R; R = NULL;
	for (int i = 0; i < maxT_round; ++i) {
		delete[] pre[i];	pre[i] = NULL;
		delete[] recall[i];	recall[i] = NULL;
	}
	delete[] pre;	 pre = NULL;
	delete[] recall; recall = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int norm_distribution(				// analyse norm distribution of data
	int   n,							// number of data points
	int   d,							// dimension of space
	const float **data,					// data set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  calc norm for all data objects
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	vector<float> norm(n, 0.0f);
	float max_norm = MINREAL;

	for (int i = 0; i < n; ++i) {
		norm[i] = sqrt(calc_inner_product(d, data[i], data[i]));
		if (norm[i] > max_norm) max_norm = norm[i];
	}

	// -------------------------------------------------------------------------
	//  get the percentage of frequency of norm
	// -------------------------------------------------------------------------
	int m = 25;
	float interval = max_norm / m;
	printf("m = %d, max_norm = %f, interval = %f\n", m, max_norm, interval);

	vector<int> freq(m, 0);
	for (int i = 0; i < n; ++i) {
		int id = (int) ceil(norm[i] / interval) - 1;
		if (id < 0) id = 0;
		if (id >= m) id = m - 1;
		freq[id]++;
	}

	// -------------------------------------------------------------------------
	//  write norm distribution
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%snorm_distribution.out", output_folder);

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	float num = 0.5f / m; 
	float step = 1.0f / m;
	for (int i = 0; i < m; ++i) {
		fprintf(fp, "%.1f\t%f\n", (num + step * i) * 100.0, freq[i] * 100.0 / n);
	}
	fprintf(fp, "\n");
	fclose(fp);

	gettimeofday(&end_time, NULL);
	float runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
		start_time.tv_usec) / 1000000.0f;
	printf("Norm distribution: %.6f Seconds\n\n", runtime);

	return 0;
}