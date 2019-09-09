#include "headers.h"

// -----------------------------------------------------------------------------
int h2_alsh(						// mip search via h2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sh2_alsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	H2_ALSH *lsh = new H2_ALSH(n, d, nn_ratio, mip_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  c-AMIP search via H2_ALSH
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of H2_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

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
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sl2_alsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	L2_ALSH *lsh = new L2_ALSH(n, d, m, U, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  c-AMIP search via L2_ALSH
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of L2_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

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
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sl2_alsh2.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	L2_ALSH2 *lsh = new L2_ALSH2(n, qn, d, m, U, nn_ratio, fp, data, query);

	// -------------------------------------------------------------------------
	//  c-AMIP search via L2_ALSH2
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of L2_ALSH2: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

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
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sxbox.out", output_folder);
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	XBox *lsh = new XBox(n, d, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  c-AMIP search via XBox
	// -------------------------------------------------------------------------
	printf("Top-k c-AMIP of XBox: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, false, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
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
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, true, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

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
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%ssign_alsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	Sign_ALSH *lsh = new Sign_ALSH(n, d, K, m, U, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  c-AMIP search via Sign_ALSH
	// -------------------------------------------------------------------------
	printf("Top-k c-AMIP of Sign_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

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
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%ssimple_lsh.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	Simple_LSH *lsh = new Simple_LSH(n, d, K, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  c-AMIP search via Simple_LSH
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of Simple_LSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, (const Result *) R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear_scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
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

	printf("Top-k MIP of Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			for (int j = 0; j < n; ++j) {
				float ip = calc_inner_product(d, data[j], query[i]);
				list->insert(ip, j + 1);
			}
			g_recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
			g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision recall curve of h2_alsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char output_set[200];
	sprintf(output_set, "%sh2_alsh_precision_recall.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	H2_ALSH *lsh = new H2_ALSH(n, d, nn_ratio, mip_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of H2_ALSH
	// -------------------------------------------------------------------------
	for (int t = 0; t < MAX_T; ++t) {
		int top_t = tMIPs[t];
		MaxK_List* list = new MaxK_List(top_t);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_t, query[i], list);

			for (int round = 0; round < MAX_ROUND; ++round) {
				int top_k = kMIPs[round];
				int hits = get_hits(top_t, top_k, R[i], list);

				pre[round][t]    += hits / (float) top_t;
				recall[round][t] += hits / (float) top_k;
			}
		}
		delete list; list = NULL;
	}
	delete lsh; lsh = NULL;

	for (int round = 0; round < MAX_ROUND; ++round) {
		int top_k = kMIPs[round];
		printf("Top-%d\t\tRecall\t\tPrecision\n", top_k);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_k);
		
		for (int t = 0; t < MAX_T; ++t) {
			int top_t = tMIPs[t];
			pre[round][t]    = pre[round][t]    * 100.0f / qn;
			recall[round][t] = recall[round][t] * 100.0f / qn;

			printf("%4d\t\t%.2f\t\t%.2f\n", top_t, recall[round][t], pre[round][t]);
			fprintf(fp, "%d\t%f\t%f\n", top_t, recall[round][t], pre[round][t]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
		g_start_time.tv_usec) / 1000000.0f;
	printf("Precision-Recall Curve: %.6f Seconds\n\n", runtime);

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
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char output_set[200];
	sprintf(output_set, "%ssign_alsh_precision_recall.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	Sign_ALSH *lsh = new Sign_ALSH(n, d, K, m, U, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Sign-ALSH
	// -------------------------------------------------------------------------	
	for (int t = 0; t < MAX_T; ++t) {
		int top_t = tMIPs[t];
		MaxK_List* list = new MaxK_List(top_t);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_t, query[i], list);

			for (int round = 0; round < MAX_ROUND; ++round) {
				int top_k = kMIPs[round];
				int hits = get_hits(top_t, top_k, R[i], list);

				pre[round][t]    += hits / (float) top_t;
				recall[round][t] += hits / (float) top_k;
			}
		}
		delete list; list = NULL;
	}
	delete lsh; lsh = NULL;

	for (int round = 0; round < MAX_ROUND; ++round) {
		int top_k = kMIPs[round];
		printf("Top-%d\t\tRecall\t\tPrecision\n", top_k);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_k);
		
		for (int t = 0; t < MAX_T; ++t) {
			int top_t = tMIPs[t];
			pre[round][t]    = pre[round][t]    * 100.0f / qn;
			recall[round][t] = recall[round][t] * 100.0f / qn;

			printf("%4d\t\t%.2f\t\t%.2f\n", top_t, recall[round][t], pre[round][t]);
			fprintf(fp, "%d\t%f\t%f\n", top_t, recall[round][t], pre[round][t]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
		g_start_time.tv_usec) / 1000000.0f;
	printf("Precision-Recall Curve: %.6f Seconds\n\n", runtime);

	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision recall curve of simple_lsh
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   K,							// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	float **pre,						// precision 
	float **recall,						// recall
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char output_set[200];
	sprintf(output_set, "%ssimple_lsh_precision_recall.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	Simple_LSH *lsh = new Simple_LSH(n, d, K, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Simple_LSH
	// -------------------------------------------------------------------------
	for (int t = 0; t < MAX_T; ++t) {
		int top_t = tMIPs[t];
		MaxK_List* list = new MaxK_List(top_t);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_t, query[i], list);

			for (int round = 0; round < MAX_ROUND; ++round) {
				int top_k = kMIPs[round];
				int hits = get_hits(top_t, top_k, R[i], list);

				pre[round][t]    += hits / (float) top_t;
				recall[round][t] += hits / (float) top_k;
			}
		}
		delete list; list = NULL;
	}
	delete lsh; lsh = NULL;

	for (int round = 0; round < MAX_ROUND; ++round) {
		int top_k = kMIPs[round];
		printf("Top-%d\t\tRecall\t\tPrecision\n", top_k);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_k);
		
		for (int t = 0; t < MAX_T; ++t) {
			int top_t = tMIPs[t];
			pre[round][t]    = pre[round][t]    * 100.0f / qn;
			recall[round][t] = recall[round][t] * 100.0f / qn;

			printf("%4d\t\t%.2f\t\t%.2f\n", top_t, recall[round][t], pre[round][t]);
			fprintf(fp, "%d\t%f\t%f\n", top_t, recall[round][t], pre[round][t]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + (g_end_time.tv_usec - 
		g_start_time.tv_usec) / 1000000.0f;
	printf("Precision-Recall Curve: %.6f Seconds\n\n", runtime);

	return 0;
}

