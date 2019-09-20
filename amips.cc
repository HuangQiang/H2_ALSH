#include "headers.h"

// -----------------------------------------------------------------------------
int linear_scan(					// k-MIP search by linear scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%slinear.mip", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  k-MIP search by linear scan
	// -------------------------------------------------------------------------
	printf("Top-k MIP of Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		
		int top_k = TOPK[num];
		Result **result = new Result*[qn];
		for (int i = 0; i < qn; ++i) {
			result[i] = new Result[top_k];
		}
		k_mip_search(n, qn, d, top_k, data, norm_d, query, norm_q, result);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			g_recall += calc_recall(top_k, R[i], (const Result *) result[i]);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / result[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		for (int i = 0; i < qn; ++i) {
			delete[] result[i]; result[i] = NULL;
		}
		delete[] result; result = NULL;

		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
int l2_alsh(						// k-MIP search by l2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   m,							// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data objects
	const float **query,				// query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sl2_alsh.mip", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	L2_ALSH *lsh = new L2_ALSH(n, d, m, U, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  k-MIP search by L2_ALSH
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of L2_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
int l2_alsh2(						// k-MIP search by l2_alsh2
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   m,							// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data objects
	const float **query,				// query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sl2_alsh2.mip", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	L2_ALSH2 *lsh = new L2_ALSH2(n, qn, d, m, U, nn_ratio, fp, data, query);

	// -------------------------------------------------------------------------
	//  k-MIP search by L2_ALSH2
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of L2_ALSH2: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
int xbox(							// k-MIP search by xbox
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	const float **data,					// data objects
	const float **query,				// query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sxbox.mip", out_path);
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	XBox *lsh = new XBox(n, d, nn_ratio, fp, data);

	// -------------------------------------------------------------------------
	//  k-MIP search by XBox
	// -------------------------------------------------------------------------
	printf("Top-k c-AMIP of XBox: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, false, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
	//  k-MIP search by H2-ALSH-
	// -------------------------------------------------------------------------	
	sprintf(output_set, "%sh2_alsh_minus.mip", out_path);
	fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	printf("Top-k c-AMIP of H2-ALSH-: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, true, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
int sign_alsh(						// k-MIP search by sign_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	const float **data,					// data objects
	const float **query,				// query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%ssign_alsh.mip", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	Sign_ALSH *lsh = new Sign_ALSH(n, d, K, m, U, fp, data);

	// -------------------------------------------------------------------------
	//  k-MIP search by Sign_ALSH
	// -------------------------------------------------------------------------
	printf("Top-k c-AMIP of Sign_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
int simple_lsh(						// k-MIP search by simple_lsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	const float **data,					// data objects
	const float **query,				// query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%ssimple_lsh.mip", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	Simple_LSH *lsh = new Simple_LSH(n, d, K, fp, data);

	// -------------------------------------------------------------------------
	//  k-MIP search by Simple_LSH
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of Simple_LSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
int h2_alsh(						// k-MIP search by h2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R,					// MIP ground truth results
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sh2_alsh.mip", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	H2_ALSH *lsh = new H2_ALSH(n, d, nn_ratio, mip_ratio, fp, data, norm_d);

	// -------------------------------------------------------------------------
	//  k-MIP search by H2_ALSH
	// -------------------------------------------------------------------------	
	printf("Top-k c-AMIP of H2_ALSH: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], norm_q[i], list);
			g_recall += calc_recall(top_k, R[i], list);
			
			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

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
