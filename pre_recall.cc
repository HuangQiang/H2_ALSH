#include "pre_recall.h"

// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision-recall curve of h2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	float **pre,						// precision 
	float **recall,						// recall
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%spre_recall_%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	H2_ALSH *lsh = new H2_ALSH(n, d, nn_ratio, mip_ratio, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n\n", indexing_time);

	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n", g_memory / 1048576.0f);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of H2_ALSH
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	for (int t = 0; t < MAX_T; ++t) {
		int top_t = tMIPs[t];
		MaxK_List* list = new MaxK_List(top_t);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_t, query[i], norm_q[i], list);

			for (int r = 0; r < MAX_ROUND; ++r) {
				int top_k = TOPK[r];
				int hits  = get_hits(top_t, top_k, R[i], list);

				pre[r][t]    += (float) hits / (float) top_t;
				recall[r][t] += (float) hits / (float) top_k;
			}
		}
		for (int r = 0; r < MAX_ROUND; ++r) {
			pre[r][t] = pre[r][t] * 100.0f / qn;
			recall[r][t] = recall[r][t] * 100.0f / qn;
		}
		delete list; list = NULL;
	}
	delete lsh; lsh = NULL;

	// write the precision-recall curves to disk
	write_pre_recall(fp, (const float **) pre, (const float **) recall);
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Precision-Recall Curve: %.6f Seconds\n\n", runtime);

	return 0;
}

// -----------------------------------------------------------------------------
int sign_alsh_precision_recall(		// precision-recall curve of sign_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	int   m,							// param of sign_alsh
	float U,							// param of sign_alsh
	float **pre,						// precision 
	float **recall,						// recall
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%spre_recall_%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Sign_ALSH *lsh = new Sign_ALSH(n, d, K, m, U, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n\n", indexing_time);

	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n", g_memory / 1048576.0f);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Sign-ALSH
	// -------------------------------------------------------------------------	
	gettimeofday(&g_start_time, NULL);
	for (int t = 0; t < MAX_T; ++t) {
		int top_t = tMIPs[t];
		MaxK_List* list = new MaxK_List(top_t);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_t, query[i], norm_q[i], list);

			for (int r = 0; r < MAX_ROUND; ++r) {
				int top_k = TOPK[r];
				int hits  = get_hits(top_t, top_k, R[i], list);

				pre[r][t]    += hits / (float) top_t;
				recall[r][t] += hits / (float) top_k;
			}
		}
		for (int r = 0; r < MAX_ROUND; ++r) {
			pre[r][t] = pre[r][t] * 100.0f / qn;
			recall[r][t] = recall[r][t] * 100.0f / qn;
		}
		delete list; list = NULL;
	}
	delete lsh; lsh = NULL;

	// write the precision-recall curves to disk
	write_pre_recall(fp, (const float **) pre, (const float **) recall);
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Precision-Recall Curve: %.6f Seconds\n\n", runtime);

	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision-recall curve of simple_lsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	float **pre,						// precision 
	float **recall,						// recall
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%spre_recall_%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Simple_LSH *lsh = new Simple_LSH(n, d, K, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n\n", indexing_time);

	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n", g_memory / 1048576.0f);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Simple_LSH
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	for (int t = 0; t < MAX_T; ++t) {
		int top_t = tMIPs[t];
		MaxK_List* list = new MaxK_List(top_t);

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_t, query[i], norm_q[i], list);

			for (int r = 0; r < MAX_ROUND; ++r) {
				int top_k = TOPK[r];
				int hits  = get_hits(top_t, top_k, R[i], list);

				pre[r][t]    += hits / (float) top_t;
				recall[r][t] += hits / (float) top_k;
			}
		}
		for (int r = 0; r < MAX_ROUND; ++r) {
			pre[r][t] = pre[r][t] * 100.0f / qn;
			recall[r][t] = recall[r][t] * 100.0f / qn;
		}
		delete list; list = NULL;
	}
	delete lsh; lsh = NULL;

	// write the precision-recall curves to disk
	write_pre_recall(fp, (const float **) pre, (const float **) recall);
	fclose(fp);
	
	gettimeofday(&g_end_time, NULL);
	float runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Precision-Recall Curve: %.6f Seconds\n\n", runtime);

	return 0;
}

