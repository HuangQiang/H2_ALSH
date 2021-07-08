#include "amips.h"

namespace mips {

// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth MIP results
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// dimensionality
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const char  *truth_set) 			// address of truth set
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) { printf("Could not create %s\n", truth_set); return 1; }

	// -------------------------------------------------------------------------
	//  calc the norm of data
	// -------------------------------------------------------------------------
	Result *order_d = new Result[n];
	for (int i = 0; i < n; ++i) {
		order_d[i].id_  = i;
		order_d[i].key_ = norm_d[i][0];
	}
	qsort(order_d, n, sizeof(Result), ResultCompDesc);

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	MaxK_List *list = new MaxK_List(MAXK);

	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		float kip = MINREAL;
		list->reset();
		for (int j = 0; j < n; ++j) {
 			int id = order_d[j].id_;
			if (norm_d[id][0] * norm_q[i][0] <= kip) break;

			float ip = calc_inner_product(d, kip, data[id], norm_d[id], 
				query[i], norm_q[i]);
			kip = list->insert(ip, id + 1);
		}
		for (int j = 0; j < list->size(); ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	delete[] order_d; 
	delete   list;
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float truth_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);
	
	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// k-MIP search by linear_scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%s%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  calc the norm of data
	// -------------------------------------------------------------------------
	Result *order_d = new Result[n];
	for (int i = 0; i < n; ++i) {
		order_d[i].id_  = i;
		order_d[i].key_ = norm_d[i][0];
	}
	qsort(order_d, n, sizeof(Result), ResultCompDesc);

	// -------------------------------------------------------------------------
	//  k-MIPS of linear_scan
	// -------------------------------------------------------------------------
	printf("k-MIPS of %s:\n", method_name);
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			float kip = MINREAL;
			list->reset();
			for (int j = 0; j < n; ++j) {
				int id = order_d[j].id_;
				if (norm_d[id][0] * norm_q[i][0] <= kip) break;
				
				float ip = calc_inner_product(d, kip, data[id], norm_d[id], 
					query[i], norm_q[i]);
				kip = list->insert(ip, id + 1);
			}
			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete[] order_d;

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
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%s%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	L2_ALSH *lsh = new L2_ALSH(n, d, m, U, nn_ratio, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time:    %f Seconds\n", g_indextime);
	printf("Estimated Memory: %f MB\n\n", g_memory);

	fprintf(fp, "%s: m=%d, U=%.2f, c0=%.2f\n", method_name, m, U, nn_ratio);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  k-MIPS of l2_alsh
	// -------------------------------------------------------------------------	
	printf("k-MIPS of %s:\n", method_name);
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

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete lsh;

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
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%s%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	L2_ALSH2 *lsh = new L2_ALSH2(n, qn, d, m, U, nn_ratio, data, norm_d, norm_q);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time:    %f Seconds\n", g_indextime);
	printf("Estimated Memory: %f MB\n\n", g_memory);
	
	fprintf(fp, "%s: m=%d, U=%.2f, c0=%.2f\n", method_name, m, U, nn_ratio);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  k-MIPS of l2_alsh2
	// -------------------------------------------------------------------------	
	printf("k-MIPS of %s:\n", method_name);
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

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete lsh;

	return 0;
}

// -----------------------------------------------------------------------------
int xbox(							// k-MIP search by xbox
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	const char *method_name1,			// name of method
	const char *method_name2,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	FILE *fp = NULL;

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	XBox *xbox = new XBox(n, d, nn_ratio, data, norm_d);
	xbox->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = xbox->get_memory_usage() / 1048576.0f;

	printf("Indexing Time:    %f Seconds\n", g_indextime);
	printf("Estimated Memory: %f MB\n\n", g_memory);

	// -------------------------------------------------------------------------
	//  k-MIPS of xbox
	// -------------------------------------------------------------------------
	sprintf(output_set, "%s%s.out", out_path, method_name1);
	fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	fprintf(fp, "%s: c0=%.2f\n", method_name1, nn_ratio);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	printf("k-MIPS of %s:\n", method_name1);
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			xbox->kmip(top_k, false, query[i], norm_q[i], list);

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	//  k-MIPS of h2_alsh-
	// -------------------------------------------------------------------------	
	sprintf(output_set, "%s%s.out", out_path, method_name2);
	fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	fprintf(fp, "%s: c0=%.2f\n", method_name2, nn_ratio);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	printf("k-MIPS of %s:\n", method_name2);
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List* list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			xbox->kmip(top_k, true, query[i], norm_q[i], list);

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete xbox;

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
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%s%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Sign_ALSH *lsh = new Sign_ALSH(n, d, K, m, U, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time:    %f Seconds\n", g_indextime);
	printf("Estimated Memory: %f MB\n\n", g_memory);
	
	fprintf(fp, "%s: K=%d, m=%d, U=%.2f\n", method_name, K, m, U);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  k-MIPS of sign_alsh
	// -------------------------------------------------------------------------
	printf("k-MIPS of %s:\n", method_name);
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

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete lsh;

	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh(						// k-MIP search by simple_lsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   K,							// number of hash tables
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%s%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Simple_LSH *lsh = new Simple_LSH(n, d, K, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time:    %f Seconds\n", g_indextime);
	printf("Estimated Memory: %f MB\n\n", g_memory);

	fprintf(fp, "%s: K=%d\n", method_name, K);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  k-MIPS of simple_lsh
	// -------------------------------------------------------------------------	
	printf("k-MIPS of %s:\n", method_name);
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kmip(top_k, query[i], norm_q[i], list);

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete lsh;

	return 0;
}

// -----------------------------------------------------------------------------
int h2_alsh(						// k-MIP search by h2_alsh
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float nn_ratio,						// approximation ratio for ANN search
	float mip_ratio,					// approximation ratio for AMIP search
	const char *method_name,			// name of method
	const char *out_path,				// output path
	const float **data,					// data objects
	const float **norm_d,				// l2-norm of data objects
	const float **query,				// query objects
	const float **norm_q,				// l2-norm of query objects
	const Result **R)					// MIP ground truth results
{
	char output_set[200];
	sprintf(output_set, "%s%s.out", out_path, method_name);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	H2_ALSH *lsh = new H2_ALSH(n, d, nn_ratio, mip_ratio, data, norm_d);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time:    %f Seconds\n", g_indextime);
	printf("Estimated Memory: %f MB\n\n", g_memory);

	fprintf(fp, "%s: c=%.2f, c0=%.2f\n", method_name, mip_ratio, nn_ratio);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  k-MIPS of h2_alsh
	// -------------------------------------------------------------------------	
	printf("k-MIPS of %s:\n", method_name);
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

			g_ratio  += calc_ratio(top_k,  R[i], list);
			g_recall += calc_recall(top_k, R[i], list);
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
	delete lsh;

	return 0;
}

} // end namespace mips
