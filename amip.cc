#include "headers.h"


// -----------------------------------------------------------------------------
int ground_truth(					// output the ground truth results
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float **data,						// data set
	float **query,						// query set
	char *truth_set)					// address of ground truth file
{
	clock_t startTime = (clock_t) -1;
	clock_t endTime   = (clock_t) -1;

	int i, j;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  output ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	int maxk = MAXK;
	float ip = -1.0F;
	MaxK_List* list = new MaxK_List(5 * MAXK);

	fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	fprintf(fp, "%d %d\n", qn, maxk);
	for (i = 0; i < qn; i++) {
		// ---------------------------------------------------------------------
		//  find k-mip points of query
		// ---------------------------------------------------------------------
		list->reset();
		for (j = 0; j < n; j++) {	
			ip = calc_inner_product(data[j], query[i], d);
			list->insert(ip, j + 1);
		}

		int real_k = maxk;
		for (j = real_k; j < 4 * MAXK; j++) {
			if (list->ith_largest_key(j) == list->ith_largest_key(maxk - 1)) {
				real_k++;
			}
		}
		
		// ---------------------------------------------------------------------
		//  output real k-mip results (consider tie)
		// ---------------------------------------------------------------------
		fprintf(fp, "%d %d", i + 1, real_k);
		for (j = 0; j < real_k; j++) {
			fprintf(fp, " %f", list->ith_largest_key(j));
		}
		int rank = 1;
		for (j = 0; j < real_k; j++) {
			if (j > 0 && list->ith_largest_key(j - 1) > list->ith_largest_key(j)) {
				rank = j + 1;
			}
			fprintf(fp, " %d %d", list->ith_largest_id(j), rank);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	endTime = clock();
	printf("Ground Truth (Linear Scan): %.6f seconds\n\n",
		((float)endTime - startTime) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int l2_alsh(						// mip search via l2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || m < 0 || U < 0.0f || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t)-1;
	clock_t endTime = (clock_t)-1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init L2_ALSH
	// -------------------------------------------------------------------------
	L2_ALSH *lsh = new L2_ALSH();
	lsh->init(n, d, m, U, nn_ratio, data);

	char output_set[200];			// output file name
	int  kMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int  maxRound = 4;				// max round
	int  top_k = -1;				// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int   thisCand = 0;
	int   allCand = 0;

	// -------------------------------------------------------------------------
	//  MIP search via L2_ALSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "l2_alsh.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of L2_ALSH: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int l2_alsh2(						// mip search via l2_alsh2
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder) 				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || m < 0 || U < 0.0f || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t)-1;
	clock_t endTime = (clock_t)-1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init L2_ALSH2
	// -------------------------------------------------------------------------
	L2_ALSH2 *lsh = new L2_ALSH2();
	lsh->init(n, qn, d, m, U, nn_ratio, data, query);

	char output_set[200];			// output file name
	int kMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int maxRound = 4;				// max round
	int top_k = -1;					// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int thisCand = 0;
	int allCand = 0;

	// -------------------------------------------------------------------------
	//  MIP search via L2_ALSH2
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "l2_alsh2.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of L2_ALSH2: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int xbox(							// mip search via xbox
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t)-1;
	clock_t endTime = (clock_t)-1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init XBox
	// -------------------------------------------------------------------------
	XBox *lsh = new XBox();
	lsh->init(n, d, nn_ratio, data);

	char output_set[200];			// output file name
	int kMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int maxRound = 4;				// max round
	int top_k = -1;					// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int thisCand = 0;
	int allCand = 0;

	// -------------------------------------------------------------------------
	//  MIP search via XBox
	// -------------------------------------------------------------------------
	strcpy(output_set, output_folder);
	strcat(output_set, "xbox.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of XBox: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, false, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  MIP search via XBox2
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "xbox2.out");

	fp = fopen(output_set, "w");	// create output file
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of XBox2: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, true, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int h2_alsh(						// mip search via h2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t)-1;
	clock_t endTime = (clock_t)-1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init H2_ALSH
	// -------------------------------------------------------------------------
	H2_ALSH *lsh = new H2_ALSH();
	lsh->init(n, d, nn_ratio, mip_ratio, data);

	char output_set[200];			// output file name
	int kMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int maxRound = 4;				// max round
	int top_k = -1;					// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int thisCand = 0;
	int allCand = 0;

	// -------------------------------------------------------------------------
	//  MIP search via H2_ALSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "h2_alsh.out");

	fp = fopen(output_set, "w");	// create output file
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of H2_ALSH: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int sign_alsh(						// mip search via sign_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	int m,								// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || K <= 0 || m < 0 || U < 0.0f || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t)-1;
	clock_t endTime = (clock_t)-1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init Sign_ALSH
	// -------------------------------------------------------------------------
	Sign_ALSH *lsh = new Sign_ALSH();
	lsh->init(n, d, K, m, U, nn_ratio, data);

	char output_set[200];			// output file name
	int kMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int maxRound = 4;				// max round
	int top_k = -1;					// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int thisCand = 0;
	int allCand = 0;

	// -------------------------------------------------------------------------
	//  MIP search via Sign_ALSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "sign_alsh.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of Sign_ALSH: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;

		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh(						// mip search via simple_lsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || K <= 0 || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t)-1;
	clock_t endTime = (clock_t)-1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init Simple_LSH
	// -------------------------------------------------------------------------
	Simple_LSH *lsh = new Simple_LSH();
	lsh->init(n, d, K, nn_ratio, data);

	char output_set[200];			// output file name
	int kMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int maxRound = 4;				// max round
	int top_k = -1;					// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int thisCand = 0;
	int allCand = 0;

	// -------------------------------------------------------------------------
	//  MIP search via Simple_LSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "simple_lsh.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-k c-AMIP of Simple_LSH: \n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock();
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);

		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			thisCand = lsh->kmip(query[i], top_k, list);
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear scan
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimensionality
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL) {
		error("parameters error in function read_file.", true);
	}

	clock_t startTime = (clock_t) -1;
	clock_t endTime   = (clock_t) -1;

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  MIP search via Linear_Scan method
	// -------------------------------------------------------------------------
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "linear.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	int kMIPs[] = {1, 2, 5, 10};	// different top-k values
	int maxRound = 4;				// max round
	int top_k = -1;					// tmp vars

	float allTime = -1.0f;
	float allRatio = -1.0f;
	float thisRatio = -1.0f;

	float allRecall = -1.0f;
	float thisRecall = -1.0f;
	int thisCand = 0;
	int allCand = 0;

	printf("Linear Scan Search:\n");
	printf("    Top-k\tRatio\t\tCandidates\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; num++) {
		startTime = clock(); 
		top_k = kMIPs[num];
		MaxK_List* list = new MaxK_List(top_k);
		
		allRatio = 0.0f;
		allRecall = 0.0f;
		allCand = 0;
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			for (j = 0; j < n; j++) {
				float ip = calc_inner_product(data[j], query[i], d);
				list->insert(ip, j + 1);
			}
			thisCand = n;
			allCand += thisCand;

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += list->ith_largest_key(j) / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;

			thisRecall = calc_recall(top_k, list, id_rank[i]);
			allRecall += thisRecall;
		}
		delete list; list = NULL;
		endTime  = clock();
		allTime = ((float)endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime = (allTime * 1000.0f) / qn;
		allRecall = (allRecall * 100.0f) / qn;
		allCand = (int)ceil((float)allCand / (float)qn);

		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\t\t%.2f%%\n",
			top_k, allRatio, allCand, allTime, allRecall);
		fprintf(fp, "%d\t%f\t%d\t%f\t%f\n",
			top_k, allRatio, allCand, allTime, allRecall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int l2_alsh_precision_recall(		// precision recall curve of l2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh
	float U,							// param of l2_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || m < 0 || U < 0.0f || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init L2_ALSH
	// -------------------------------------------------------------------------
	L2_ALSH *lsh = new L2_ALSH();
	lsh->init(n, d, m, U, nn_ratio, data);

	char output_set[200];			// output file name
	int tMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxTRound = 4;				// max round
	int maxKRound = 16;
	int maxK = kMIPs[maxKRound - 1];
	int top_t = -1;					// tmp vars

	float** allPrecision = NULL;
	float** allRecall = NULL;

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of L2_ALSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "l2_alsh_precision_recall.out");

	fp = fopen(output_set, "w");	// create output file
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of L2_ALSH: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];

		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
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
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;
	delete[] allPrecision;	allPrecision = NULL;
	delete[] allRecall;	allRecall = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int l2_alsh2_precision_recall(		// precision recall curve of l2_alsh2
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int m,								// param of l2_alsh2
	float U,							// param of l2_alsh2
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || m < 0 || U < 0.0f || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Init L2_ALSH2
	// -------------------------------------------------------------------------
	L2_ALSH2 *lsh = new L2_ALSH2();
	lsh->init(n, qn, d, m, U, nn_ratio, data, query);

	char output_set[200];			// output file name
	int tMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxTRound = 4;				// max round
	int maxKRound = 16;
	int maxK = kMIPs[maxKRound - 1];
	int top_t = -1;					// tmp vars

	float** allPrecision = NULL;
	float** allRecall = NULL;

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of L2_ALSH2
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "l2_alsh2_precision_recall.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of L2_ALSH2: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
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
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;
	delete[] allPrecision;	allPrecision = NULL;
	delete[] allRecall;	allRecall = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int xbox_precision_recall(			// precision recall curve of via xbox
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL) {
		error("parameters error in function read_file.", true);
	}

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Init Xbox
	// -------------------------------------------------------------------------
	XBox *lsh = new XBox();
	lsh->init(n, d, nn_ratio, data);

	char output_set[200];			// output file name
	int tMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxTRound = 4;				// max round
	int maxKRound = 16;
	int maxK = kMIPs[maxKRound - 1];
	int top_t = -1;					// tmp vars

	float** allPrecision = NULL;
	float** allRecall = NULL;

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of XBox
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "xbox_precision_recall.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of XBox: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, false, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of XBox2
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "xbox2_precision_recall.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of XBox2: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, true, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
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
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	for (int i = 0; i < maxTRound; i++) {
		delete[] allPrecision[i];	allPrecision[i] = NULL;
		delete[] allRecall[i];	allRecall[i] = NULL;
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;
	delete[] allPrecision;	allPrecision = NULL;
	delete[] allRecall;	allRecall = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int h2_alsh_precision_recall(		// precision recall curve of h2_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	float nn_ratio,						// approximation ratio for nn search
	float mip_ratio,					// approximation ratio for mip search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL) {
		error("parameters error in function read_file.", true);
	}

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Init H2_ALSH
	// -------------------------------------------------------------------------
	H2_ALSH *lsh = new H2_ALSH();
	lsh->init(n, d, nn_ratio, mip_ratio, data);

	char output_set[200];			// output file name
	int tMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxTRound = 4;				// max round
	int maxKRound = 16;
	int maxK = kMIPs[maxKRound - 1];
	int top_t = -1;					// tmp vars

	float** allPrecision = NULL;
	float** allRecall = NULL;

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of H2_ALSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "h2_alsh_precision_recall.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of H2_ALSH: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
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
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	for (int i = 0; i < maxTRound; i++) {
		delete[] allPrecision[i];	allPrecision[i] = NULL;
		delete[] allRecall[i];	allRecall[i] = NULL;
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;
	delete[] allPrecision;	allPrecision = NULL;
	delete[] allRecall;	allRecall = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int sign_alsh_precision_recall(		// precision recall curve of sign_alsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	int m,								// param of sign_alsh
	float U,							// param of sign_alsh
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || K <= 0 || m < 0 || U < 0.0f || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init Sign_ALSH
	// -------------------------------------------------------------------------
	Sign_ALSH *lsh = new Sign_ALSH();
	lsh->init(n, d, K, m, U, nn_ratio, data);

	char output_set[200];			// output file name
	int tMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxTRound = 4;				// max round
	int maxKRound = 16;
	int maxK = kMIPs[maxKRound - 1];
	int top_t = -1;					// tmp vars

	float** allPrecision = NULL;
	float** allRecall = NULL;

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Sign-ALSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "sign_alsh_precision_recall.out");

	fp = fopen(output_set, "w");
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of Sign_ALSH: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
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
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	for (int i = 0; i < maxTRound; i++) {
		delete[] allPrecision[i];	allPrecision[i] = NULL;
		delete[] allRecall[i];	allRecall[i] = NULL;
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;
	delete[] allPrecision;	allPrecision = NULL;
	delete[] allRecall;	allRecall = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int simple_lsh_precision_recall(	// precision recall curve of simple_lsh
	int n,								// number of data points
	int qn,								// number of query points
	int d,								// dimension of space
	int K,								// number of hash tables
	float nn_ratio,						// approximation ratio for nn search
	float **data,						// data set
	float **query,						// query set
	char *truth_set,					// address of ground truth file
	char *output_folder)				// output folder
{
	if (n < 0 || qn < 0 || d < 0 || K <= 0 || nn_ratio < 1.0f) {
		error("parameters error in function read_file.", true);
	}
	if (truth_set == NULL || output_folder == NULL || data == NULL || query == NULL) {
		error("parameters error in function read_file.", true);
	}

	int i, j, maxk = MAXK;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	fp = fopen(truth_set, "r");
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	float** R = new float*[qn];
	map<int, int>* id_rank = new map<int, int>[qn];
	int real_k, id, rank;
	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (i = 0; i < qn; i++) {
		fscanf(fp, "%d %d", &j, &real_k);
		R[i] = new float[real_k];

		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %f", &R[i][j]);
		}
		for (j = 0; j < real_k; j++) {
			fscanf(fp, " %d %d", &id, &rank);
			id_rank[i].insert(make_pair(id, rank));
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  init Simple_LSH
	// -------------------------------------------------------------------------
	Simple_LSH *lsh = new Simple_LSH();
	lsh->init(n, d, K, nn_ratio, data);

	char output_set[200];			// output file name
	int tMIPs[] = { 1, 2, 5, 10 };	// different top-k values
	int kMIPs[] = { 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000 };
	int maxTRound = 4;				// max round
	int maxKRound = 16;
	int maxK = kMIPs[maxKRound - 1];
	int top_t = -1;					// tmp vars

	float** allPrecision = NULL;
	float** allRecall = NULL;

	// -------------------------------------------------------------------------
	//  Precision Recall Curve of Simple_LSH
	// -------------------------------------------------------------------------	
	strcpy(output_set, output_folder);
	strcat(output_set, "simple_lsh_precision_recall.out");

	fp = fopen(output_set, "w");	// create output file
	if (!fp) error("Could not create output file.", true);

	printf("Top-t c-AMIP of Simple_LSH: \n");

	allPrecision = new float*[maxTRound];
	allRecall = new float*[maxTRound];
	for (int t_round = 0; t_round < maxTRound; t_round++) {
		allPrecision[t_round] = new float[maxKRound];
		allRecall[t_round] = new float[maxKRound];
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = 0;
			allRecall[t_round][k_round] = 0;
		}
	}

	for (int k_round = 0; k_round < maxKRound; k_round++) {
		int top_k = kMIPs[k_round];
		MaxK_List* list = new MaxK_List(top_k);
		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  find topk-mip points of query
			// -----------------------------------------------------------------
			list->reset();
			lsh->kmip(query[i], top_k, list);

			for (int t_round = 0; t_round < maxTRound; t_round++) {
				int hits = get_hits(top_k, tMIPs[t_round], list, id_rank[i]);

				allPrecision[t_round][k_round] += hits / (float)kMIPs[k_round];
				allRecall[t_round][k_round] += hits / (float)tMIPs[t_round];
			}
		}
		delete list; list = NULL;
	}

	for (int t_round = 0; t_round < maxTRound; t_round++) {
		top_t = tMIPs[t_round];
		printf("    Top-%d\tRecall\t\tPrecision\n", top_t);
		fprintf(fp, "Top-%d\tRecall\t\tPrecision\n", top_t);
		
		for (int k_round = 0; k_round < maxKRound; k_round++) {
			allPrecision[t_round][k_round] = (allPrecision[t_round][k_round] * 100.0f) / qn;
			allRecall[t_round][k_round] = (allRecall[t_round][k_round] * 100.0f) / qn;

			printf("    %4d\t%.2f\t\t%.2f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
			fprintf(fp, "%d\t%f\t%f\n", kMIPs[k_round],
				allRecall[t_round][k_round], allPrecision[t_round][k_round]);
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
	for (int i = 0; i < qn; i++) {
		delete[] R[i];  R[i] = NULL;
		id_rank[i].clear();
	}
	for (int i = 0; i < maxTRound; i++) {
		delete[] allPrecision[i];	allPrecision[i] = NULL;
		delete[] allRecall[i];	allRecall[i] = NULL;
	}
	delete[] R; R = NULL;
	delete[] id_rank; id_rank = NULL;
	delete[] allPrecision;	allPrecision = NULL;
	delete[] allRecall;	allRecall = NULL;

	return 0;
}
