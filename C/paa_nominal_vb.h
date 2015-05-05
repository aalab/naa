/* compute ELBO for nominal paa */
int computeCost(double **nFeatSam, double **matFeatSam, double *matSamLat, double *matLatSam, double priorMatSamLat, double priorMatLatSam, int *nd, int K, int *noCat, double *obj)
{
    double *psiMatLatSam = new double[nd[0] * K];
    double *psiMatSamLat = new double[nd[0] * K];
    double *sumMatLatSam = new double[nd[0]];
    double *sumMatSamLat = new double[K];

    double eps = 0.0000000000000001;
    for (int i = 0; i < nd[0] * K; i++) {
	psiMatSamLat[i] = digamma(matSamLat[i]);
	psiMatLatSam[i] = digamma(matLatSam[i]);
    }
    //printMatrix(psiMatSamLat, nd[0], K);
    //printMatrix(psiMatLatSam, K, nd[0]);

    double value = 0;
    double *sum = new double;

    colSum(matSamLat, nd[0], K, sumMatSamLat);
    colSum(matLatSam, K, nd[0], sumMatLatSam);

    colSub(psiMatSamLat, nd[0], K, sumMatSamLat);
    colSub(psiMatLatSam, K, nd[0], sumMatLatSam);

    colSum(psiMatLatSam, K, nd[0], sumMatLatSam);
    colSum(psiMatSamLat, nd[0], K, sumMatSamLat);
    colSum(sumMatLatSam, nd[0], 1, sum);
    value = value + (priorMatLatSam - 1) * *sum;
    //cout << priorMatLatSam << " " << value << endl;
    colSum(sumMatSamLat, K, 1, sum);
    value = value + (priorMatSamLat - 1) * *sum;
    //cout << "value " << value << endl;

    for (int i = 0; i < nd[0] * K; i++) {
	psiMatSamLat[i] = exp(psiMatSamLat[i]);
	psiMatLatSam[i] = exp(psiMatLatSam[i]);
    }

    for (int i = 0; i < nd[1]; i++) {
	double *tmp1 = new double[noCat[i] * K];
	double *tmp2 = new double[noCat[i] * nd[0]];

	cblas_dgemm(CblasRowMajor,
		    CblasNoTrans, CblasNoTrans, noCat[i], K, nd[0],
	       1.0, matFeatSam[i], nd[0], psiMatSamLat, K, 0.0, tmp1, K);

	cblas_dgemm(CblasRowMajor,
		    CblasNoTrans, CblasNoTrans, noCat[i], nd[0], K,
		    1.0, tmp1, K, psiMatLatSam, nd[0], 0.0, tmp2, nd[0]);

	for (int j = 0; j < noCat[i] * nd[0]; j++) {
	    value = value + nFeatSam[i][j] * log(eps + tmp2[j]);
	}

	delete[] tmp1;
	delete[] tmp2;
    }
    //cout << "value " << value << endl;

    entDir(matLatSam, K, nd[0], sumMatLatSam);
    //cout << "ent sam " << printMatrix(sumMatLatSam, 1, nd[0]);
    colSum(sumMatLatSam, nd[0], 1, sum);
    //cout << "ent sam sum " << *sum << endl;
    value = value + *sum;
    entDir(matSamLat, nd[0], K, sumMatSamLat);
    //cout << "ent lat " << printMatrix(sumMatSamLat, 1, K);
    colSum(sumMatSamLat, K, 1, sum);
    //cout << "ent lat sum " << *sum << endl;
    value = value + *sum;
    //cout << "value " << value << endl;

    *obj = value;

    delete[] psiMatLatSam;
    delete[] psiMatSamLat;
    delete[] sumMatLatSam;
    delete[] sumMatSamLat;
    delete sum;
    return 0;
}

/* VB iterations for nominal valued paa */
int run_paa(double **nFeatSam, double **matFeatSam, int *noCat, int *nd, int K, double *matLatSam, double *matSamLat, double *obj, double thresh, double priorMatLatSam, double priorMatSamLat, int maxIter)
{
/* temporary storage for intermediate values tmp1 nxcat tmp2 Kxcat
       tmp3 nxK */
    double **tmp1 = new double *[nd[1]];
    for (int i = 0; i < nd[1]; ++i) {
	tmp1[i] = new double[nd[0] * noCat[i]];
    }
    double **tmp2 = new double *[nd[1]];
    for (int i = 0; i < nd[1]; ++i) {
	tmp2[i] = new double[K * noCat[i]];
    }
    double *tmp3 = new double[nd[0] * K];

    double *psiMatLatSam = new double[nd[0] * K];
    double *psiMatSamLat = new double[nd[0] * K];
    double *sumMatLatSam = new double[nd[0]];
    double *sumMatSamLat = new double[K];

    cout << "Starting iterations, ";
    for (int iter = 0; iter < maxIter; iter++) {
	if (iter % 100 == 0) {
	    cout << ".";
	}
	//cout << "matLatSam" << endl;
	//printMatrix(matLatSam, K, nd[0], 12, 12);
	//cout << "matSamLat" << endl;
	//printMatrix(matSamLat, nd[0], K, 12, 12);

	for (int i = 0; i < nd[0] * K; i++) {
	    psiMatSamLat[i] = digamma(matSamLat[i]);
	    psiMatLatSam[i] = digamma(matLatSam[i]);
	}
	//printMatrix(psiMatLatSam, K, nd[0]);
	//printMatrix(psiMatSamLat, nd[0], K);

	colSum(matSamLat, nd[0], K, sumMatSamLat);
	colSum(matLatSam, K, nd[0], sumMatLatSam);
	//printMatrix(sumMatLatSam, 1, nd[0]);
	//printMatrix(sumMatSamLat, 1, K);

	colSub(psiMatSamLat, nd[0], K, sumMatSamLat);
	colSub(psiMatLatSam, K, nd[0], sumMatLatSam);
	for (int i = 0; i < nd[0] * K; i++) {
	    psiMatSamLat[i] = exp(psiMatSamLat[i]);
	    psiMatLatSam[i] = exp(psiMatLatSam[i]);
	}
	//cout << "psiMatLatSam" << endl;
	//printMatrix(psiMatLatSam, K, nd[0]);
	//cout << "psiMatSamLat" << endl;
	//printMatrix(psiMatSamLat, nd[0], K);


	for (int i = 0; i < nd[0] * K; i++) {
	    matSamLat[i] = 0;
	    matLatSam[i] = 0;
	}

	for (int i = 0; i < nd[1]; i++) {
	    /* ...(.., .., .., nRow_final_matrix, nCol_final_matrix,
	       mult_dim, alpha, A, nCol_A, B, nCol_B, beta, C, nCol_C) */

	    /* tmp2 = matFeatSam * psiMatSamLat */
	    cblas_dgemm(CblasRowMajor,
			CblasNoTrans, CblasNoTrans, noCat[i], K, nd[0],
	    1.0, matFeatSam[i], nd[0], psiMatSamLat, K, 0.0, tmp2[i], K);
	    //cout << "nFeatSam[i]" << endl;
	    //printMatrix(nFeatSam[i], noCat[i], nd[0]);
	    //cout << "matFeatSam[i] * psiMatSamLat" << endl;
	    //printMatrix(tmp2[i], noCat[i], K);

	    /* tmp1 = tmp2 * psiMatLatSam */
	    cblas_dgemm(CblasRowMajor,
			CblasNoTrans, CblasNoTrans, noCat[i], nd[0], K,
	      1.0, tmp2[i], K, psiMatLatSam, nd[0], 0.0, tmp1[i], nd[0]);
	    //printMatrix(tmp1[i], noCat[i], nd[0]);

	    for (int j = 0; j < noCat[i] * nd[0]; j++) {
		tmp1[i][j] = nFeatSam[i][j] / (0.0000000000000001 + tmp1[i][j]);
	    }
	    //cout << "tmp" << endl;
	    //printMatrix(tmp1[i], noCat[i], nd[0]);

	    /* tmp2 = tmp1 * psiMatLatSam' */
	    cblas_dgemm(CblasRowMajor,
			CblasNoTrans, CblasTrans, noCat[i], K, nd[0],
	      1.0, tmp1[i], nd[0], psiMatLatSam, nd[0], 0.0, tmp2[i], K);

	    /* tmp3 = matFeatSam' * tmp2 */
	    cblas_dgemm(CblasRowMajor,
			CblasTrans, CblasNoTrans, nd[0], K, noCat[i],
		    1.0, matFeatSam[i], nd[0], tmp2[i], K, 0.0, tmp3, K);

	    for (int j = 0; j < K * nd[0]; j++) {
		matSamLat[j] = matSamLat[j] + tmp3[j] * psiMatSamLat[j];
	    }

	    /* tmp2 = psiMatSamLat' * matFeatSam' */
	    cblas_dgemm(CblasRowMajor,
			CblasTrans, CblasTrans, K, noCat[i], nd[0],
			1.0, psiMatSamLat, K, matFeatSam[i], nd[0], 0.0, tmp2[i], noCat[i]);

	    /* tmp3 = tmp2 * tmp1 */
	    cblas_dgemm(CblasRowMajor,
			CblasNoTrans, CblasNoTrans, K, nd[0], noCat[i],
	       1.0, tmp2[i], noCat[i], tmp1[i], nd[0], 0.0, tmp3, nd[0]);

	    for (int j = 0; j < K * nd[0]; j++) {
		matLatSam[j] = matLatSam[j] + tmp3[j] * psiMatLatSam[j];
	    }

	}

	for (int j = 0; j < K * nd[0]; j++) {
	    matSamLat[j] = priorMatSamLat + matSamLat[j];
	    matLatSam[j] = priorMatLatSam + matLatSam[j];
	}

	if (thresh != -1) {
	    computeCost(nFeatSam, matFeatSam, matSamLat, matLatSam,
	       priorMatSamLat, priorMatLatSam, nd, K, noCat, &obj[iter]);
	    if (iter > 1) {
		if (fabs(obj[iter] - obj[iter -
				    1]) / fabs(obj[iter - 1]) < thresh) {
		    cout << iter + 1 << ", required threshold reached with obj: " << obj[iter] << endl;
		    break;
		}
	    }
	}
	if (iter == maxIter - 1) {
	    if (thresh == -1) {
		computeCost(nFeatSam, matFeatSam, matSamLat, matLatSam,
		  priorMatSamLat, priorMatLatSam, nd, K, noCat, &obj[0]);
		cout << iter + 1 << ", maximum iteration reached with obj: " << obj[0] << endl;
	    } else {
		cout << iter + 1 << ", maximum iteration reached with obj: " << obj[iter] << endl;
	    }
	}
	//cout << "obj" << obj[iter] << endl;

    }
    cout << endl;
    delete[] tmp3;

    delete[] psiMatLatSam;
    delete[] psiMatSamLat;
    delete[] sumMatLatSam;
    delete[] sumMatSamLat;

    for (int i = 0; i < nd[1]; ++i) {
	delete[] tmp1[i];
    }
    delete[] tmp1;
    for (int i = 0; i < nd[1]; ++i) {
	delete[] tmp2[i];
    }
    delete[] tmp2;
    return 1;
}
