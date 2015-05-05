#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string.h>
#include<math.h>
#include<fstream>
#include </opt/local/include/gsl/gsl_cblas.h>
#include "supportingFunctions.h"
#include "paa_nominal_vb.h"
using namespace std;

/* this script runs variational Bayes' update rule for archetypal analysis on nominal observations
example: time ./paa --maxIter 100 --k 20 --priorMatLatSam 0.3 --priorMatSamLat 0.1 --thresh -1 --randSeed 4 --d <d> --n <n>--input input-file --output output-file
where input-file is a space separated text file with each column a feature and each row a sample
the script produces three files output-file.* with extensions W, H and O for W, H, and objective values
description of optional parameters:
--d number of features, if n or d are missing they are calculated from input-file
--n number of samples, if n or d are missing they are calculated from input-file
--k maximum number of archetypes, default 20
--maxIter maximum number of iterations, default 1000
--thresh convergence criteria, stop if relative change between successive objective value is less than thresh, default 0.0001, (for speed) if thresh is - 1 then objective function is not stored at each iteration, but only the final value
--randSeed random seed, default 1
--priorMatLatSam and --priorMatSamLat priors over H and W matrices respectively
*/

/* copyright Sohan Seth sohan.seth@gmail.com */


int main(int argc, char **argv)
{
    double *parameters = new double[8];
    parameters[0] = 20;
    //K, number of latent vectors
	parameters[1] = 1000;
    //maxIter, maximum number of iterations
	parameters[2] = 0.0001;
    //thresh, breaking point
	parameters[3] = 0.3;
    //priorMatLatSam
	parameters[4] = 0.1;
    //priorMatSamLat
	parameters[5] = 1;
    //random seed
	parameters[6] = 0;
    parameters[7] = 0;
    //d and n
	parse(argc, argv, parameters);
    int K = int (parameters[0]);
    int maxIter = int (parameters[1]);
    double thresh = parameters[2];
    double priorMatLatSam = parameters[3];
    double priorMatSamLat = parameters[4];
    int randSeed = int (parameters[5]);
    int nd[2];
    nd[1] = int (parameters[6]);
    nd[0] = int (parameters[7]);

    char file[100];		/* = argv[argc - 1] */
    char outputfile[100], opfile[100];
    parse2(argc, argv, file, outputfile);

    cout << K << " " << maxIter << " " << thresh << " " << priorMatLatSam << " " << priorMatSamLat << " " << randSeed << " " << file << " " << outputfile << endl;

    double *obj = NULL;
    if (thresh == -1) {
	obj = new double[1];
    } else {
	obj = new double[maxIter];
	for (int i = 0; i < maxIter; i++) {
	    obj[i] = 0;
	}
    }

    delete[] parameters;

    if (nd[0] == 0 || nd[1] == 0) {
	retnnd(file, nd);
    }
    //cout << "nr of samples " << nd[0] << " nr of dimension " << nd[1] << endl;

    /* load original data in a 2D matrix, each (i,j)-th entry is a
       category of j-th feature */
    int **mat = new int *[nd[0]];
    for (int i = 0; i < nd[0]; ++i) {
	mat[i] = new int[nd[1]];
    }
    readDataFile(file, mat);
    //cout << "File read complete" << endl;

    /* find nr of categories for each feature */
    int *noCat = new int[nd[1]];
    compCat(mat, nd[0], nd[1], noCat);
    //cout << "nr of categories" << endl;
    //printMatrix(noCat, 1, nd[1]);

    /* convert input matrix in binary format */
    double **nFeatSam = new double *[nd[1]];
    for (int i = 0; i < nd[1]; ++i) {
	nFeatSam[i] = new double[nd[0] * noCat[i]];
    }
    mat2bin(mat, nFeatSam, nd[0], nd[1]);
    double **matFeatSam = new double *[nd[1]];
    for (int i = 0; i < nd[1]; ++i) {
	matFeatSam[i] = new double[nd[0] * noCat[i]];
    }
    mat2bin(mat, matFeatSam, nd[0], nd[1]);
    for (int i = 0; i < nd[1]; ++i) {
	colDiv(matFeatSam[i], noCat[i], nd[0]);
	//printMatrix(matFeatSam[i], noCat[i], nd[0]);
	//cout << i << endl;
    }

    int *assignment = new int[nd[0]];
    unique(mat, nd[0], nd[1], assignment);
    /* for (int i = 0; i < nd[0]; ++i) { printMatrix(mat[i], 1, nd[1]); }
       printMatrix(assignment, nd[0], 1); */

    /* for (int i = 0; i < nd[1]; ++i) { cout << "feature " << i << "
       categories " << noCat[i] << endl; printMatrix(nFeatSam[i],
       noCat[i], nd[0]); } */
    for (int i = 0; i < nd[0]; ++i) {
	delete[] mat[i];
    }
    delete[] mat;

    double *matLatSam = new double[nd[0] * K];
    double *matSamLat = new double[nd[0] * K];

    /* initialization */
    srand(randSeed);
    for (int i = 0; i < 1000; i++) {
	rand();
    }
    for (int i = 0; i < nd[0] * K; i++) {
	matSamLat[i] = priorMatSamLat * (1 + 0.01 * rand() / (double)(RAND_MAX));
	matLatSam[i] = priorMatLatSam * (1 + 0.01 * rand() / (double)(RAND_MAX));
    }
    for (int i = nd[0] - 1; i >= 0; --i) {
	for (int j = 0; j < K; ++j) {
	    matSamLat[i * K + j] = matSamLat[assignment[i] * K + j];
	    //(1 + j);
	    //
		matLatSam[j * nd[0] + i] = matLatSam[j * nd[0] + assignment[i]];
	    //(1 + j);
	    //
	}
    }

    //printMatrix(matLatSam, K, nd[0], 5, 5);
    //cout << endl;
    //printMatrix(matSamLat, nd[0], K, 5, 5);

    run_paa(nFeatSam, matFeatSam, noCat, nd, K, matLatSam, matSamLat, obj, thresh, priorMatLatSam, priorMatSamLat, maxIter);

    //cout << "matLatSam" << endl;
    colDiv(matLatSam, K, nd[0]);
    colDiv(matSamLat, nd[0], K);

    //printMatrix(matLatSam, K, nd[0], 5, 5);
    double *maxMatLatSam = new double[K];
    maxRow(matLatSam, K, nd[0], maxMatLatSam);
    //printMatrix(maxMatLatSam, K, 1);
    cout << "active components (generating observations):";
    int totActComp = 0;
    int activeIndex[K];
    for (int i = 0; i < K; i++) {
	if (maxMatLatSam[i] > 0.15) {
	    activeIndex[i] = 1;
	    cout << endl << i + 1 << ": ";
	    for (int j = 0; j < nd[0]; j++) {
		if (matSamLat[K * j + i] > 0.01) {
		    cout << j << ",";
		}
	    }
	    totActComp = totActComp + 1;
	} else {
	    activeIndex[i] = 0;
	}
    } cout << endl << "total " << totActComp << ": ";
    printMatrix(activeIndex, 1, K);
    cout << endl;

    sprintf(opfile, "%s.%s", outputfile, "H");
    writeMatrix(matLatSam, K, nd[0], opfile, 0, activeIndex);
    //printMatrix(matLatSam, K, nd[0]);
    //cout << "matSamLat" << endl;
    sprintf(opfile, "%s.%s", outputfile, "W");
    writeMatrix(matSamLat, nd[0], K, opfile, 1, activeIndex);

    //printMatrix(matSamLat, nd[0], K);
    //printMatrix(obj, maxIter, 1);
    //cout << "final objective: " << obj[maxIter] << endl;

    /* write obj function to *.O file */
    sprintf(opfile, "%s.%s", outputfile, "O");
    if (thresh == -1) {
	writeMatrix(obj, 1, 1, opfile);
    } else {
	writeMatrix(obj, maxIter, 1, opfile);
    }

    delete[] matLatSam;
    delete[] matSamLat;

    for (int i = 0; i < nd[1]; ++i) {
	delete[] nFeatSam[i];
    }
    delete[] nFeatSam;
    for (int i = 0; i < nd[1]; ++i) {
	delete[] matFeatSam[i];
    }
    delete[] matFeatSam;
    delete[] assignment;

    delete[] noCat;
    delete[] obj;
    return 0;
}
