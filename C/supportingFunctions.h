using namespace std;

/* digamma function by Antti Honkela */
#define DIGAMMA_S 1e-5
#define DIGAMMA_C 8.5
#define DIGAMMA_S3 1.0/12
#define DIGAMMA_S4 1.0/120
#define DIGAMMA_S5 1.0/252
#define DIGAMMA_D1 -0.5772156649015328605195174

double digamma(double x)
{
    double y = 0.0;
    double r = 0.0;
    double xn = x;

    if (xn < 0) {
	cout << xn << " Argument of digamma must be positve" << endl;
	return 0;
	//returns nan in the original code, *modified by Sohan Seth *
    }
    if (xn <= DIGAMMA_S)
	y = DIGAMMA_D1 - 1.0 / xn;
    else {
	while (xn < DIGAMMA_C) {
	    y = y - 1.0 / xn;
	    xn = xn + 1;
	}
	r = 1.0 / xn;
	y = y + log(xn) - .5 * r;
	r = r * r;
	y = y - r * (DIGAMMA_S3 - r * (DIGAMMA_S4 - r * DIGAMMA_S5));
    }

    return y;
}

/* max over each column */
int maxCol(int **array, int nRow, int nCol, int *max)
{
    for (int i = 0; i < nCol; i++) {
	max[i] = array[0][i];
    }
    for (int i = 0; i < nCol; i++) {
	for (int j = 1; j < nRow; j++) {
	    if (array[j][i] > max[i]) {
		max[i] = array[j][i];
	    }
	}
    }
    for (int i = 0; i < nCol; i++) {
	max[i] = max[i] + 1;
    }
    return 0;
}

/* max over each row */
int maxRow(double *array, int nRow, int nCol, double *max)
{
    for (int i = 0; i < nRow; i++) {
	max[i] = array[nCol * i];
    }
    for (int i = 0; i < nRow; i++) {
	for (int j = 1; j < nCol; j++) {
	    if (array[nCol * i + j] > max[i]) {
		max[i] = array[nCol * i + j];
	    }
	}
    }
    return 0;
}

/* The following function returns the dimensionality and number of samples from a file */
int retnnd(char *file, int *nd)
{
    FILE *f;
    char *tmp = new char[100];
    snprintf(tmp, 100, "%s %s %s", "wc -l ", file, " | grep -E -o '[0-9]+'");
    f = (FILE *) (popen(tmp, "r"));
    fgets(tmp, sizeof tmp, f);
    nd[0] = atoi(tmp);
    pclose(f);
    snprintf(tmp, 100, "%s %s %s", "head -n1", file, " | tr ' ' '\n' | wc -l | grep -E -o '[0-9]+'");
    f = (FILE *) (popen(tmp, "r"));
    fgets(tmp, sizeof tmp, f);
    nd[1] = atoi(tmp);
    pclose(f);
    delete tmp;
    return 1;
}

/* read nominal values in memory */
int readDataFile(char *file, int **mat)
{
    char buffer[1024];
    char *record, *line;
    int i = 0, j = 0;
    FILE *f = fopen(file, "r");
    if (f == NULL) {
	printf("\n file opening failed ");
	return -1;
    }
    while ((line = fgets(buffer, sizeof(buffer), f)) != NULL) {
	record = strtok(line, " \n");
	while (record != NULL) {
	    mat[i][j] = atoi(record);
	    j++;
	    record = strtok(NULL, " \n");
	}
	j = 0;
	i++;
    }
    return 0;
}

/* convert each column to presence absence vector */
int mat2bin(int **mat, double **nFeatSam, int n, int d)
{
    for (int i = 0; i < d; i++) {
	for (int j = 0; j < n; j++) {
	    nFeatSam[i][n * mat[j][i] + j] = 1;
	}
    }
    return 0;
}

/* compute number of categories in each column */
int compCat(int **array, int nRow, int nCol, int *max)
{
    maxCol(array, nRow, nCol, max);
    for (int i = 0; i < nCol; i++) {
	max[i] = max[i] + 1;
    }
    return 0;
}

/* sum over columns of a matrix */
int colSum(double *mat, int nRow, int nCol, double *sum)
{
    for (int i = 0; i < nCol; i++) {
	sum[i] = 0;
	for (int j = 0; j < nRow; j++) {
	    sum[i] = sum[i] + mat[nCol * j + i];
	}
    }
    return 0;
}

/* normalize each column */
int colDiv(double *mat, int nRow, int nCol)
{
    double *sum = new double[nCol];
    colSum(mat, nRow, nCol, sum);
    for (int i = 0; i < nCol; i++) {
	for (int j = 0; j < nRow; j++) {
	    mat[nCol * j + i] = mat[nCol * j + i] / sum[i];
	}
    }
    delete[] sum;
    return 0;
}

/* perform operation digamma(a_j) - digamma(sum a_j) */
int colSub(double *mat, int nRow, int nCol, double *sum)
{
    for (int i = 0; i < nCol; i++) {
	sum[i] = digamma(sum[i]);
	//cout << sum[i] << endl;
	for (int j = 0; j < nRow; j++) {
	    mat[nCol * j + i] = mat[nCol * j + i] - sum[i];
	}
    }
    return 0;
}

/* entropy of Dirichlet distribution
      sum(gammaln(alpha)) - gammaln(sum(alpha)) ...
    + (sum(alpha) - size(alpha, 1)) .* psi(sum(alpha)) ...
    - sum((alpha-1) .* psi(alpha));*/
int entDir(double *alpha, int nRow, int nCol, double *ent)
{
    double *sumAlpha = new double[nCol];
    colSum(alpha, nRow, nCol, sumAlpha);
    for (int j = 0; j < nCol; j++) {
	ent[j] = 0;
	for (int i = 0; i < nRow; i++) {
	    ent[j] = ent[j] + lgamma(alpha[nCol * i + j]) - (alpha[nCol * i + j] - 1) * digamma(alpha[nCol * i + j]);
	}
	ent[j] = ent[j] - lgamma(sumAlpha[j]) + (sumAlpha[j] - nRow) * digamma(sumAlpha[j]);
    }
    delete[] sumAlpha;
    return 0;
}

/* print a matrix */
int printMatrix(double *mat, int nRow, int nCol)
{
    for (int i = 0; i < nRow; i++) {
	for (int j = 0; j < nCol; j++) {
	    //cout << i << " " << j << endl;
	    cout << mat[nCol * i + j] << " ";
	}
	cout << endl;
    }
    return 0;
}

int printMatrix(double *mat, int nRow, int nCol, int nRow_disp, int nCol_disp)
{
    for (int i = 0; i < min(nRow_disp, nRow); i++) {
	for (int j = 0; j < min(nCol_disp, nCol); j++) {
	    //cout << i << " " << j << endl;
	    cout << mat[nCol * i + j] << " ";
	}
	cout << endl;
    }
    return 0;
}

int printMatrix(int *mat, int nRow, int nCol)
{
    for (int i = 0; i < nRow; i++) {
	for (int j = 0; j < nCol; j++) {
	    //cout << i << " " << j << endl;
	    cout << mat[nCol * i + j] << " ";
	}
	cout << endl;
    }
    return 0;
}

int printMatrix(int *mat, int nRow, int nCol, int nRow_disp, int nCol_disp)
{
    for (int i = 0; i < min(nRow_disp, nRow); i++) {
	for (int j = 0; j < min(nCol_disp, nCol); j++) {
	    //cout << i << " " << j << endl;
	    cout << mat[nCol * i + j] << " ";
	}
	cout << endl;
    }
    return 0;
}

/* write matrix to a file */
int writeMatrix(double *mat, int nRow, int nCol, char *file)
{
    ofstream f;
    f.open(file);
    for (int i = 0; i < nRow; i++) {
	for (int j = 0; j < nCol; j++) {
	    //cout << i << " " << j << endl;
	    f << mat[nCol * i + j] << " ";
	}
	f << endl;
    }
    f.close();
    return 0;
}

/* indRowCol = 0: print selected rows, = 1: print selected columns */
/* index is a binary matrix for which 1 values are printed */
int writeMatrix(double *mat, int nRow, int nCol, char *file, int indRowCol, int *index)
{
    ofstream f;
    f.open(file);
    for (int i = 0; i < nRow; i++) {
	for (int j = 0; j < nCol; j++) {
	    if (indRowCol == 0) {
		if (index[i] == 1) {
		    f << mat[nCol * i + j] << " ";
		}
	    }
	    if (indRowCol == 1) {
		if (index[j] == 1) {
		    f << mat[nCol * i + j] << " ";
		}
	    }
	}
	if (indRowCol == 0) {
	    if (index[i] == 1) {
		f << endl;
	    }
	}
	if (indRowCol == 1) {
	    f << endl;
	}
    }
    f.close();
    return 0;
}

/* parse optional arguments */
int parse(int argc, char **argv, double *parameters)
{
    for (int i = 1; i < argc; i = i + 2) {
	if (strcmp(argv[i], "--k") == 0) {
	    parameters[0] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--maxIter") == 0) {
	    parameters[1] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--thresh") == 0) {
	    parameters[2] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--priorMatLatSam") == 0) {
	    parameters[3] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--priorMatSamLat") == 0) {
	    parameters[4] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--randSeed") == 0) {
	    parameters[5] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--d") == 0) {
	    parameters[6] = atof(argv[i + 1]);
	}
	if (strcmp(argv[i], "--n") == 0) {
	    parameters[7] = atof(argv[i + 1]);
	}
    }
    return 1;
}

int parse2(int argc, char **argv, char *file, char *opfile)
{
    for (int i = 1; i < argc; i = i + 2) {
	if (strcmp(argv[i], "--input") == 0) {
	    strcpy(file, argv[i + 1]);
	}
	if (strcmp(argv[i], "--output") == 0) {
	    strcpy(opfile, argv[i + 1]);
	}
    }
    return 1;
}

/* find non-unique rows, assignment holds the first index matching a non-unique row */
int unique(int **mat, int nRow, int nCol, int *assignment)
{
    int match = 1;
    for (int i = nRow - 1; i >= 0; i--) {
	for (int j = 0; j < i; j++) {
	    match = 1;
	    for (int k = 0; k < nCol; k++) {
		match = match * int (mat[i][k] == mat[j][k]);
		if (match == 0) {
		    break;
		}
	    }
	    if (match == 0) {
		assignment[i] = i;
	    } else {
		assignment[i] = j;
		break;
	    }
	}
    }
    return 1;
}
