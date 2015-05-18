#include <cstdio>
#include <cstdlib>
#include <vector>
#include <limits>
#include <iostream>
using std::vector;
using std::numeric_limits;
using std::cout;
using std::endl;

#define BLAS_UNDERSCORE

#ifdef BLAS_UNDERSCORE
#define BLAS(name) name ## _
#else
#define BLAS(name) name
#endif

extern "C" {

int BLAS(dlaneg)( int* n, double* D, double* LLD, double* sigma, double* tol, int* r );

}

/*
dlaneg.c
A C implementation of dlaneg, as found in the LAPACK Library

Computes the Strum count encountered while factoring tridiagonal T - sigma I = L D L^T
This is equivalent to the number of eigenvalues of T less than sigma

Performs dtwqds transform, counts negative diagonal pivots

BJN 5/18/15

INPUTS:
	n (int) - size of matrix
	D (double *) - size n array of diagonal entries in D matrix in LDL^T factorization
	L (double *) - size (n-1) array of lower diagonal entries in L matrix in LDL^T factorization
	sigma (double) - shift of T matrix
	r (int) - which twisted factorization to use
OUTPUTS:
	negcnt (int) - strum count of T - sigma I

*/

// generate NaN
template<typename Real>
Real generate_nan(){
	if (numeric_limits<Real>::has_quiet_NaN) {
		return numeric_limits<Real>::quiet_NaN();
	} else {
		return Real(0)/Real(0);
	}
}

template<typename Real>
int laneg(int n, vector<Real>& D, vector<Real>& L, Real sigma, int r) {
	int negcnt = 0;

	int neg1, neg2;
	Real dminus, dplus, gamma, p, t, tmp;

	// upper part of twisted factorizationi LDL^T - sigma I = L+D+L+^T
	neg1 = 0;
	t = -sigma;
	for (int j = 0; j < r-1; j++) {
		dplus = D[j] + t;
		if (dplus < 0) {
			neg1++;
		}
		tmp = t / dplus;
		t = tmp * L[j] * L[j] * D[j] - sigma;
	}
	negcnt += neg1;

	// lower poart of twisted factorization LDL^T - sigma I = U-D-U-^T
	neg2 = 0;
	p = D[n-1] - sigma;
	for ( int j = n-1; j >= r; j--) {
		dminus = L[j] * L[j] * D[j] + p;
		if (dminus < 0) {
			neg2++;
		}
		tmp = p / dminus;
		p = tmp * D[j] - sigma;
	}
	negcnt += neg2;

	// twist in factorization
	gamma = (t + sigma) + p;
	if (gamma < 0) {
		negcnt++;
	}

	return negcnt;
}




int main(int argc, char *argv[]) {
	
	if (argc != 4) {
		printf("argc = %d\n", argc);
		printf("please supply n, sigma, and r\n");
		return 1;
	}
	int n = atoi(argv[1]);
	double sigma = atof(argv[2]);
	int r = atoi(argv[3]);

	cout << generate_nan<double>() << endl;
	cout << sizeof(generate_nan<double>()) << endl;
	cout << generate_nan<float>() << endl;
	cout << sizeof(generate_nan<float>()) << endl;

	printf("---------------------\n");
	printf("testing dlaneg.c\n");
	int negcnt;
	//int n = 3;
	//double sigma = 3.0;
	//int r = 3;
	//double * L = (double *) malloc((n-1) * sizeof(double));
	//double * D = (double *) malloc(n * sizeof(double));
	vector<double> L(n-1);
	vector<double> D(n);
	for (int i = 0; i < n-1; i++) {
		L[i] = 1.0;
	}
	for (int i = 0; i < n; i++) {
		D[i] = 1.0;
	}

	negcnt = laneg(n, D, L, sigma, r);
	printf("n = %d, sigma = %.3f, r = %d \n", n, sigma, r);
	cout << "my routine" << endl;
	printf("Sturm count is %d\n", negcnt);

	for (int i = 0; i < n-1; i++) {
		L[i] = 1.0;
	}
	for (int i = 0; i < n; i++) {
		D[i] = 1.0;
	}
	double pivmin = 1e-16;
	for (int i = 0; i < n-1; i++) {
		L[i] = L[i]*L[i]*D[i];
	}
	negcnt = BLAS(dlaneg)(&n, D.data(), L.data(), &sigma, &pivmin, &r);
	printf("n = %d, sigma = %.3f, r = %d \n", n, sigma, r);
	cout << "using LAPACK" << endl;
	printf("Sturm count is %d\n", negcnt);

	//free(L);
	//free(D);
	// vectors are freed at end of scope
	printf("---------------------\n");
	return 0;
}