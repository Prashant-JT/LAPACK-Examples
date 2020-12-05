#include <cstdio>
#include <cstdlib>
#include <string>

#include <mkl.h>

using namespace std;

class Utils {
public:
	void printMatrix(float* mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (mat[i] < 0) {
				printf("%4.5f\t", mat[i]);
			}
			else {
				printf("%4.6f\t", mat[i]);
			}
			if ((i + 1) % n == 0) {
				printf("\n\t\t");
			}
		}
		printf("\n");
	}

	float* clone(float* matrix, int size) {
		float* res = new float[size];
		memcpy(res, matrix, sizeof(float) * size);
		return res;
	}

	float* processMatrix(string mat, int size) {
		float* vect = new float[size];
		int p = 0;
		string temp = "";

		for (int i = 0; i < mat.size(); i++) {
			if (mat[i] != ' ') {
				temp += mat[i];
			}
			else {
				vect[p] = stof(temp);
				p++;
				temp = "";
			}
		}
		return vect;
	}
};

class funcion2 {
public:
	funcion2(float* matrix, int m_, int n_) {
		layout = LAPACK_ROW_MAJOR;
		m = m_;
		n = n_;
		lda = m;
		matrixC = matrix;		
	}
	
	void cholesky(char uplo) {
		float* matrixOver = cloneMatrix();
		lapack_int res;
		res = LAPACKE_spotrf(layout, uplo, n, matrixOver, lda);
		utils.printMatrix(matrixOver, m, n);
	}

	void QR() {
		float* matrixOver = cloneMatrix();
		lapack_int res;
		float * tau = new float[n];
		res = LAPACKE_sgeqrf(layout, m, n, matrixOver, lda, tau);
		utils.printMatrix(matrixOver, m, n);
	}

	void SVD(char jobu, char jobvt) {
		float* matrixOver = cloneMatrix();	
		float* s = new float[m];
		float* superb = new float[n];
		lapack_int res, ldvt, ldu;
		ldvt = lda;
		ldu = lda;
		float* u = new float[ldu];
		float* vt = new float[ldvt];
		res = LAPACKE_sgesvd(layout, jobu, jobvt, m, n, matrixOver, lda, s, u, ldu, vt, ldvt, superb);
		utils.printMatrix(matrixOver, m, n);
	}

	void autos(char jobz, char uplo) {
		float* matrixOver = cloneMatrix();
		float* w = new float[m];
		lapack_int res;
		res = LAPACKE_ssyev(layout, jobz, uplo, n, matrixOver, lda, w);
		printf("\t  Autovectores:\n");
		utils.printMatrix(matrixOver, m, n);
		printf("\t  Autovalores:\n");
		utils.printMatrix(w, 1, m);
	}
private:
	lapack_int m, n, lda, layout;
	float* matrixC;	
	Utils utils;

	float* cloneMatrix() {
		return utils.clone(matrixC, m * n);
	}
};

int main(int argc, char* argv[]) {
	Utils utils;
	int m = 6; int n = 6;
	string matSym =
		"43 59 51 50 39 54 "
		"59 118 103 101 61 97 "
		"51 103 121 93 56 102 "
		"50 101 93 111 69 99 "
		"39 61 56 69 77 86 "
		"54 97 102 99 86 134 ";

	float* matSxS = utils.processMatrix(matSym, m * n);
	printf("Matriz a operar:\n");
	utils.printMatrix(matSxS, m, n);

	funcion2 funcion2(matSxS, m, n);
	printf("\n\nEjercicio 3_A (CHOLESKY) --> L es traspuesta de U:\n");
	printf("\tLower (Triangulo inferior se sobreescribe):\n");
	funcion2.cholesky('L');
	printf("\tUpper (Triangulo superior se sobreescribe):\n");
	funcion2.cholesky('U');

	printf("\n\nEjercicio 3_B (QR):\n");
	funcion2.QR();

	/*
	JOBU is CHARACTER*1
          Specifies options for computing all or part of the matrix U:
          = 'A':  all M columns of U are returned in array U:
          = 'S':  the first min(m,n) columns of U (the left singular
                  vectors) are returned in the array U;
          = 'O':  the first min(m,n) columns of U (the left singular
                  vectors) are overwritten on the array A;
          = 'N':  no columns of U (no left singular vectors) are
                  computed.
	
	JOBVT is CHARACTER*1
		  Specifies options for computing all or part of the matrix
		  V**T:
		  = 'A':  all N rows of V**T are returned in the array VT;
		  = 'S':  the first min(m,n) rows of V**T (the right singular
				  vectors) are returned in the array VT;
		  = 'O':  the first min(m,n) rows of V**T (the right singular
				  vectors) are overwritten on the array A;
		  = 'N':  no rows of V**T (no right singular vectors) are
				  computed.

		  JOBVT and JOBU cannot both be 'O'.

	funcion2.SVD(jobu, jobvt); Defecto ambos igual a 'A'.
	Mucho texto: 
	http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html
	*/

	printf("\n\nEjercicio 3_C (SVD):\n");
	funcion2.SVD('A', 'A');

	/*
	JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.

	UPLO is CHARACTER*1
		  = 'U':  Upper triangle of A is stored;
		  = 'L':  Lower triangle of A is stored.

	Se pide ambos, por tanto se pasa jobz='V', uplo=L o U (traspuestas).
	Como jobz='V', la matriz contiene los autovectores.
	Los autovalores se sitúan en el vector "w".
	*/

	printf("\n\nEjercicio 3_D (Autovalores y Autovectores):\n");
	printf("\tLower (Triangulo inferior se sobreescribe):\n");
	funcion2.autos('V', 'L');
	printf("\tUpper (Triangulo superior se sobreescribe):\n");
	funcion2.autos('V', 'U');

	char a = std::getchar();
	return 0;
}
