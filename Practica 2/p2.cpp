#include <cstdio>
#include <cstdlib>
#include <string>
#include<random>

#include <mkl.h>

using namespace std;
std::default_random_engine generator;
std::normal_distribution<double> random(0.0, 1.0);

class Utils {

public:
	void printMatrix(double* mat, int m, int n) {
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

	double* clone(double* matrix, int size) {
		double* res = new double[size];
		memcpy(res, matrix, sizeof(double) * size);
		return res;
	}

	double* processMatrix(string mat, int size) {
		double* vect = new double[size];
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

class funcion1 {

public:
	funcion1(double* matrix, int m_, int n_) {
		layout = LAPACK_ROW_MAJOR;
		nhrs = 1;
		n = n_;
		lda = n; 
		ldb = m_;
		matrixC = matrix;
		matrixB = new double[ldb * nhrs];

		for (int i = 0; i < ldb * nhrs; i++) {
			matrixB[i] = 0;
		}
	}

	void LU() {
		lapack_int res;
		double *matrixOver = utils.clone(matrixC, n * n);
		lapack_int *ipiv = new lapack_int[n];
		res = LAPACKE_dgesv(layout, n, nhrs, matrixOver, lda, ipiv, matrixB, ldb);
		utils.printMatrix(matrixOver, n, n);
	}

	void inversa() {
		double *matrixOver = utils.clone(matrixC, n * n);
		lapack_int* ipiv = new lapack_int[n];
		double *matrixI = new double[n * n];

		for (int i = 0; i < n * n; i++) {
			matrixI[i] = 0.0;
		}

		for (int i = 0; i < n; i++) {
			matrixI[n*i+i] = 1.0;
		}

		lapack_int res;
		nhrs = 6;
		res = LAPACKE_dgesv(layout, n, nhrs, matrixOver, lda, ipiv, matrixI, ldb);
		utils.printMatrix(matrixI, n, n);
	}

private:
	lapack_int nhrs, n, lda, ldb, layout;
	double* matrixC;
	double* matrixB;
	Utils utils;
};

class funcion2 {

public:
	funcion2(int m_, int n_) {
		layout = LAPACK_ROW_MAJOR;
		nhrs = 1;
		lda = n_;
		ldb = m_;
		m = m_;
		n = n_;
		matrixA = new double[m * n];
		matrixB = new double[m * n];
		matrixZeros = new double[m * n];

		for (int i = 0; i < m*n; i++) {
			matrixA[i] = 0.0;
			matrixB[i] = random(generator);
			matrixZeros[i] = 0.0;
		}

	}

	void generateMatrix() {

		for (int i = 0; i < n; i++) {
			matrixA[n * i + i] = random(generator);

			if (i < n) {
				matrixA[n * i + i + 1] = random(generator);
			}
			
			if (i > 0) {
				matrixA[n * i + i - 1] = random(generator);
			}
			
		}
		
		printf("\n Matriz A (tridiagonal): \n");
		utils.printMatrix(matrixA, m, n);
		printf("\n Matriz B (general): \n");
		utils.printMatrix(matrixB, m, n);
	}
	
	void matrixCompact() {
		compactA = new double[5 * n];

		for (int i = 0; i < 5 * n; i++) {
			compactA[i] = 0;
		}

		int count = 1;
		for (int i = n + 1; i < 2 * n; i++) {
			compactA[i] = matrixA[count];
			count += n + 1;
		}

		count = 0;
		for (int i = 2 * n; i <= 3 * n; i++) {
			compactA[i] = matrixA[count];
			count += n + 1;
		}

		count = n;
		for (int i = 3 * n; i < 4 * n; i++) {
			compactA[i] = matrixA[count];
			count += n + 1;
		}

		utils.printMatrix(compactA, n, n);
	}

	void dgesvGeneral() {
		lapack_int* ipiv;
		double* matrixOver;
		double* matrixZ;
		int res;

		double start, fin = dsecnd();
		start = dsecnd();
		for(int i = 0; i < 1000000; i++) {
			matrixOver = utils.clone(matrixB, m * n);
			matrixZ = utils.clone(matrixZeros, m * n);
			ipiv = new lapack_int[n];
			res = LAPACKE_dgesv(layout, n, nhrs, matrixOver, lda, ipiv, matrixZ, ldb);
		}
		fin = (dsecnd() - start);
		double totalTime = fin / 1000000;

		printf("\n Matriz general (INFO = %d)", res);
		printf(" | Tiempo promedio en segundos = %f\n", totalTime);

		res = LAPACKE_dgesv(layout, n, nhrs, matrixOver, lda, ipiv, matrixZ, ldb);
		utils.printMatrix(matrixOver, n, n);
	}

	void dgesvBanda() {
		lapack_int* ipiv;
		double* matrixOver;
		double* matrixZ;
		int res;
		int ldab = 4;

		double start, fin = dsecnd();
		start = dsecnd();
		for (int i = 0; i < 1000000; i++) {
			matrixOver = utils.clone(compactA, m * n);
			matrixZ = utils.clone(matrixZeros, m * n);
			ipiv = new lapack_int[n];
			res = LAPACKE_dgbsv(layout, n, 1, 1, nhrs, matrixOver, ldab, ipiv, matrixZ, ldb);	
		}
		fin = (dsecnd() - start);
		double totalTime = fin / 1000000;

		printf("\n Matriz tridiagonal (INFO = %d)", res);
		printf(" | Tiempo promedio en segundos = %f\n", totalTime);
		utils.printMatrix(matrixOver, n, n);
	}

private:
	lapack_int nhrs, m, n, lda, ldb, layout;
	double* matrixA;
	double* matrixB;
	double* matrixZeros;
	double* compactA;
	Utils utils;
};

int main(int argc, char* argv[]) {
	Utils utils;
	int m = 6; int n = 6;
	string mat =
		"9 3 5 3 7 2 "
		"6 10 2 6 2 1 "
		"10 5 9 5 9 8 "
		"5 3 6 8 2 8 "
		"4 5 4 1 7 5 "
		"5 8 3 9 6 1 ";

	double* matSxS = utils.processMatrix(mat, m * n);
	printf("Matriz a operar:\n");
	utils.printMatrix(matSxS, m, n);

	funcion1 funcion1(matSxS, m, n);
	printf("\n\nEjercicio 1 con dgesv (LU):\n");
	funcion1.LU();

	printf("\n\nEjercicio 1 con dgesv (inversa):\n");
	funcion1.inversa();

	m = 4; n = 4;
	funcion2 funcion2(m, n);
	printf("\n\nEjercicio 2_A (matrices tridiagonales):\n");
	funcion2.generateMatrix();

	printf("\n\nEjercicio 2_B (matriz A compacta):\n");
	funcion2.matrixCompact();

	printf("\n\nEjercicio 2_C y 2_D (matriz A (LAPACKE_dgesv)):\n");
	funcion2.dgesvGeneral();
	funcion2.dgesvBanda();

	char a = std::getchar();
	return 0;
}