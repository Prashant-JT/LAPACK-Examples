#include <cstdio>
#include <cstdlib>
#include <string>

#include <mkl.h>

using namespace std;

class Utils {
public:
	void printMatrix(double* mat, int m, int n) {
		printf("\t\t");
		for (int i = 0; i < m * n; i++)
		{
			if (mat[i] < 0) {
				printf("%4.5f\t", mat[i]);
			} else {
				printf("%4.6f\t", mat[i]);
			}			
			if ((i + 1) % n == 0) {
				printf("\n\t\t");
			}
		}
		printf("\n");
	}

	double * clone(double* matrix, int size) {
		double* res = new double[size];
		memcpy(res, matrix, sizeof(double) * size);
		return res;
	}

	lapack_int * cloneInt(lapack_int* matrix, int size) {
		lapack_int* res = new lapack_int[size];
		memcpy(res, matrix, sizeof(lapack_int) * size);
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

	double* getIdentity(int size) {
		double* matrix = new double[(int)(size)];
		for (int i = 0; i < size*size; i++) {
			matrix[i] = 0.0;
		}
		for (int i = 0; i < size; i++) {
			matrix[size * i + i] = 1.0;
		}
		return matrix;
	}
};

class funcion2 {
public:
	funcion2(double* matrix, int m_, int n_) {
		layout = LAPACK_ROW_MAJOR;
		m = m_;
		n = n_;
		lda = m; //n si fuese por columnas. (cuadrada)
		ipiv = new lapack_int[n];
		matrixC = matrix;
		matrixOver = utils.clone(matrixC, m * n);
	}

	void LU() {
		lapack_int res;		
		res = LAPACKE_dgetrf(layout, m, n, matrixOver, lda, ipiv);
		utils.printMatrix(matrixOver, m, n);
	}

	void determinante() {
		double det = 1;
		for (int i = 0; i < m; i++) {
			det *= matrixOver[m*i+i];
		}
		printf("\t\tDeterminante = %f", det);
	}

	void inversa(char trans) {
		lapack_int res, nrhs, * ipivB, ldb;
		nrhs = n;
		ldb = n;
		ipivB = utils.cloneInt(ipiv, n);
		double* matrixB = utils.getIdentity(ldb);	
		res = LAPACKE_dgetrs(layout, trans, n, nrhs, matrixOver, lda, ipivB, matrixB, ldb);
		utils.printMatrix(matrixB, m, n);
	}

	void inversa_dgetri() {
		lapack_int res;
		res = LAPACKE_dgetri(layout, n, matrixOver, lda, ipiv);
		utils.printMatrix(matrixOver, m, n);
	}

private:
	lapack_int m, n, lda, * ipiv, layout;
	double* matrixC;
	double* matrixOver;
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

	double* matSxS = utils.processMatrix(mat, m*n);
	printf("Matriz a operar:\n");
	utils.printMatrix(matSxS, m, n);

	funcion2 funcion2(matSxS, m, n);
	printf("\n\nEjercicio 2_A (LU):\n");	
	funcion2.LU();

	printf("\n\nEjercicio 2_B (Determinante):\n");
	funcion2.determinante();

	/*
	TRANS is CHARACTER*1
          Specifies the form of the system of equations:
          = 'N':  A * X = B  (No transpose)
          = 'T':  A**T* X = B  (Transpose)
          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
	*/
	printf("\n\nEjercicio 2_C (Inversa AX=I):\n");
	funcion2.inversa('N');

	printf("\n\nEjercicio 2_D (Inversa con DGETRI):\n");
	funcion2.inversa_dgetri();

	char a = std::getchar();
	return 0;
}