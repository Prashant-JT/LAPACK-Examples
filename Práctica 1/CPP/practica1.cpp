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
			} else {
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
		lda = m; //n si fuese por columnas. (cuadrada)
		ipiv = new lapack_int[n];
		matrixC = matrix;
		matrixOver = utils.clone(matrixC, m * n);
	}

	void LU() {
		lapack_int res;		
		res = LAPACKE_sgetrf(layout, m, n, matrixOver, lda, ipiv);
		utils.printMatrix(matrixOver, m, n);
	}

	void determinante() {
		float det = 1;
		for (int i = 0; i < m; i++) {
			det *= matrixOver[m*i+i];
		}
		printf("\t\tDeterminante = %f", det);
	}

	void inversa() {
		//NO SE COMO RESOLVER AX=I
	}

	void inversa_sgetri() {
		lapack_int res;
		res = LAPACKE_sgetri(layout, n, matrixOver, lda, ipiv);
		utils.printMatrix(matrixOver, m, n);
	}

private:
	lapack_int m, n, lda, * ipiv, layout;
	float* matrixC;
	float* matrixOver;
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

	float* matSxS = utils.processMatrix(mat, m*n);
	printf("Matriz a operar:\n");
	utils.printMatrix(matSxS, m, n);

	funcion2 funcion2(matSxS, m, n);
	printf("\n\nEjercicio 2_A (LU):\n");	
	funcion2.LU();

	printf("\n\nEjercicio 2_B (Determinante):\n");
	funcion2.determinante();

	printf("\n\nEjercicio 2_C (Inversa):\n");
	funcion2.inversa();

	printf("\n\nEjercicio 2_D (Inversa con SGETRI):\n");
	funcion2.inversa_sgetri();

	//printf("\n\nEjercicio 3:\n");
	//funcion3();

	char a = std::getchar();
	return 0;
}