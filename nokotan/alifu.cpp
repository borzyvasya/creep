#include "mkl.h"
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

const int M = 10, N = 10, LDA = M, LDB = N, NRHS = N;

void print_mat(string, MKL_INT, MKL_INT, double[], MKL_INT);
void mult(double[], MKL_INT, double[], MKL_INT, double[]);
void create_mat(double[], int, int);

int main() {
    MKL_INT m = M, n = N, lda = LDA, ldb = LDB, info,
        nrhs = NRHS, ipiv[N];

    double b[LDB * NRHS], a[LDA * M], asv[LDB * M], bsv[LDB * NRHS];

    
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            b[i * N + j] = (i == j) ? 1.0 : 0.0; 
        }
    }

    create_mat(a, m, n); // A

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < LDA; ++j) {
            asv[i * LDA + j] = a[i * LDA + j];
        }
    }

    
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            bsv[i * NRHS + j] = b[i * NRHS + j];
        }
    }

    print_mat("Matrix [A]", m, n, asv, lda);

    // Вычисление обратной матрицы с помощью LAPACKE_dgesv
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);

    if (info > 0) {
        cout << "Computing failed cuz of conversion..." << endl;
        return EXIT_FAILURE; 
    }

    if (n <= 20) {
        print_mat("Matrix [LU] : ", m, n, a, ldb);
        print_mat("Inverse of Matrix [A]", m, nrhs, b, ldb);

        mult(asv, lda, b, ldb, bsv); // Проверка A * A^-1 = E

        print_mat("Matrix [A] applying to its inverse form", m, nrhs, bsv, ldb);

        cout << endl << "Success! press any kek to end..." << endl;
        cin.get();
        return EXIT_SUCCESS;
    }

    return EXIT_SUCCESS;
}

void print_mat(string str, MKL_INT rows, MKL_INT cols, double matrix[], MKL_INT ld) {
    cout << str << endl;
    cout << fixed << setprecision(3);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = matrix[i * ld + j];
            if (abs(value) < 1e-04) value = 0.0;
            cout << setw(10) << value;
        }
        cout << endl;
    }
}

void mult(double a[], MKL_INT lda, double b[], MKL_INT ldb, double result[]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k) 
                sum += a[i * lda + k] * b[k * ldb + j];
  
            result[i * ldb + j] = sum;
        }
    }
}

void create_mat(double matrix[], int rows, int cols) {
    srand(time(NULL));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (i == j) matrix[i * cols + j] = rand() / 1.0e+4; 
            else matrix[i * cols + j] = 1.0; 
            
        }
    }
}