#include <iostream>
#include <omp.h>
#include "mkl.h"
#include <iomanip>

const int N = 10;

using namespace std;

void generateMatrix(double[][N], double[]);
void checkDiagonalDominance(double[][N]);
// Метод простых итераций с oneMKL (cblas_ddot)
void simpleIteration(double[][N], double[], double[], double);
// Проверка решения с oneMKL (cblas_dgemv)
void checkSolution(double[][N], double[], double[]);


int main() {
    double A[N][N], b[N], x[N];

    generateMatrix(A, b);

    cout << fixed << setprecision(3) << "Matrix A and b: " << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << "\t";
        cout << "| " << b[i] << endl;
    }

    simpleIteration(A, b, x, 1e-6);

    cout << "\nSolution: " << endl;
    for (int i = 0; i < N; i++)
        cout << setprecision(7) << "x[" << i << "] = " << x[i] << endl;

    return EXIT_SUCCESS;
}

void generateMatrix(double A[][N], double b[]) {
    srand(time(0));
    double CF = 1.0e04;

    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                A[i][j] = rand() / CF;
                sum += fabs(A[i][j]);
            }
        }
        A[i][i] = sum + (rand() / CF);
    }

    for (int i = 0; i < N; i++) b[i] = rand() / CF;
}

void checkDiagonalDominance(double A[][N]) {
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) 
                sum += fabs(A[i][j]);
        }
        if (fabs(A[i][i]) <= sum) {
            cout << "Matrix does not satisfy diagonal dominance." << endl;
            exit(1);
        }
    }
}

void simpleIteration(double A[][N], double b[], double x[], double epsilon) {
    double x_new[N];
    int iter = 0;
    double error = epsilon + 1;

    for (int i = 0; i < N; i++)
        x[i] = 0;

    while (error > epsilon) {
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            // Используем cblas_ddot для вычисления A[i] * x
            double sum = cblas_ddot(N, &A[i][0], 1, x, 1);
            // Вычитаем диагональный вклад A[i][i] * x[i]
            x_new[i] = (b[i] - (sum - A[i][i] * x[i])) / A[i][i];
        }

        error = 0;
        #pragma omp parallel for reduction(max:error)
        for (int i = 0; i < N; i++) {
            error = max(error, fabs(x_new[i] - x[i]));
            x[i] = x_new[i];
        }

        iter++;
        cout << setprecision(7) << "Iteration " << iter << ", error: " << error << endl;
    }
}

void checkSolution(double A[][N], double b[], double x[]) {
    double Ax[N];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, &A[0][0], N, x, 1, 0.0, Ax, 1);

    double residual = 0;
    for (int i = 0; i < N; i++)
        residual = max(residual, fabs(Ax[i] - b[i]));
 
    cout << "Residual: " << residual << endl;
}