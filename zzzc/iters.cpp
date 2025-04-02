#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include "mkl.h"
#include <iomanip>

const int N = 10;

using namespace std;

void generateMatrix(double A[N][N], double b[N]) {
    srand(time(0));
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                A[i][j] = rand() % 10 + 1;
                sum += fabs(A[i][j]);
            }
        }
        A[i][i] = sum + (rand() % 10 + 1);
    }
    for (int i = 0; i < N; i++) {
        b[i] = rand() % 50 + 1;
    }
}

bool checkDiagonalDominance(double A[N][N]) {
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                sum += fabs(A[i][j]);
            }
        }
        if (fabs(A[i][i]) <= sum) {
            return false;
        }
    }
    return true;
}

// Метод простых итераций с oneMKL (cblas_ddot)
void simpleIteration(double A[N][N], double b[N], double x[N], double epsilon) {
    double x_new[N];
    int iter = 0;
    double error = epsilon + 1;

    for (int i = 0; i < N; i++)
        x[i] = 0;

    if (!checkDiagonalDominance(A)) {
        cout << "Matrix does not satisfy diagonal dominance.\n";
        return;
    }

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
        cout << "Iteration " << iter << ", error: " << error << "\n";
    }
}

// Проверка решения с oneMKL (cblas_dgemv)
void checkSolution(double A[N][N], double b[N], double x[N]) {
    double Ax[N];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, &A[0][0], N, x, 1, 0.0, Ax, 1);

    double residual = 0;
    for (int i = 0; i < N; i++) {
        residual = max(residual, fabs(Ax[i] - b[i]));
    }
    cout << "Residual: " << residual << "\n";
}

int main() {
    double A[N][N], b[N], x[N];

    generateMatrix(A, b);

    cout << "Matrix A and b:\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << "\t";
        cout << "| " << b[i] << "\n";
    }

    simpleIteration(A, b, x, 1e-6);

    cout << "\nSolution:\n";
    for (int i = 0; i < N; i++)
        cout << setprecision(10) << fixed << "x[" << i << "] = " << x[i] << "\n";

    checkSolution(A, b, x);

    return 0;
}