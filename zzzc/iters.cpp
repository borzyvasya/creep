#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include "mkl.h" // Заголовочный файл oneMKL
#include <iomanip>

const int N = 10;

using namespace std;

void generateMatrix(double A[N][N], double b[N]) {
    std::srand(std::time(0));
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                A[i][j] = std::rand() % 10 + 1;
                sum += fabs(A[i][j]);
            }
        }
        A[i][i] = sum + (std::rand() % 10 + 1);
    }
    for (int i = 0; i < N; i++) {
        b[i] = std::rand() % 50 + 1;
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

// Метод простых итераций с oneMKL
void simpleIteration(double A[N][N], double b[N], double x[N], double epsilon) {
    double x_new[N];
    double temp[N]; // Временный массив для хранения Ax
    int iter = 0;
    double error = epsilon + 1;

    // Инициализация начального приближения
    for (int i = 0; i < N; i++)
        x[i] = 0;

    if (!checkDiagonalDominance(A)) {
        cout << "Matrix does not satisfy the diagonal dominance condition.\n";
        return;
    }

    while (error > epsilon) {
        // Вычисление Ax с помощью oneMKL
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, &A[0][0], N, x, 1, 0.0, temp, 1);

        // Вычисление нового приближения x_new
        #pragma omp parallel for
        for (int i = 0; i < N; i++) 
            x_new[i] = (b[i] - temp[i] + A[i][i] * x[i]) / A[i][i];

        // Вычисление погрешности
        error = 0;
        #pragma omp parallel for reduction(max:error)
        for (int i = 0; i < N; i++) {
            error = std::max(error, fabs(x_new[i] - x[i]));
            x[i] = x_new[i];
        }

        iter++;
        cout << "Iteration " << iter << ", error: " << error << "\n";
    }

    cout << "Solution found after " << iter << " iterations.\n";
}

// Проверка решения с использованием oneMKL
void checkSolution(double A[N][N], double b[N], double x[N]) {
    double Ax[N];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, &A[0][0], N, x, 1, 0.0, Ax, 1);

    double residual = 0.0;
    for (int i = 0; i < N; i++) {
        residual = std::max(residual, fabs(Ax[i] - b[i]));
    }
    cout << "Residual (max |Ax - b|): " << residual << "\n";
}

int main() {
    double A[N][N];
    double b[N];
    double x[N];

    generateMatrix(A, b);

    cout << "Generated matrix A:\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << A[i][j] << "\t";
        }
        cout << "| " << b[i] << "\n";
    }

    double epsilon = 1e-6;
    simpleIteration(A, b, x, epsilon);

    cout << "\nSolution of the system (Simple Iteration):\n";
    for (int i = 0; i < N; i++) {
        cout << setprecision(10) << fixed << "x[" << i << "] = " << x[i] << "\n";
    }

    checkSolution(A, b, x);


    return EXIT_SUCCESS;
}