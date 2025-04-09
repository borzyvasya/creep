#include <iostream>
#include <iomanip>
#include <omp.h>
#include "mkl.h"
#include <fstream>

const int N = 10;

using namespace std;

void generateMatrix(double[][N], double[]);
void printMatrix(double[][N], double[]);
void checkDiagonalDominance(double[][N]);
void simpleIteration(double[][N], double[], double[], double);
void checkSolution(double[][N], double[], double[]);
void gaussSeidel(double[][N], double[], double[], double);
void sor(double[][N], double[], double[], double, double);
void outputSolutions(double[][N], double[], double[], string);

int main() {
    double A[N][N], b[N], x[N];
    double eps = 1e-6;
    generateMatrix(A, b);

    cout << fixed << setprecision(3) << "Matrix A and b: " << endl;
    printMatrix(A, b);

    cout << "\nSimple iterations (count of iterations): " << endl;
    simpleIteration(A, b, x, eps);
    outputSolutions(A, b, x, "Simple iterations");

    cout << "\nGauss-Seidel (count of iterations): " << endl;
    gaussSeidel(A, b, x, eps);
    outputSolutions(A, b, x, "Gauss-Seidel");

    cout << "\nSuccessive Over-Relaxation (SOR) (count of iterations): " << endl;
    sor(A, b, x, eps, 1.3); // omega = 1.3 как пример
    outputSolutions(A, b, x, "Successive Over-relaxation");

    return EXIT_SUCCESS;
}

void generateMatrix(double A[][N], double b[]) {
    double fixedA[N][N] = {
        {10, 1, 2, 0, 0, 0, 0, 0, 0, 0},
        {1, 10, 1, 2, 0, 0, 0, 0, 0, 0},
        {2, 1, 10, 1, 2, 0, 0, 0, 0, 0},
        {0, 2, 1, 10, 1, 2, 0, 0, 0, 0},
        {0, 0, 2, 1, 10, 1, 2, 0, 0, 0},
        {0, 0, 0, 2, 1, 10, 1, 2, 0, 0},
        {0, 0, 0, 0, 2, 1, 10, 1, 2, 0},
        {0, 0, 0, 0, 0, 2, 1, 10, 1, 2},
        {0, 0, 0, 0, 0, 0, 2, 1, 10, 1},
        {0, 0, 0, 0, 0, 0, 0, 2, 1, 10}
    };

    double fixedB[N] = { 13, 16, 20, 22, 24, 24, 23, 22, 19, 16 };

    for (int i = 0; i < N; i++) {
        b[i] = fixedB[i];
        for (int j = 0; j < N; j++) {
            A[i][j] = fixedA[i][j];
        }
    }
    //srand(time(0));
    //double CF = 1.0e04;

    //for (int i = 0; i < N; i++) {
    //    double sum = 0;
    //    for (int j = 0; j < N; j++) {
    //        if (i != j) {
    //            A[i][j] = rand() / CF;
    //            sum += fabs(A[i][j]);
    //        }
    //    }
    //    A[i][i] = sum + (rand() / CF);
    //}

    //for (int i = 0; i < N; i++) b[i] = rand() / CF;
}

void printMatrix(double A[][N], double b[]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << "\t";
        cout << "| " << b[i] << endl;
    }
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
        //#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            // Используем cblas_ddot для вычисления A[i] * x
            double sum = cblas_ddot(N, &A[i][0], 1, x, 1);
            // Вычитаем диагональный вклад A[i][i] * x[i]
            x_new[i] = (b[i] - (sum - A[i][i] * x[i])) / A[i][i];
        }

        error = 0;
        //#pragma omp parallel for reduction(max:error)
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

void gaussSeidel(double A[][N], double b[], double x[], double epsilon) {
    int iter = 0;
    double error = epsilon + 1;

    for (int i = 0; i < N; i++)
        x[i] = 0;

    while (error > epsilon) {
        error = 0;
        #pragma omp parallel for reduction(max:error)
        for (int i = 0; i < N; i++) {
            double x_old = x[i];
            // Вычисляем A[i] * x с использованием cblas_ddot
            double sum = cblas_ddot(N, &A[i][0], 1, x, 1);

            // Обновляем x[i], исключая старый диагональный вклад
            x[i] = (b[i] - (sum - A[i][i] * x_old)) / A[i][i];
            error = max(error, fabs(x[i] - x_old));
        }

        iter++;
        cout << setprecision(7) << "Iteration " << iter << ", error: " << error << endl;
    }
}

void sor(double A[][N], double b[], double x[], double epsilon, double omega) {
    int iter = 0;
    double error = epsilon + 1;

    for (int i = 0; i < N; i++)
        x[i] = 0;

    while (error > epsilon) {
        error = 0;
        for (int i = 0; i < N; i++) {
            double x_old = x[i];
            double sum = cblas_ddot(N, &A[i][0], 1, x, 1);
            double gs_update = (b[i] - (sum - A[i][i] * x_old)) / A[i][i]; // Гаусс-Зейдель шаг
            x[i] = (1 - omega) * x_old + omega * gs_update; // Применение релаксации
            error = max(error, fabs(x[i] - x_old));
        }

        iter++;
        cout << setprecision(7) << "Iteration " << iter << ", error: " << error << endl;
    }
}

void outputSolutions(double A[][N], double b[], double x[], string nameOfAlgorithm) {
    cout << "\nSolution(s) " << "(" << nameOfAlgorithm << "):" << endl;
    for (int i = 0; i < N; i++)
        cout << setprecision(7) << "x[" << i << "] = " << x[i] << endl;

    checkSolution(A, b, x);
}


