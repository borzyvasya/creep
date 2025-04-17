#include <iostream>
#include <iomanip>
#include "mkl.h"

using namespace std;

const int N = 10;
const double EPS   = 0.00001;
const double EPS1  = 0.0001;
const double OMEGA = 1.067;

void generate_matrix(double[][N], double[]);
void print_matrix(double[][N], double[]); 
void print_LDU(double[][N], double[], double[][N]);
void LDU_decomposition(double[][N], double[][N], double[], double[][N]);
void simple_iteration_LDU(double[][N], double[], double[], double[], double[][N], double[][N], double, double);
void gauss_seidel_LDU(double[][N], double[], double[], double[], double[][N], double[][N], double, double);
void sor_LDU(double[][N], double[], double[], double[], double[][N], double[][N], double, double, double);
void output_solutions(double[][N], double[], double[], const string&);

int main() {
    double A[N][N], b[N], x[N], L[N][N], U[N][N], D[N];

    generate_matrix(A, b);
    cout << "Matrix A and b:" << endl;
    print_matrix(A, b);

    LDU_decomposition(A, L, D, U);

    cout << "--------------------------------------------------" << endl;
    print_LDU(L, D, U);

    cout << "--------------------------------------------------" << endl;
    cout << "Simple Iterations (Jacobi Method) using LDU:" << endl;
    simple_iteration_LDU(A, b, x, D, L, U, EPS, EPS1);
    output_solutions(A, b, x, "Simple Iterations");

    cout << "--------------------------------------------------" << endl;
    cout << "Gauss-Seidel Method using LDU:" << endl;
    gauss_seidel_LDU(A, b, x, D, L, U, EPS, EPS1);
    output_solutions(A, b, x, "Gauss-Seidel");

    cout << "--------------------------------------------------" << endl;
    cout << "Successive Over-Relaxation (SOR) Method using LDU:" << endl;
    sor_LDU(A, b, x, D, L, U, EPS, EPS1, OMEGA);
    output_solutions(A, b, x, "SOR");

    cout << "--------------------------------------------------" << endl;

    return EXIT_SUCCESS;
}

void print_LDU(double L[][N], double D[], double U[][N]) {
    cout << "Matrix L (Lower triangular):" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(8) << fixed << setprecision(3) << L[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Matrix D (Diagonal):" << endl;
    for (int i = 0; i < N; ++i) {
        cout << setw(8) << fixed << setprecision(3) << D[i] << " ";
        cout << endl;
    }

    cout << "Matrix U (Upper triangular):" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(8) << fixed << setprecision(3) << U[i][j] << " ";
        }
        cout << endl;
    }
}


void generate_matrix(double A[][N], double b[]) {
    double tempA[][N] = {
        {10, -1, 2, 0, 0, 0, 0, 0, 0, 0},
        {-1, 11, -1, 3, 0, 0, 0, 0, 0, 0},
        {2, -1, 10, -1, 0, 0, 0, 0, 0, 0},
        {0, 3, -1, 8, -2, 0, 0, 0, 0, 0},
        {0, 0, 0, -2, 7, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 12, -1, 0, 0, 0},
        {0, 0, 0, 0, 0, -1, 9, -1, 0, 0},
        {0, 0, 0, 0, 0, 0, -1, 7, -2, 0},
        {0, 0, 0, 0, 0, 0, 0, -2, 8, 3},
        {0, 0, 0, 0, 0, 0, 0, 0, 3, 10}
    };

    double tempB[] = { 6, 25, -11, 15, 9, 10, -6, 5, 2, 3 };

    for (int i = 0; i < N; ++i) {
        b[i] = tempB[i];
        for (int j = 0; j < N; ++j) {
            A[i][j] = tempA[i][j];
        }
    }
}

void print_matrix(double A[][N], double b[]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(8) << fixed << setprecision(3) << A[i][j] << " ";
        }
        cout << "| " << setw(8) << b[i] << endl;
    }
}

void LDU_decomposition(double A[][N], double L[][N], double D[], double U[][N]) {
    double LU[N * N];
    int ipiv[N], info;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            LU[i * N + j] = A[i][j];

    dgetrf(&N, &N, LU, &N, ipiv, &info);

    for (int i = 0; i < N; ++i) {
        D[i] = LU[i * N + i];  // диагональ

        for (int j = 0; j < N; ++j) {
            if (i > j) {
                L[i][j] = LU[i * N + j];
                U[i][j] = 0.0;
            }
            else if (i == j) {
                L[i][j] = 0.0;
                U[i][j] = 0.0;
            }
            else {
                L[i][j] = 0.0;
                U[i][j] = LU[i * N + j] / D[i];  
            }
        }
    }
}


void simple_iteration_LDU(double A[][N], double b[], double x[], double D[], double L[][N], double U[][N], double epsilon, double epsilon1) {
    double x_new[N];
    double error = epsilon + 1.0;
    int iter = 0;

    for (int i = 0; i < N; ++i) x[i] = 0.0;

    while (error > epsilon && iter < 1000) {
        ++iter;
        error = 0.0;

        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                if (j != i) sum += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum) / D[i];
            error = max(error, fabs(x_new[i] - x[i]));
        }

        for (int i = 0; i < N; ++i) x[i] = x_new[i];

        double residual = 0.0;
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j)
                sum += A[i][j] * x[j];

            residual = max(residual, fabs(sum - b[i]));
        }

        if (residual < epsilon1) break;

        cout << "Iteration " << setw(5) << iter << ", error: " << error << ", residual: " << residual << endl;
    }
}

void gauss_seidel_LDU(double A[][N], double b[], double x[], double D[], double L[][N], double U[][N], double epsilon, double epsilon1) {
    double error = epsilon + 1.0;
    int iter = 0;

    for (int i = 0; i < N; ++i) x[i] = 0.0;

    while (error > epsilon && iter < 1000) {
        ++iter;
        error = 0.0;

        for (int i = 0; i < N; ++i) {
            double sum = b[i];
            for (int j = 0; j < N; ++j) {
                if (j != i)
                    sum -= A[i][j] * x[j];
            }
            double new_xi = sum / D[i];
            error = max(error, fabs(new_xi - x[i]));
            x[i] = new_xi;
        }

        double residual = 0.0;
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j)
                sum += A[i][j] * x[j];

            residual = max(residual, fabs(sum - b[i]));
        }

        if (residual < epsilon1) break;

        cout << "Iteration " << setw(3) << iter << ", error: " << error << setprecision(10) <<
            ", residual: " << residual << endl;
    }
}

void sor_LDU(double A[][N], double b[], double x[], double D[], double L[][N], double U[][N], double epsilon, double epsilon1, double omega) {
    double error = epsilon + 1.0;
    int iter = 0;

    for (int i = 0; i < N; ++i) x[i] = 0.0;

    while (error > epsilon && iter < 1000) {
        ++iter;
        error = 0.0;

        for (int i = 0; i < N; ++i) {
            double x_old = x[i];
            double sum = b[i];
            for (int j = 0; j < N; ++j) {
                if (j != i)
                    sum -= A[i][j] * x[j];
            }
            double gs_update = sum / D[i];
            x[i] = (1 - omega) * x_old + omega * gs_update;

            error = max(error, fabs(x[i] - x_old));
        }

        double residual = 0.0;
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                sum += A[i][j] * x[j];
            }
            residual = max(residual, fabs(sum - b[i]));
        }

        if (residual < epsilon1) break;

        cout << "Iteration " << setw(3) << iter << ", error: " << error << setprecision(10) <<
            ", residual: " << residual << endl;
    }
}

void output_solutions(double A[][N], double b[], double x[], const string& method_name) {
    double Ax[N];
    double residual = 0.0;

    for (int i = 0; i < N; ++i) {
        Ax[i] = 0.0;
        for (int j = 0; j < N; ++j)
            Ax[i] += A[i][j] * x[j];
        residual = max(residual, fabs(Ax[i] - b[i]));
    }

    cout << "Solution (" << method_name << "):" << endl;
    for (int i = 0; i < N; ++i)
        cout << "x[" << i + 1 << "] = " << fixed << setprecision(6) << x[i] << endl;

    cout << "Residual: " << residual << endl;
}
