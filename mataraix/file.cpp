#include <iostream>
#include <iomanip>

using namespace std;

const int N = 3;
const double EPS = 0.0003;
const double OMEGA = 1.15;

void generate_matrix(double[][N], double[]);
void print_matrix(double[][N], double b[]);
void simple_iteration(double[][N], double[], double[], double);
void gauss_seidel(double[][N], double[], double[], double);
void sor(double[][N], double[], double[], double, double);
void output_solutions(double[][N], double[], double[], const string&);

int main() {
    double A[N][N], b[N], x[N];

    generate_matrix(A, b);
    cout << "Matrix A and b:" << endl;
    print_matrix(A, b);

    cout << "--------------------------------------------------" << endl;
    cout << "Simple Iterations (Jacobi Method):" << endl;
    simple_iteration(A, b, x, EPS);
    output_solutions(A, b, x, "Simple Iterations");

    cout << "--------------------------------------------------" << endl;
    cout << "Gauss-Seidel Method:" << endl;
    gauss_seidel(A, b, x, EPS);
    output_solutions(A, b, x, "Gauss-Seidel");

    cout << "--------------------------------------------------" << endl;
    cout << "Successive Over-Relaxation (SOR) Method:" << endl;
    sor(A, b, x, EPS, OMEGA);
    output_solutions(A, b, x, "SOR");

    cout << "--------------------------------------------------" << endl;

    return EXIT_SUCCESS;
}

void generate_matrix(double A[][N], double b[]) {
    double tempA[N][N] = {
        {3.0,  1.0, -1.0},
        {2.0,  6.0,  3.0},
        {-1.0, 1.0,  4.0}
    };
    double tempB[N] = { 7.0, -2.0, 4.0 };

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
            cout << setw(10) << fixed << setprecision(3) << A[i][j] << " ";
        }
        cout << "| " << setw(10) << b[i] << endl;
    }
}

void simple_iteration(double A[][N], double b[], double x[], double epsilon) {
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

            if (fabs(A[i][i]) < 1e-10) {
                cout << "Error: Near-zero diagonal element at row " << i << endl;
                return;
            }

            x_new[i] = (b[i] - sum) / A[i][i];
            error = max(error, fabs(x_new[i] - x[i]));
        }

        for (int i = 0; i < N; ++i)
            x[i] = x_new[i];

        cout << "Iteration " << setw(3) << iter << ", error: " << error << endl;
    }

}

void gauss_seidel(double A[][N], double b[], double x[], double epsilon) {
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
            if (fabs(A[i][i]) < 1.0e-10) {
                cout << "Error: Zero or near-zero diagonal element at " << i << endl;
                return;
            }

            double new_xi = sum / A[i][i];
            error = max(error, fabs(new_xi - x[i]));
            x[i] = new_xi;
        }

        cout << "Iteration " << setw(3) << iter << ", error: " << error << endl;
    }

    if (iter == 1000)
        cout << "Warning: Reached max iterations in Gauss-Seidel method" << endl;
}

void sor(double A[][N], double b[], double x[], double epsilon, double omega) {
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

            if (fabs(A[i][i]) < 1.0e-10) {
                cout << "Error: Zero or near-zero diagonal element at " << i << endl;
                return;
            }

            double gs_update = sum / A[i][i];
            x[i] = (1 - omega) * x_old + omega * gs_update;

            error = max(error, fabs(x[i] - x_old));
        }

        cout << "Iteration " << setw(3) << iter << ", error: " << error << endl;
    }

    if (iter == 1000)
        cout << "Warning: Reached max iterations in SOR method" << endl;
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