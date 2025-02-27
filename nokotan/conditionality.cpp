#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>

using namespace std;

const int Q = 10;
constexpr int Q_WITH_ONE = Q + 1;
const int Q_EXT = Q * 2; // Для создания расширенный матрицы в которой исходная и единичная

void randMatrix(float[][Q_WITH_ONE]);
void printMatrix(float[][Q_WITH_ONE]);
void forwardSubstitution(float[][Q_WITH_ONE], float[]);
float calculateMatrixNorm(float[][Q_WITH_ONE]);
float calculateInverseMatrixNorm(float[][Q_WITH_ONE]);
void checkSolution(float[][Q_WITH_ONE], float[]);

int main() {
    float matrix[Q][Q_WITH_ONE];
    float solution[Q];
    randMatrix(matrix);

    cout << "Original matrix:" << endl;
    printMatrix(matrix);

    forwardSubstitution(matrix, solution);

    cout << "Upper triangular matrix:" << endl;
    printMatrix(matrix);

    cout << "Checking solution: "; checkSolution(matrix, solution);

    float conditionNumber = calculateMatrixNorm(matrix) * calculateInverseMatrixNorm(matrix);
    cout << "Condition number: " << conditionNumber << endl;

    return EXIT_SUCCESS;
}

void randMatrix(float matrix[][Q_WITH_ONE]) {
    float CF = 1.0e+03;
    srand(time(NULL));
    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < Q_WITH_ONE; ++j) {
            matrix[i][j] = rand() / CF;
        }
    }
}

void printMatrix(float matrix[][Q_WITH_ONE]) {
    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < Q; ++j) {
            cout << setw(8) << fixed << setprecision(2) << matrix[i][j] << " ";
        }
        cout << " | " << setw(8) << fixed << setprecision(2) << matrix[i][Q];
        cout << endl;
    }
    cout << endl;
}

void forwardSubstitution(float matrix[][Q_WITH_ONE], float solution[]) {
    for (int i = 0; i < Q; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < Q; ++k) {
            if (fabs(matrix[k][i]) > fabs(matrix[maxRow][i]))
                maxRow = k;
        }

        if (maxRow != i) {
            for (int j = 0; j < Q_WITH_ONE; ++j) {
                float temp = matrix[maxRow][j];
                matrix[maxRow][j] = matrix[i][j];
                matrix[i][j] = temp;
            }
        }

        if (fabs(matrix[i][i]) < 1e-6) {
            cout << "Matrix is singular!" << endl;
            return;
        }

        // Приведение к верхнетреугольному виду
        for (int k = i + 1; k < Q; ++k) {
            float factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j < Q_WITH_ONE; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
                if (fabs(matrix[k][j]) < 1e-4) matrix[k][j] = 0.0;
            }
        }
    }

    // Обратный ход (Back Substitution) для solutions[]
    for (int i = Q - 1; i >= 0; --i) {
        solution[i] = matrix[i][Q];
        for (int j = i + 1; j < Q; ++j) 
            solution[i] -= matrix[i][j] * solution[j];
 
        solution[i] /= matrix[i][i];
    }
}

float calculateMatrixNorm(float matrix[][Q_WITH_ONE]) {
    float norm = 0.0;
    for (int j = 0; j < Q; ++j) {
        float colSum = 0.0;
        for (int i = 0; i < Q; ++i) {
            colSum += fabs(matrix[i][j]);
        }
        norm = max(norm, colSum);
    }
    return norm;
}

float calculateInverseMatrixNorm(float matrix[][Q_WITH_ONE]) {
    float augmented[Q][Q_EXT];

    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < Q; ++j) {
            augmented[i][j] = matrix[i][j];
        }
        for (int j = 0; j < Q; ++j) {
            augmented[i][j + Q] = (i == j) ? 1.0f : 0.0f; 
        }
    }

    for (int i = 0; i < Q; i++) {
        if (fabs(augmented[i][i]) < 1e-9) {
            cout << "Matrix is singular" << endl;
            return -1.0f;
        }
        for (int k = i + 1; k < Q; k++) {
            float factor = augmented[k][i] / augmented[i][i];
            for (int j = 0; j < 2 * Q; ++j) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Обратный ход для augmented
    for (int i = Q - 1; i >= 0; i--) {
        for (int j = Q; j < 2 * Q; j++) {
            augmented[i][j] /= augmented[i][i];
        }
        for (int k = i - 1; k >= 0; k--) {
            float factor = augmented[k][i];
            for (int j = Q; j < 2 * Q; j++) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    float norm = 0.0f;
    for (int j = Q; j < 2 * Q; j++) {
        float colSum = 0.0f;
        for (int i = 0; i < Q; i++) {
            colSum += fabs(augmented[i][j]);
        }
        norm = max(norm, colSum);
    }
    return norm;
}

void checkSolution(float matrix[][Q_WITH_ONE], float solution[]) {
    for (int i = 0; i < Q; ++i) {
        float sum = 0;
        for (int j = 0; j < Q; ++j) sum += matrix[i][j] * solution[j];

        if (fabs(sum - matrix[i][Q]) > 1e-4) {
            cout << "Solution is incorrect" << endl;
            return;
        }
    }
    cout << "Solution is correct" << endl;
}
