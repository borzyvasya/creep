#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

const int Q = 10;
const int Q_EXT = Q * 2;
constexpr int Q_WITH_ONE = Q + 1;

void randMatrix(float[][Q_WITH_ONE]);          // Генерация матрицы Qx(Q+1)
void printMatrix(float[][Q_WITH_ONE]);         // Вывод матрицы Qx(Q+1)
void forwardSubstitution(float[][Q_WITH_ONE]); // Прямой ход для Qx(Q+1)
void backSubstitution(float[][Q_WITH_ONE], float[]); // Обратный ход
float calculateMatrixNorm(float[][Q_WITH_ONE]);    // Норма матрицы QxQ
float calculateInverseMatrixNorm(float[][Q_WITH_ONE]); // Норма обратной матрицы
void checkSolution(float[][Q_WITH_ONE], float[]);

int main() {
    float matrix[Q][Q_WITH_ONE];
    float solution[Q];
    randMatrix(matrix);

    cout << "Original matrix with b:" << endl;
    printMatrix(matrix);

    // Вычисляем решение, изменяя matrix напрямую
    forwardSubstitution(matrix);
    backSubstitution(matrix, solution);


    cout << "Checking solution:" << endl;
    checkSolution(matrix, solution);

    // Число обусловленности вычисляем отдельно с копированием
    float obuslov = calculateMatrixNorm(matrix) * calculateInverseMatrixNorm(matrix);
    cout << "Condition number: " << obuslov << endl;

    return EXIT_SUCCESS;
}

void randMatrix(float matrix[][Q_WITH_ONE]) {
    float CF = 1.0e+04;
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

void printExtendedMatrix(float matrix[][Q_EXT]) {
    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < Q_EXT; ++j) {
            cout << setw(8) << fixed << setprecision(2) << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void forwardSubstitution(float matrix[][Q_WITH_ONE]) {
    for (int i = 0; i < Q; i++) {
        float maxEl = fabs(matrix[i][i]);
        int maxRow = i;

        for (int k = i + 1; k < Q; k++) {
            if (fabs(matrix[k][i]) > maxEl) {
                maxEl = fabs(matrix[k][i]);
                maxRow = k;
            }
        }

        if (maxRow != i) {
            for (int j = 0; j < Q_WITH_ONE; j++) { // Q+1 столбцов
                float temp = matrix[maxRow][j];
                matrix[maxRow][j] = matrix[i][j];
                matrix[i][j] = temp;
            }
        }

        for (int k = i + 1; k < Q; k++) {
            float factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j < Q_WITH_ONE; j++) {
                matrix[k][j] -= factor * matrix[i][j];
                if (fabs(matrix[k][j]) < 1e-4)
                    matrix[k][j] = 0.0;
            }
        }
    }
}

void backSubstitution(float matrix[][Q_WITH_ONE], float solution[]) {
    for (int i = Q - 1; i >= 0; i--) {
        solution[i] = matrix[i][Q]; // b находится в столбце Q
        for (int j = i + 1; j < Q; j++) {
            solution[i] -= matrix[i][j] * solution[j];
        }
        if (fabs(matrix[i][i]) < 1e-9) {
            cout << "Matrix is singular in back substitution!" << endl;
            return;
        }
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
    float augmented[Q][Q * 2]; // Расширенная матрица (A | I)

    // Заполняем расширенную матрицу (A | I)
    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < Q; ++j) {
            augmented[i][j] = matrix[i][j];
        }
        for (int j = 0; j < Q; ++j) {
            augmented[i][j + Q] = (i == j) ? 1.0f : 0.0f;
        }
    }

    // Прямой ход метода Гаусса без перестановки строк
    for (int i = 0; i < Q; i++) {
        if (fabs(augmented[i][i]) < 1e-9) {
            cout << "Matrix is singular, cannot compute inverse!" << endl;
            return -1.0f;
        }
        for (int k = i + 1; k < Q; k++) {
            float factor = augmented[k][i] / augmented[i][i];
            for (int j = 0; j < 2 * Q; j++) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Обратный ход метода Гаусса
    for (int i = Q - 1; i >= 0; i--) {
        float pivot = augmented[i][i];
        if (fabs(pivot) < 1e-9) {
            cout << "Matrix is singular, cannot compute inverse!" << endl;
            return -1.0f;
        }
        for (int j = 0; j < 2 * Q; j++) {
            augmented[i][j] /= pivot;
        }
        for (int k = 0; k < i; k++) {
            float factor = augmented[k][i];
            for (int j = 0; j < 2 * Q; j++) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Вычисление нормы обратной матрицы
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
    for (int i = 0; i < Q; i++) {
        float sum = 0;
        for (int j = 0; j < Q; j++) sum += matrix[i][j] * solution[j];

        if (fabs(sum - matrix[i][Q]) > 1e-6) {
            cout << "Equation " << i + 1 << " is not satisfied: " << sum
                << " != " << matrix[i][Q] << endl;
            cout << "Solution is incorrect" << endl;
            return;
        }
    }
    cout << "Solution is correct" << endl;
}

