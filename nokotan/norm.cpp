//‖x‖₁ = Σ |xᵢ|
// Σ i = 1 to m

#include <iostream>
  
using namespace std;

const int SIZE = 4;

// L1
float normL1(float arr[]) {
    float sum = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        sum += abs(arr[i]); // Модуль каждого элемента
    }
    return sum;
}

// L2
float normL2(float arr[]) {
    float sum = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        sum += arr[i] * arr[i]; // Квадраты элементов
    }
    return sqrt(sum); // Квадратный корень из суммы квадратов
}

// L infinite for max
float l_inf_norm(float arr[]) {
    float max_val = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        float abs_val = abs(arr[i]);
        if (abs_val > max_val) max_val = abs_val;
    }

    return max_val;
}

int main() {
    float x[SIZE] = { -1.0, 6.5, 3.4, 8.1 };

    cout << "Norm of the array (L1): " << normL1(x) << endl;
    cout << "Norm of the array (L2): " << normL2(x) << endl;
    cout << "Norm of the array (Linf): " << l_inf_norm(x) << endl;


    return EXIT_SUCCESS;
}