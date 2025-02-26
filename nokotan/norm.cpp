#include <iostream>
  
using namespace std;

const int SIZE = 4;

float normL1(float[]); // L1
float normL2(float[]); // L2 
float l_inf_norm(float[]); // L inf for max

int main() {
    float x[SIZE] = { 12.7, 4.4, 2.3, 5.5 };

    cout << "Norm of the array (L1): " << normL1(x) << endl;
    cout << "Norm of the array (L2): " << normL2(x) << endl;
    cout << "Norm of the array (Linf): " << l_inf_norm(x) << endl;

    return EXIT_SUCCESS;
}

float normL1(float arr[]) {
    float sum = 0.0;
    for (int i = 0; i < SIZE; ++i)
        sum += abs(arr[i]); // Модуль каждого элемента

    return sum;
}

float normL2(float arr[]) {
    float sum = 0.0;
    for (int i = 0; i < SIZE; ++i) 
        sum += arr[i] * arr[i]; // Квадраты элементов

    return sqrt(sum); // Квадратный корень из суммы квадратов
}

float l_inf_norm(float arr[]) {
    float max_val = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        float abs_val = abs(arr[i]);
        if (abs_val > max_val) max_val = abs_val;
    }

    return max_val;
}