#include <iostream>
  
using namespace std;

const int SIZE = 4;

float normL1(float[]); // L1
float normL2(float[]); // L2 
float l_inf_norm(float[]); // L inf for max

int main() {
    float x[SIZE] = { -12.7, 4.4, 2.3, 5.5 };

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
        sum += sqrt(arr[i] * arr[i]); 

    return sum; 
}

float l_inf_norm(float arr[]) {
    float max_value = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        float abs_value = abs(arr[i]);
        if (abs_value > max_value) max_value = abs_value;
    }

    return max_value;
}