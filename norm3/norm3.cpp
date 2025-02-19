//  ‖x‖∞ = max | xᵢ |
//  1≤i≤m


#include <iostream>
#include <cmath>  

using namespace std;
const int SIZE = 4;

float l_inf_norm(float[], int);

int main() {
    float x[SIZE] = { 3.0, -4.5, 2.1, -6.7 }; 
    cout << "Norm of the array: " << l_inf_norm(x, SIZE) << endl;

    return EXIT_SUCCESS;
}

float l_inf_norm(float arr[], int SIZE) {
    float max_val = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        float abs_val = abs(arr[i]);
        if (abs_val > max_val) max_val = abs_val;
    }

    return max_val;
}