//  ‖x‖₂ = (Σ |xᵢ|²)^(1/2)
//  i = 1 to m


#include <iostream>
#include <cmath>  

using namespace std;

const int SIZE = 4;
float norm(float[]);

int main() {
    float x[SIZE] = { 3.0, -4.5, 2.1, -6.7 };
    cout << "Norm of the array: " << norm(x) << endl;

    return EXIT_SUCCESS;
}

float norm(float arr[]) {
    float sum = 0.0;
    for (int i = 0; i < SIZE; ++i) {
        sum += arr[i] * arr[i];
    }

    return sqrt(sum);
}