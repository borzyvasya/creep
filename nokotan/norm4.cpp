//‖x‖₂ = (⟨x, x⟩) ^ (1 / 2)


#include <iostream>
#include <iomanip>

using namespace std;

const int SIZE = 4;

float product(float[]);
float norm(float[]);

int main() {
    float x[SIZE] = { 3.0, -4.5, 2.1, -6.7 };  
    cout << setprecision(4) << "Norm of the vector: " << norm(x) << endl;

    return EXIT_SUCCESS;
}

float product(float arr[]) {
    float sum = 0.0;
    for (int i = 0; i < SIZE; ++i) 
        sum += arr[i] * arr[i];  // x_i * x_i

    return sum;
}

float norm(float arr[]) {
    return sqrt(product(arr));
}

