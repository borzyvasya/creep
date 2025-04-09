#include <iostream>

using namespace std;

float f(float x) {
    return pow(2, x) - 5 * x - 3;
}

int main() {
    float a = 0;        
    float b = 5;        
    float epsilon = 1e-6; 
    float step = 0.0001;  

    for (float x = a; x <= b; x += step) {
        if (fabs(f(x)) < epsilon) {
            cout << "Корень уравнения: " << x << endl;
            break;
        }
    }

    return EXIT_SUCCESS;
}
