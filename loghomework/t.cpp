#include <iostream>

using namespace std;

int main() {
    // Перебираем возможные значения признаков сосуда (греческий, финикийский)
    for (int greek = 0; greek <= 1; greek++) {
        for (int phoenician = 0; phoenician <= 1; phoenician++) {
            // Пропускаем, если сосуд одновременно и греческий, и финикийский, т.к. это невозможно
            if (greek && phoenician) continue;
            // Перебираем возможные века (3-5)
            for (int century = 3; century <= 5; century++) {
                // Проверка, сколько утверждений верно (только одно утверждение должно быть верным)
                int correctStatements = 0;
                correctStatements += (greek == 1 && century == 5);
                correctStatements += (phoenician == 1 && century == 3);
                correctStatements += (greek == 0 && century == 4);

                if (correctStatements == 1) {
                    // Выводим правильный результат
                    cout << "The vessel is " << (greek ? "Greek" : "Phoenician")
                        << " and was made in the " << century << "th century." << endl;
                    return EXIT_SUCCESS;
                }
            }
        }
    }
}
