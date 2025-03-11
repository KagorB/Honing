#include <iostream>
#include <cmath>
#define PI 3.141592653589793
using namespace std;

int main() {
  setlocale(LC_ALL, "Russian");
  // Данные
  double l = 0.175;          // Длина хвостовика (м)
  double R = 0.035;     // Радиус хвостовика (м)
  double r = R / 2;        // Радиус головки (м)
  double G = 80e9;         // Модуль сдвига (Па)
  double P = 9090;          // Давление (Па)
  double v = 100;          // Скорость (м/с)
  double p = 7850;         // Плотность (кг/м^3)
  double f = 0.05;         // Коэффициент трения
  double f_min = 0.05;
  double f_max = 0.15;
  double N = 100;
  double N_min = 100;
  double N_max = 500;
  double A = 2 * PI * R * 0.1;
  double F = 200;
  double F_min = 200;

  int count = 1;

  for (N = N_min, F = F_min; N <= N_max; N += 25, F += 10) {
    for (f = f_min; f <= f_max; f += 0.01)
    {
      P = F / A;
      double omega = 2 * PI * N / 60; // радианы в секунду
      double J0 = (1.0 / 2) * p * PI * pow(r, 4) * l;        // J0
      double J1 = J0 + (PI * pow(r, 4) * p * l / 6.0);       // J1
      double X = (PI * pow(r, 4) * G) / (2.0 * l);           // X
      double Tetta = sqrt(X * J1) / (R * P * f) * omega;
      // Вывод результатов
      if (Tetta >= 1 && Tetta <= 2.1) {
        cout << count << ": ";
        cout << "F: " << F; cout << " P: " << P;
        cout << " N: " << N;
        cout << " T: " << Tetta << endl;
        count++;
      }
    }
  }
  return 0;
}


/*
int main() {
  setlocale(LC_ALL, "Russian");
  // Данные
  double l = 0.175;          // Длина хвостовика (м)
          // Радиус хвостовика (м)
  double R = 0.035;
  double r = R/2;        // Радиус головки (м)
  double G = 80e9;         // Модуль сдвига (Па)
  double P = 9090;          // Давление (Па)
  double v = 100;          // Скорость (м/с)
  double p = 7850;         // Плотность (кг/м^3)
  double f1 = 0.03;         // Коэффициент трения
  double f1_min = 0.05;
  double f1_max = 0.15;

  // Вычисления
  double N = 100; //обороты в минуту
  double N_min = 100;
  double N_max = 500;
  double A = 2 * PI * R * 0.1;
  double F = 100;
  double F_min = 200;
  double F_max = 500;
  P = F / A;
  cout << A << endl;
  P = 50000;
  double omega = 2 * PI * N / 60; // радианы в секунду
  double J0 = (1.0 / 2) * p * PI * pow(r, 4) * l;        // J0
  double J1 = J0 + (PI * pow(r, 4) * p * l / 6.0);       // J1
  double X = (PI * pow(r, 4) * G) / (2.0 * l);           // X
  double Tetta = sqrt(X * J1) / (R * P * f1)* omega;            // Tetta

  // Вывод результатов
  cout << "Результаты вычислений:" << endl;
  cout << "Давление:" << P << endl;
  cout << "Угловая скорость (omega): " << omega << " рад/с" << endl;
  cout << "J0: " << J0 << " кг·м^2" << endl;
  cout << "J1: " << J1 << " кг·м^2" << endl;
  cout << " (X): " << X << " Н·м/рад" << endl;
  cout << "Tetta: " << Tetta << endl;

  return 0;
}
*/