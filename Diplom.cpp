#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#define PI 3.141592653589793
using namespace std;
using namespace chrono;

//Многопоточка
const int NUM_THREADS = 16;
mutex file_mutex;
double iter = 0.01;

double sign(double x) {
  if (x > 0) return 1.0;
  if (x < 0) return -1.0;
  return 0.0;
}

void derivatives(double t, double y[], double dy[], double A) {
  dy[0] = y[1];
  dy[1] = -y[0] - sign(y[1] - A);
}

void rungeKutta(double t, double y[], double dt, double A) {
  int n = 2;
  double k1[2], k2[2], k3[2], k4[2], y_temp[2];

  derivatives(t, y, k1, A);
  for (int i = 0; i < n; i++) y_temp[i] = y[i] + dt / 2 * k1[i];
  derivatives(t + dt / 2, y_temp, k2, A);
  for (int i = 0; i < n; i++) y_temp[i] = y[i] + dt / 2 * k2[i];
  derivatives(t + dt / 2, y_temp, k3, A);
  for (int i = 0; i < n; i++) y_temp[i] = y[i] + dt * k3[i];
  derivatives(t + dt, y_temp, k4, A);

  for (int i = 0; i < n; i++) {
    y[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  }
}

void worker(double A, double dt, double epsilon_start, double epsilon_end, double epsilon_step, ofstream& file) {
  double t = 0.0;

  for (double epsilon_star = epsilon_start; epsilon_star <= epsilon_end; epsilon_star += epsilon_step) {
    iter += epsilon_step;
    cout << iter << endl;
    vector<double> T;
    double y[2] = { 1.0, 2.0 };

    while (T.size() < 1000) {
      if (fabs(y[1] - A) < dt && -1 < y[0] && y[0] < 1) {
        double y0 = y[0];
        double t0 = 0.01;
        double epsilon = 0.0;
        double new_y = 0.0;

        while (true) {
          new_y = A * t0 + y0;
          epsilon = (t0 < epsilon_star) ? t0 : epsilon_star;

          if (new_y > 1 + epsilon) {
            T.push_back(t0);
            y[0] = new_y;
            y[1] = A;
            break;
          }
          t0 += dt;
        }
      }
      else {
        rungeKutta(t, y, dt, A);
        t += dt;
      }
    }

    lock_guard<mutex> lock(file_mutex);
    file << epsilon_star << endl;
    for (double val : T) {
      file << val << " ";
    }
    file << endl;
  }
}

int main() {
  setlocale(LC_ALL, "Russian");

  auto start_time = high_resolution_clock::now(); 

  double A = 1.7;
  double dt = 0.001;
  //double dt = 0.00001;
  double epsilon_star_0 = 0.01;
  double epsilon_star_end = 4.0;
  //double epsilon_step = 0.005;
  double epsilon_step = 0.005;

  ofstream file("T_values.txt");
  if (!file) {
    cerr << "Не удалось открыть файл для записи T_values.txt!" << endl;
    return 1;
  }
  file << A << endl;

  vector<thread> threads;
  int num_steps = (epsilon_star_end - epsilon_star_0) / epsilon_step;
  int steps_per_thread = num_steps / NUM_THREADS;

  for (int i = 0; i < NUM_THREADS; i++) {
    double start = epsilon_star_0 + i * steps_per_thread * epsilon_step;
    double end = (i == NUM_THREADS - 1) ? epsilon_star_end : start + (steps_per_thread - 1) * epsilon_step;
    threads.emplace_back(worker, A, dt, start, end, epsilon_step, ref(file));
  }

  for (auto& th : threads) {
    th.join();
  }

  file.close();
  cout << "Значения массива T записаны в файл T_values.txt" << endl;

  auto end_time = high_resolution_clock::now();
  duration<double> elapsed = end_time - start_time;
  cout << "Время выполнения: " << elapsed.count() << " секунд" << endl;

  return 0;
}



//Однопоточка
//double sign(double x) {
//  if (x > 0) return 1.0;
//  if (x < 0) return -1.0;
//  return 0.0;
//}
//
//// Функция для вычисления производных
//void derivatives(double t, double y[], double dy[], double A) {
//  dy[0] = y[1];  // y' = v
//  dy[1] = -y[0] - sign(y[1] - A);  // v' = -y - sign(v - A)
//}
//
//// Метод Рунге-Кутты 4-го порядка
//void rungeKutta(double t, double y[], double dt, double A) {
//  int n = 2;  // количество уравнений
//  double k1[2], k2[2], k3[2], k4[2], y_temp[2];
//
//  derivatives(t, y, k1, A);
//  for (int i = 0; i < n; i++) y_temp[i] = y[i] + dt / 2 * k1[i];
//  derivatives(t + dt / 2, y_temp, k2, A);
//  for (int i = 0; i < n; i++) y_temp[i] = y[i] + dt / 2 * k2[i];
//  derivatives(t + dt / 2, y_temp, k3, A);
//  for (int i = 0; i < n; i++) y_temp[i] = y[i] + dt * k3[i];
//  derivatives(t + dt, y_temp, k4, A);
//
//  for (int i = 0; i < n; i++) {
//    y[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
//  }
//}
//
//int main() {
//  // Начальные условия
//  setlocale(LC_ALL, "Russian");
//  //double l = 0.175;          // Длина хвостовика (м)
//  //double r = 0.0175;        // Радиус хвостовика (м)
//  //double R = 0.035;        // Радиус головки (м)
//  //double G = 80e9;         // Модуль сдвига (Па)
//  //double P = 1e6;          // Давление (Па)
//  //double v = 100;          // Скорость (м/с)
//  //double p = 7850;         // Плотность (кг/м^3)
//  //double f1 = 0.2;         // Коэффициент трения
//
//  // Вычисления
//  //double omega = v / R;                                  // Угловая скорость
//  //double J0 = (1.0 / 2) * p * PI * pow(r, 4) * l;        // J0
//  //double J1 = J0 + (PI * pow(r, 4) * p * l / 6.0);       // J1
//  //double X = (PI * pow(r, 4) * G) / (2.0 * l);           // X
//  //double Tetta = sqrt(X * J1) / (R * P * f1) * omega;    // Tetta
//  //cout << Tetta << endl;
//  //
//  ////double A = 1.6;
//  //double A = round(Tetta * 10) / 10;
//  //cout << A << endl;
//  double A = 1.6;
//  double n = A + 0.2;
//  double y[2] = { -0.5, n };  // y и y'
//    // значение A
//  double dt = 0.001;  // шаг по времени
//  int steps = 1000;  // количество шагов
//  double epsilon_star_end = 4.0;  // e*
//  double epsilon_star_0 = 0.01;
//
//  // Открытие файла для записи
//  ofstream file("T_values.txt");
//  if (!file) {
//    cerr << "Не удалось открыть файл для записи T_values.txt!" << endl;
//    return 1;
//  }
//  file << A << endl;
//  for (double epsilon_star = epsilon_star_0; epsilon_star <= epsilon_star_end; epsilon_star += 0.001) {
//    cout << epsilon_star<<endl;
//    vector<double> T;  // массив для значений t
//    double t = 0.0;  // начальное время
//
//    // Основной цикл интегрирования
//    while (T.size() < 10000) {
//      if (round(y[1] * 10) / 10 == A && -1 < y[0] && y[0] < 1) {
//        double y0 = y[0];
//        double t0 = 0.01;
//        double h = 0.001;
//        double epsilon = 0.0;
//        double new_y = 0.0;
//
//        while (true) {
//          new_y = A * t0 + y0;
//          epsilon = (t0 < epsilon_star) ? t0 : epsilon_star;
//
//          if (new_y > 1 + epsilon) {
//            T.push_back(t0);
//            y[0] = new_y;
//            y[1] = A;
//            break;
//          }
//          t0 += h;
//        }
//      }
//      else {
//        rungeKutta(t, y, dt, A);
//        t += dt;
//      }
//    }
//
//    // Запись значений массива T в файл
//    file << epsilon_star << endl;
//    for (double val : T) {
//      file << val << " ";
//    }
//    file << endl;
//  }
//
//  // Закрытие файла
//  file.close();
//
//  cout << "Значения массива T записаны в файл T_values.txt" << endl;
//
//  return 0;
//}
