#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FIREFLY_COUNT 200
#define N 100
#define MAX_ITER 4000

#define ALPHA 0.2
#define BETA0 1.0
#define GAMMA 0.01

#define MINVAL -2.0
#define MAXVAL 2.0

double drand() { return ((double)rand() / RAND_MAX); } // 0.0-1.0

double firefly_brightness(double f_x) {
  return 1 / (1 + f_x);
} // szukamy minimum, wiec jasnosc wieksza, im mniejsza wartosc funkcji

void move(double *a, double *b) {
  // odleglosc
  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    double d = a[i] - b[i];
    sum += d * d;
  }
  double distance = sqrt(sum);

  double beta = BETA0 * exp(-GAMMA * distance * distance);
}

double quadratic(double *x, int n) {
  double sum = 0.0;
  for (int i = 2; i < n; i++) {
    sum += 100.0 * (x[i] * x[i] + x[i - 1] * x[i - 1]) + x[i - 2] * x[i - 2];
  }
  return sum;
}

int main() {
  srand(time(NULL));

  double fireflies[FIREFLY_COUNT][N];
  double brightness[FIREFLY_COUNT];

  // losowe wartosci na poczatku i wyliczenie jasnosci
  for (int i = 0; i < FIREFLY_COUNT; i++) {
    for (int k = 0; k < N; k++) {
      fireflies[i][k] = MINVAL + (MAXVAL - MINVAL) * drand();
    }
    brightness[i] = firefly_brightness(quadratic(fireflies[i], N));
  }

  clock_t t_start = clock();

  for (int iter = 0; iter < MAX_ITER; iter++) {
    for (int i = 0; i < FIREFLY_COUNT; i++) {
      for (int j = 0; j < FIREFLY_COUNT; j++) {
        if (brightness[j] < brightness[i]) {
          double pos[N];
          double r = distance(fireflies[i], fireflies[j], N);
          double new_beta = BETA0 * exp(-GAMMA * r * r);

          for (int n = 0; n < N; n++) {
            fireflies[i][n] = new_beta * (fireflies[j][n] - fireflies[i][n]) +
                              ALPHA * (drand() - 0.5);
          }
        }
      }
      brightness[i] = quadratic(fireflies[i], N);
    }
    // co 200 iteracji wypisz wynik
    if (iter % 200 == 0) {
      int best = 0;
      for (int i = 0; i < FIREFLY_COUNT; i++) {
        if (brightness[i] < brightness[best])
          best = i;
        printf("iter %4d | best f(x) = %.6e\n", iter, brightness[best]);
      }
    }
  }

  clock_t t_end = clock();

  int best = 0;
  for (int i = 0; i < FIREFLY_COUNT; i++)
    if (brightness[i] < brightness[best])
      best = i;
  printf("\n ~~~~~~~~WYNIK~~~~~~~~ \n");
  printf("Best f(x) = %.10e\n", brightness[best]);
  // printf("x = ");
  // for (int n = 0; n < N; n++) {
  //   printf("%.6f ", fireflies[best][n]);
  // }
  // printf("\n");
  double time_sec = (double)(t_end - t_start) / CLOCKS_PER_SEC;
  printf("Czas sekwencyjny: %.2f s\n", time_sec);
  return 0;

  return 0;
}
