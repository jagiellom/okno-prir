#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define FIREFLY_COUNT 200
#define N 100
#define MAX_ITER 4000

#define ALPHA 0.2
#define BETA0 1.0
#define GAMMA 0.01

#define MINVAL -2.0
#define MAXVAL 2.0

/* ===================== RNG ===================== */

double drand_r(unsigned int *seed) { return (double)rand_r(seed) / RAND_MAX; }

/* ===================== NORM ===================== */

double norm2(double *x, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    sum += x[i] * x[i];
  return sqrt(sum);
}

/* ===================== FA CORE ===================== */

double firefly_brightness(double f_x) { return 1.0 / (1.0 + f_x); }

void move(double *a, double *b, double alpha, unsigned int *seed) {
  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    double d = a[i] - b[i];
    sum += d * d;
  }

  double distance = sqrt(sum);
  double beta = BETA0 * exp(-GAMMA * distance * distance);

  for (int n = 0; n < N; n++) {
    a[n] += beta * (b[n] - a[n]) + alpha * (drand_r(seed) - 0.5);

    if (a[n] < MINVAL)
      a[n] = MINVAL;
    if (a[n] > MAXVAL)
      a[n] = MAXVAL;
  }
}

/* ===================== TEST FUNCTIONS ===================== */

double quadratic(double *x, int n) {
  double sum = 0.0;
  for (int i = 2; i < n; i++)
    sum += 100.0 * (x[i] * x[i] + x[i - 1] * x[i - 1]) + x[i - 2] * x[i - 2];
  return sum;
}

double arrowhead(double *x, int n) {
  double sum = 0.0;
  double xn2 = x[n - 1] * x[n - 1];

  for (int i = 0; i < n - 1; i++) {
    double xi2 = x[i] * x[i];
    double t = xi2 + xn2;
    sum += t * t - 4.0 * x[i] + 3.0;
  }
  return sum;
}

/* ===================== FIREFLY ALGORITHM ===================== */

void firefly(double (*fun)(double *, int)) {

  static double fireflies[FIREFLY_COUNT][N];
  static double brightness[FIREFLY_COUNT];

  /* initialization */
  for (int i = 0; i < FIREFLY_COUNT; i++) {
    unsigned int seed = 1234 + i;
    for (int k = 0; k < N; k++)
      fireflies[i][k] = MINVAL + (MAXVAL - MINVAL) * drand_r(&seed);
    brightness[i] = firefly_brightness(fun(fireflies[i], N));
  }

  double t_start = omp_get_wtime();

  for (int iter = 0; iter < MAX_ITER; iter++) {
    double alpha = ALPHA * (1.0 - (double)iter / MAX_ITER);

#pragma omp parallel
    {
      unsigned int seed = 4321 + omp_get_thread_num();

#pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < FIREFLY_COUNT; i++) {

        double new_pos[N];
        for (int n = 0; n < N; n++)
          new_pos[n] = fireflies[i][n];

        double current = brightness[i];

        for (int j = 0; j < FIREFLY_COUNT; j++) {
          if (brightness[j] > current) {
            move(new_pos, fireflies[j], alpha, &seed);
          }
        }

        for (int n = 0; n < N; n++)
          fireflies[i][n] = new_pos[n];

        brightness[i] = firefly_brightness(fun(fireflies[i], N));
      }
    }

    if (iter % 200 == 0) {
      int best = 0;
      for (int i = 1; i < FIREFLY_COUNT; i++)
        if (brightness[i] > brightness[best])
          best = i;

      printf("iter %4d | min f(x) = %.6e | ||x||_2 = %.6e\n", iter,
             fun(fireflies[best], N), norm2(fireflies[best], N));
    }
  }

  double t_end = omp_get_wtime();

  int best = 0;
  for (int i = 1; i < FIREFLY_COUNT; i++)
    if (brightness[i] > brightness[best])
      best = i;

  printf("\n===== WYNIK =====\n");
  printf("Min f(x) = %.10e\n", fun(fireflies[best], N));
  printf("||x||_2  = %.10e\n", norm2(fireflies[best], N));
  printf("Czas OpenMP = %.3f s\n\n", t_end - t_start);
}

/* ===================== MAIN ===================== */

int main() {

  printf("\n==== QUADRATIC ====\n");
  firefly(quadratic);

  printf("\n==== ARROWHEAD ====\n");
  firefly(arrowhead);

  return 0;
}
