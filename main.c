#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FIREFLY_COUNT 300
#define N 200
#define MAX_ITER 4000

#define ALPHA 0.2
#define BETA0 1.0
#define GAMMA 0.01

#define MINVAL -2.0
#define MAXVAL 2.0

double drand(uint *seed) { return ((double)rand_r(seed) / RAND_MAX); }

double norm2(double *x, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  return sqrt(sum);
}

double firefly_brightness(double f_x) { return 1 / (1 + f_x); }

void move(double *a, double *b, double alpha, unsigned int *seed) {

  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    double d = a[i] - b[i];
    sum += d * d;
  }
  double distance = sqrt(sum);
  double beta = BETA0 * exp(-GAMMA * distance * distance);
  for (int n = 0; n < N; n++) {
    a[n] += beta * (b[n] - a[n]) + alpha * (drand(seed) - 0.5);

    if (a[n] < MINVAL) {
      a[n] = MINVAL;
    }
    if (a[n] > MAXVAL) {
      a[n] = MAXVAL;
    }
  }
}

double quadratic(double *x, int n) {
  double sum = 0.0;
  for (int i = 2; i < n; i++) {
    sum += 100.0 * (x[i] * x[i] + x[i - 1] * x[i - 1]) + x[i - 2] * x[i - 2];
  }
  return sum;
}

double arrowhead(double *x, int n) {
  double sum = 0.0;
  double xn2 = x[n - 1] * x[n - 1];

  for (int i = 0; i < n - 1; i++) {
    double xi2 = x[i] * x[i];
    double sum2 = (xi2 + xn2);
    sum += sum2 * sum2 - 4.0 * x[i] + 3.0;
  }

  return sum;
}

void firefly(double (*fun)(double *, int)) {
  double fireflies[FIREFLY_COUNT][N];
  double brightness[FIREFLY_COUNT];

  for (int i = 0; i < FIREFLY_COUNT; i++) {
    unsigned int seed = 1234 + i;
    for (int k = 0; k < N; k++) {
      fireflies[i][k] = MINVAL + (MAXVAL - MINVAL) * drand(&seed);
    }
    brightness[i] = firefly_brightness(fun(fireflies[i], N));
  }

  for (int iter = 1; iter <= MAX_ITER; iter++) {
    double alpha = ALPHA * (1.0 - (double)iter / MAX_ITER);
    unsigned int seed = 4321;
    for (int i = 0; i < FIREFLY_COUNT; i++) {
      double new_pos[N];
      for (int n = 0; n < N; n++) {
        new_pos[n] = fireflies[i][n];
      }
      double current = brightness[i];
      for (int j = 0; j < FIREFLY_COUNT; j++) {
        if (brightness[j] > current) {
          move(new_pos, fireflies[j], alpha, &seed);
        }
      }
      for (int n = 0; n < N; n++) {
        fireflies[i][n] = new_pos[n];
      }
      brightness[i] = firefly_brightness(fun(fireflies[i], N));
    }
    if (iter % 200 == 0) {
      int best = 0;
      for (int i = 0; i < FIREFLY_COUNT; i++) {
        if (brightness[i] > brightness[best]) {
          best = i;
        }
      }
      double f_x = fun(fireflies[best], N);
      double x_norm = norm2(fireflies[best], N);
      printf("iter %4d | max brightness = %.6e | min f(x) = %.10e | ||x||_2 = "
             "%.10e\n",
             iter, brightness[best], f_x, x_norm);
    }
  }

  int best = 0;
  for (int i = 0; i < FIREFLY_COUNT; i++) {
    if (brightness[i] > brightness[best]) {
      best = i;
    }
  }

  double f_x = fun(fireflies[best], N);
  double x_norm = norm2(fireflies[best], N);
  printf("\n ~~~~~~~~WYNIK~~~~~~~~ \n");
  printf("Max(Brightness) = %.10e\n", brightness[best]);
  printf("Min(f(x)) = %.10e\n", f_x);
  printf("||x||_2 = %.10e\n", x_norm);
}

int main() {
  printf("\n ~~~~~~~~QUADRATIC~~~~~~~~ \n");
  firefly(quadratic);
  printf("\n ~~~~~~~~ARROWHEAD~~~~~~~~ \n");
  firefly(arrowhead);
  return 0;
}
