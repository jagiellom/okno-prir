#include <math.h>
#include <mpi.h>
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

double drand_r(unsigned int *seed) { return (double)rand_r(seed) / RAND_MAX; }

double norm2(double *x, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    sum += x[i] * x[i];
  return sqrt(sum);
}

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

void firefly(double (*fun)(double *, int), int rank, int size) {
  static double fireflies[FIREFLY_COUNT][N];
  static double brightness[FIREFLY_COUNT];

  int _n = FIREFLY_COUNT / size;
  int start = rank * _n;
  int end = start + _n;

  unsigned int seed = 1234 + rank;

  for (int i = start; i < end; i++) {
    for (int k = 0; k < N; k++)
      fireflies[i][k] = MINVAL + (MAXVAL - MINVAL) * drand_r(&seed);

    brightness[i] = firefly_brightness(fun(fireflies[i], N));
  }

  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fireflies, _n * N,
                MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, brightness, _n, MPI_DOUBLE,
                MPI_COMM_WORLD);

  for (int iter = 1; iter <= MAX_ITER; iter++) {

    double alpha = ALPHA * (1.0 - (double)iter / MAX_ITER);

    for (int i = start; i < end; i++) {

      double new_pos[N];
      for (int n = 0; n < N; n++)
        new_pos[n] = fireflies[i][n];

      for (int j = 0; j < FIREFLY_COUNT; j++) {
        if (brightness[j] > brightness[i]) {
          move(new_pos, fireflies[j], alpha, &seed);
        }
      }

      for (int n = 0; n < N; n++)
        fireflies[i][n] = new_pos[n];

      brightness[i] = firefly_brightness(fun(fireflies[i], N));
    }

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fireflies, _n * N,
                  MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, brightness, _n,
                  MPI_DOUBLE, MPI_COMM_WORLD);

    if (rank == 0 && iter % 200 == 0) {
      int best = 0;
      for (int i = 1; i < FIREFLY_COUNT; i++)
        if (brightness[i] > brightness[best])
          best = i;

      double f_x = fun(fireflies[best], N);
      double x_norm = norm2(fireflies[best], N);
      printf("iter %4d | max brightness = %.6e | min f(x) = %.10e | ||x||_2 = "
             "%.10e\n",
             iter, brightness[best], f_x, x_norm);
    }
  }

  if (rank == 0) {
    int best = 0;
    for (int i = 1; i < FIREFLY_COUNT; i++)
      if (brightness[i] > brightness[best])
        best = i;

    printf("\n===== WYNIK =====\n");
    printf("Min f(x) = %.10e\n", fun(fireflies[best], N));
    printf("Max(Brightness) = %.10e\n", brightness[best]);
    printf("||x||_2  = %.10e\n", norm2(fireflies[best], N));
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (FIREFLY_COUNT % size != 0) {
    if (rank == 0)
      fprintf(stderr, "FIREFLY_COUNT must be divisible by MPI size\n");
    MPI_Finalize();
    return 1;
  }

  if (rank == 0)
    printf("\n==== QUADRATIC ====\n");
  firefly(quadratic, rank, size);

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    printf("\n==== ARROWHEAD ====\n");
  firefly(arrowhead, rank, size);

  MPI_Finalize();
  return 0;
}
