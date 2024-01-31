#include "../numan_utils.h"
#include <assert.h>
#include <fenv.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double det_T(const unsigned int n, double T[n][n]) {
  double res = 0;
  for (unsigned int i = 0; i < n; ++i)
    res += T[i][i];
  printf("det(T): %lf\n", res);
  return res;
}

double *back_sub(const unsigned int n, double U[n][n], double b[n]) {
  if (fabs(det_T(n, U)) > P) {
    double *x = calloc(n, sizeof(double));
    double res;
    for (unsigned int i = n - 1; i != UINT_MAX; i--) {
      res = 0;
      for (unsigned int j = i + 1; j < n; j++) {
        res += U[i][j] * x[j];
      }
      x[i] = (b[i] - res) / U[i][i];
    }
    return x;
  } else {
    printf("det(U) == 0");
    return NULL;
  }
}

double det(const unsigned n, double A[n][n]) {
  return A[n - 1][n - 1]; // TODO: implement
}

double *solve(const unsigned int n, double A[n][n], double b[n]) {
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      printf("%lf ", A[i][j]);
    }
    printf("|%lf| ", b[i]);
    printf("\n");
  }
  printf("\n");

  if (fabs(det(n, A)) > P) {
    for (unsigned int i = 0; i < n - 1; ++i) {
      double a = A[i][i];
      for (unsigned int j = i + 1; j < n; ++j) {
        double d = A[j][i];
        double a_d = d / a;
        b[j] -= a_d * b[j];
        for (unsigned int k = i; k < n; k++) {
          A[j][k] -= a_d * A[i][k];
        }
      }
    }

    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        printf("%lf ", A[i][j]);
      }
      printf("|%lf| ", b[i]);
      printf("\n");
    }

    return back_sub(n, A, b);

  } else {
    return NULL;
  }
}

double *forward_sub(const unsigned int n, double L[n][n], double b[n]) {
  if (fabs(det_T(n, L)) > P) {
    double *x = calloc(n, sizeof(double));
    double res;
    for (unsigned int i = 0; i < n; i++) {
      res = 0;
      for (unsigned int j = 0; j < i; j++) {
        res += L[i][j] * x[j];
      }
      x[i] = (b[i] - res) / L[i][i];
    }
    return x;
  } else {
    printf("det(L) == 0");
    return NULL;
  }
}

void test_back_sub(void) {
  unsigned int n = 3;
  double U[n][n];
  double b[n];
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      U[i][j] = 0;
    }
  }
  int count = 0;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i; j < n; j++) {
      U[i][j] = ++count;
    }
  }
  b[0] = b[2] = 6;
  b[1] = 9;
  double *x = back_sub(n, U, b);
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      printf("%lf ", U[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  for (unsigned int i = 0; i < n; i++) {
    printf("%lf ", b[i]);
  }
  printf("\n\n");
  for (unsigned int i = 0; i < n; i++) {
    assert(fabs(x[i] - 1) < P);
    printf("%lf ", x[i]);
  }
  free(x);
  printf("\n\n\n");
}

void test_for_sub(void) {
  unsigned int n = 3;
  double L[n][n];
  double b[n];
  int count = 0;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      L[i][j] = 0;
    }
  }
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j <= i; j++) {
      L[i][j] = ++count;
    }
  }
  b[0] = 1;
  b[1] = 5;
  b[2] = 15;
  double *x = forward_sub(n, L, b);
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      printf("%lf ", L[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  for (unsigned int i = 0; i < n; i++) {
    printf("%lf ", b[i]);
  }
  printf("\n\n");
  for (unsigned int i = 0; i < n; i++) {
    assert(fabs(x[i] - 1) < P);
    printf("%lf ", x[i]);
  }
  free(x);
  printf("\n");
}

void test_solver(void) {
  unsigned int n = 3;
  double A[n][n];
  double b[n];
  int count = 0;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  double prev = 1;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i; j < n; j++) {
      prev = A[i][j] = prev * ++count;
    }
    b[i] = 1;
  }

  double *x = back_sub(n, A, b);
  for (unsigned int i = 0; i < n; i++) {
    printf("%lf ", x[i]);
  }
}

int main() {
  fesetround(FE_TONEAREST);
  test_solver();
}