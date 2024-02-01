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
    printf("det(U) == 0\n");
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

double det() {
  return 1; // TODO: implement
}

void gauss_elim(const unsigned int n, double A[n][n], double b[n]) {
  for (unsigned int i = 0; i < n - 1; ++i) {
    double a = A[i][i];
    for (unsigned int j = i + 1; j < n; ++j) {
      double d = A[j][i];
      double a_d = d / a;
      b[j] -= a_d * b[i];
      for (unsigned int k = i; k < n; k++) {
        A[j][k] -= a_d * A[i][k];
      }
    }
  }
}

void factor_LU(const unsigned int n, double A[n][n], double L[n][n]) {
  //assumes L is a 0 matrix
  for (unsigned int i = 0; i < n - 1; ++i) {
    double a = A[i][i];
    if (fabs(a) < P){
      printf("pivot is zero\n");
      break;
    }
    for (unsigned int j = i + 1; j < n; ++j) {
      double d = A[j][i];
      double a_d = d / a;
      printf("a/d: %lf\n", a_d);
      L[j][i] = a_d;
      for (unsigned int k = i; k < n; k++) {
        L[k][k] = 1;
        A[j][k] -= a_d * A[i][k];
      }
    }
  }
}

double *solve_LU(const unsigned int n, double A[n][n], double b[n]) {
  double L[n][n];
  for (unsigned int i = 0; i < n; i++){
    for (unsigned int j = 0; j < n; j++){
      L[i][j] = 0;
    }
  }
  printf("A starts as:\n");
  print_matrix(n, A);
  if (fabs(det()) > P) {

    factor_LU(n, A, L);
    printf("A after LU (so U):\n");
    print_matrix(n, A);
    printf("L after LU:\n");
    print_matrix(n, L);

    double *y = forward_sub(n, L, b);
    printf("y\n");
    print_array(n, y);
    double *x = back_sub(n, A, y);

    free(y);
    return x;

  } else {
    return NULL;
  }
}

double *solve(const unsigned int n, double A[n][n], double b[n]) {
  printf("solving with Gauss\n");
  printf("A starts as:\n");
  print_matrix(n, A);
  printf("b starts as: \n");
  print_array(n, b);
  printf("\n");
  if (fabs(det()) > P) {
    gauss_elim(n, A, b);
    printf("A after gauss (so U):\n");
    print_matrix(n, A);
    printf("b after gauss: \n");
    print_array(n, b);
    printf("\n");
    return back_sub(n, A, b);
  } else {
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
    assert(fabs(x[i] - 1) < P);
  }
  free(x);
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
    assert(fabs(x[i] - 1) < P);
  }
  free(x);
}

void test_solver(void) {
  unsigned int n = 2;
  double A[n][n];
  double b[n];
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }

  A[0][0] = 1;
  A[0][1] = 1;
  b[0] = 2;
  A[1][0] = 0;
  A[1][1] = 1;
  b[1] = 1;

  double* res = back_sub(n, A, b);
  printf("back_sub:\n");
  print_array(n, res);
  free(res);
  printf("\n");

  A[0][0] = 1;
  A[0][1] = 0;
  b[0] = 1;
  A[1][0] = 1;
  A[1][1] = 1;
  b[1] = 2;

  res = forward_sub(n, A, b);
  printf("forward_sub:\n");
  print_array(n, res);
  free(res);
  printf("\n");

  A[0][0] = 2;
  A[0][1] = 4;
  b[0] = 1;
  A[1][0] = 7;
  A[1][1] = 1;
  b[1] = 8;

  double *x;
  x = solve(n, A, b);
  printf("solution from Gauss: \n");
  print_array(n, x);
  printf("\n");
  free(x);

  A[0][0] = 2;
  A[0][1] = 4;
  b[0] = 1;
  A[1][0] = 7;
  A[1][1] = 1;
  b[1] = 8;

  printf("b: \n");
  print_array(n, b);
  printf("\n");

  x = solve_LU(n, A, b);
  printf("solution: \n");
  print_array(n, x);
  free(x);
}

int main() {
  fesetround(FE_TONEAREST);
  test_solver();
}