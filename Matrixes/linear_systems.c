#include "../numan_utils.h"
#include "../Inverse-matrix.c"
#include <assert.h>
#include <fenv.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define ITER 50

// determinant for Triangular and Diagonal matrixes
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
  // assumes L is a 0 matrix
  for (unsigned int i = 0; i < n - 1; ++i) {
    double a = A[i][i];
    if (fabs(a) < P) {
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
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
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

void jacobi_iteration(const unsigned int n, double D_1_LpU[n][n], double xk[n],
                      double xk_1[n], double D_1b[n]) {
  double D_1_LpUxk[n];
  print_array_wt(n, xk, "jacobi iteration for xk:");

  for (unsigned int i = 0; i < n; ++i) {
    D_1_LpUxk[i] = 0.0;
    for (unsigned int j = i + 1; j < n; ++j) {
      D_1_LpUxk[i] += (-D_1_LpU[i][j]) * xk[j];
    }
  }
  for (unsigned int i = 1; i < n; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      D_1_LpUxk[i] += (-D_1_LpU[i][j]) * xk[j];
    }
  }

  double left;
  double right;
  for (unsigned int i = 0; i < n; ++i) {
    left = D_1_LpUxk[i];
    right = D_1b[i];
    printf("xk[%u] = %lf + %lf = %lf\n", i, left, right,
           xk_1[i] = left + right);
  }
}

double *jacobi_method(const unsigned int n, int iterations, double A[n][n],
                      double xk[n], double b[n]) {
  double D_1[n];
  double D_1b[n];
  double D_1_LpU[n][n];
  print_matrix_wt(n, A, "A:");
  print_array_wt(n, xk, "x0:");
  print_array_wt(n, b, "b:");
  for (unsigned int i = 0; i < n; ++i) {
    D_1[i] = 1 / A[i][i];
    A[i][i] = 0;
  }

  print_matrix_wt(n, A, "L+U");
  print_array_wt(n, D_1, "D^1");

  for (unsigned int i = 0; i < n; ++i) {
    D_1b[i] = D_1[i] * b[i];
  }

  for (unsigned int i = 0; i < n; ++i) { // optimize
    for (unsigned int j = 0; j < n; ++j) {
      D_1_LpU[i][j] = D_1[i] * A[i][j];
    }
  }

  print_matrix_wt(n, D_1_LpU, "D^-1 * (L+U):");

  print_array_wt(n, D_1b, "D^-1 * b:");
  double *xk_1 = calloc(n, sizeof(double));
  double *temp;

  const int it = iterations; // FOR DEBUG
  while (iterations > 0) {
    jacobi_iteration(n, D_1_LpU, xk, xk_1, D_1b);
    printf("x%d:\n", it - iterations);
    print_array(n, xk);
    printf("x%d:\n", it - iterations + 1);
    print_array(n, xk_1);
    temp = xk;
    xk = xk_1;
    xk_1 = temp;
    temp = NULL;
    for (unsigned int i = 0; i < n; ++i) {
      xk_1[i] = 0.0;
    }
    iterations--;
  }

  double *res = calloc(n, sizeof(double));
  memcpy(res, xk, n * sizeof(double));
  free(xk);
  free(xk_1);
  return res;
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

  double *res = back_sub(n, A, b);
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

void testJacobi2x2(void) {
  unsigned int n = 2;
  double A[n][n];
  double b[n];

  A[0][0] = 2;
  A[0][1] = 1;
  A[1][0] = 3;
  A[1][1] = 4;

  b[0] = 3;
  b[1] = 5;

  double *x0 = calloc(n, sizeof(double));
  x0[0] = 0;
  x0[1] = 1;

  double *res = jacobi_method(n, ITER, A, x0, b);
  assert(fabs(res[0] - 1.4) < P);
  assert(fabs(res[1] - 0.2) < P);
  free(res);
}

void testJacobi3x3(void) {
  unsigned int n = 3;
  double A[n][n];
  double b[n];

  A[0][0] = 3;
  A[0][1] = 1;
  A[0][2] = 1;
  A[1][0] = 1;
  A[1][1] = 3;
  A[1][2] = 1;
  A[2][0] = 1;
  A[2][1] = 1;
  A[2][2] = 3;

  b[0] = 3;
  b[1] = 5;
  b[2] = 7;

  double *x0 = calloc(n, sizeof(double));
  x0[0] = 0.1;
  x0[1] = 1.1;
  x0[2] = 2.1;

  double *res = jacobi_method(n, ITER, A, x0, b);
  for (unsigned int i = 0; i < n; ++i) {
    assert(fabs(res[i] - i) < P);
  }
  free(res);
}

void testJacobi4x4(void) {

  unsigned int n = 4;
  double A[n][n];
  double b[n];

  A[0][0] = 3;
  A[0][1] = 1;
  A[0][2] = 1;
  A[0][3] = 0;

  A[1][0] = 1;
  A[1][1] = 6;
  A[1][2] = 3;
  A[1][3] = -1;

  A[2][0] = 6;
  A[2][1] = 0;
  A[2][2] = 9;
  A[2][3] = -2;

  A[3][0] = 1;
  A[3][1] = 0;
  A[3][2] = -1;
  A[3][3] = -7;

  b[0] = 1;
  b[1] = 1;
  b[2] = 1;
  b[3] = 1;


  double *x0 = calloc(n, sizeof(double));
  for (unsigned int i = 0; i < n; ++i) {
    x0[i] = 1.1;
  }

  double *res = jacobi_method(n, ITER, A, x0, b);
  assert(fabs(norm(n, res)) - 0.385746 < P*10);
  free(res);
}

void testJacobi6x6(void) {

  unsigned int n = 6;
  double A[n][n];
  double b[n];

  A[0][0] = 3;
  A[0][1] = -1;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = 0;
  A[0][5] = 1 / 2;
  A[1][0] = -1;
  A[1][1] = 3;
  A[1][2] = -1;
  A[1][3] = 0;
  A[1][4] = 1 / 2;
  A[1][5] = 0;
  A[2][0] = 0;
  A[2][1] = -1;
  A[2][2] = 3;
  A[2][3] = -1;
  A[2][4] = 0;
  A[2][5] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -1;
  A[3][3] = 3;
  A[3][4] = -1;
  A[3][5] = 0;
  A[4][0] = 0;
  A[4][1] = 1 / 2;
  A[4][2] = 0;
  A[4][3] = -1;
  A[4][4] = 3;
  A[4][5] = -1;
  A[5][0] = 1 / 2;
  A[5][1] = 0;
  A[5][2] = 0;
  A[5][3] = 0;
  A[5][4] = -1;
  A[5][5] = 3;

  b[0] = 5 / 2;
  b[1] = 3 / 2;
  b[2] = 1;
  b[3] = 1;
  b[4] = 3 / 2;
  b[5] = 5 / 2;

  double *x0 = calloc(n, sizeof(double));
  for (unsigned int i = 0; i < n; ++i) {
    x0[i] = 1.1;
  }

  double *res = jacobi_method(n, ITER, A, x0, b);
  for (unsigned int i = 0; i < n; ++i) {
    assert(fabs(res[i] - 1) < P);
  }
  free(res);
}

void testInvers2x2(void){
  unsigned int n = 2;
  double A[n][n];
  double A_1[n][n];
  double A_T[n][n];

  A[0][0] = 1;
  A[0][1] = 0;
  A[1][0] = 3;
  A[1][1] = 4;

  inverse(n, A, determinant(n,A), A_1, A_T);

  printf("A:\n");
  for (unsigned int i = 0; i < n; ++i){
    for (unsigned int j = 0; j < n; ++j){
      printf("%lf ", A[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("A^-1:\n");
  for (unsigned int i = 0; i < n; ++i){
    for (unsigned int j = 0; j < n; ++j){
      printf("%lf ", A_1[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("A^T:\n");
  for (unsigned int i = 0; i < n; ++i){
    for (unsigned int j = 0; j < n; ++j){
      printf("%lf ", A_T[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

int main() {
  fesetround(FE_TONEAREST);
  testJacobi2x2();
  testJacobi3x3();  
  testJacobi6x6();
  testJacobi4x4();
  testInvers2x2(); 
}