#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double *back_sub(const unsigned int n, double U[n][n], double b[n]) {
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
}

double *forward_sub(const unsigned int n, double L[n][n], double b[n]) {
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
}

int main() {
  unsigned int n = 10;
  double U[n][n];
  double L[n][n]; 
  double b[n];
  for (unsigned int i = 0; i < n; i++) {
    b[i] = (double)1;
    for (unsigned int j = i; j < n; j++) {
      U[i][j] = (double)1;
    }
  }
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
    printf("%lf ", x[i]);
  }
  free(x); 
  printf("\n\n\n");
  
  for (unsigned int i = 0; i < n; i++) {
    b[i] = (double)1;
    for (unsigned int j = 0; j <= i; j++) {
      L[i][j] = (double)1;
    }
  }
  double *x1 = forward_sub(n, L, b);
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
    printf("%lf ", x1[i]);
  }
  free(x1); 
  printf("\n");
  
}