#include <math.h>
#define P 1E-8
#define POW 1000000
#include <stdio.h>
#include <strings.h>

inline int min(int a, int b) { return a < b ? a : b; }

double derivative(double (*f)(double), double value, double precision) {
  double x0 = (*f)(value - precision);
  double x1 = (*f)(value + precision);
  return (x1 - x0) / (2 * precision);
}

void print_array_wt(unsigned int n, double b[n], char* str) {
  printf("%s\n", str);
  for (unsigned int i = 0; i < n; ++i) {
    printf("%lf ", b[i]);
  }
  printf("\n\n");
}

void print_array(unsigned int n, double b[n]) {

  for (unsigned int i = 0; i < n; ++i) {
    printf("%lf ", b[i]);
  }
  printf("\n\n");
}

void print_matrix(unsigned int n, double A[n][n]) {
  printf("\n");
  for (unsigned int i = 0; i < n; ++i) {
    print_array(n, A[i]); 
  }
  printf("\n\n");
}

void print_matrix_wt(unsigned int n, double A[n][n], char* str) {
  printf("%s\n", str);
  for (unsigned int i = 0; i < n; ++i) {
    print_array(n, A[i]);
  }
  printf("\n\n");
}

double norm(unsigned int n, double vector[n]){
  double sum = 0;
  double curr;
  for (unsigned int i = 0; i < n; ++i){
    curr = vector[i];
    sum += curr*curr;
  }
  return sqrt(sum);
}
