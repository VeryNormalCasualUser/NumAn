#include "../test_functions.h"
#include "../utils.h"
#include <fenv.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/_types/_sigaltstack.h>

#define P 0.0000001
#define POW 1000000

typedef struct properties_t {
  bool is_contraction;
  bool is_monotonic;
  bool is_self_mapping;
} properties_t;

properties_t *find_properties(double (*f)(double), double min, double max,
                              const double precision) {
  double max_derivative;
  double res;
  properties_t *props = malloc(sizeof(properties_t));
  props->is_monotonic = false;
  double x = (*f)(min);
  min += precision;
  double x1 = (*f)(min);
  res = ((x1 - x) / precision);
  max_derivative = fabs(res);
  props->is_monotonic = !(max_derivative < 0);
  props->is_self_mapping = res > min && res < max;
  while (max_derivative < 1 && (min += precision) < max) {
    x = x1;
    x1 = (*f)(min);
    res = ((x1 - x) / precision);
    props->is_monotonic =
        props->is_monotonic ? !(res < 0) : props->is_monotonic;
    props->is_self_mapping = res > min && res < max && props->is_self_mapping;
    res = fabs(res);
    max_derivative = res > max_derivative ? res : max_derivative;
  }
  props->is_contraction = max_derivative < 1;
  return props;
}

double find_std_fixed_point(double (*f)(double), double starting_value,
                            const int iterations, const double precision) {

  double min = starting_value - (2 * precision);
  double max = starting_value + (2 * precision);
  properties_t *props = find_properties(f, min, max, precision);
  if (!props->is_contraction) {
    printf("Not a contraction within region\n");
    return NAN;
  }
  int count = 0;
  double prev = starting_value;
  starting_value = (*f)(starting_value);
  while (fabs(starting_value - prev) > precision && count < iterations - 1) {
    prev = starting_value;
    starting_value = (*f)(starting_value);
    ++count;
  }
  free(props);
  printf("count: %d \n", count + 1);
  return starting_value;
}

extern double derivative(double (*)(double), double, double);

bool is_newton_raphson_contraction(double (*f)(double), double starting_value,
                                   const double precision) {
  double min = starting_value - (POW * precision);
  double max = starting_value + (POW * precision);
  double x0 = (min - (*f)(min) / derivative(f, min, precision));
  min += precision;
  double x1 = (min - (*f)(min) / derivative(f, min, precision));
  double max_derivative = fabs((x1 - x0) / precision);
  double res;
  while (max_derivative < 1 && (min += precision) < max) {
    x0 = x1;
    x1 = (min - (*f)(min) / derivative(f, min, precision));
    res = fabs((x1 - x0) / precision);
    max_derivative = max_derivative > res ? max_derivative : res;
  }
  return !(max_derivative > 1);
}

double newton_raphson_phi(double (*f)(double), double starting_value,
                          const double precision) {
  double fx = (*f)(starting_value);
  double fprimx = derivative(f, starting_value, precision);
  // printf("x=%lf; f(x)=%lf; f'(x)=%lf;\n", starting_value, fx, fprimx);
  return starting_value - (fx / fprimx);
}

double newton_raphson_phi_opt(double (*f)(double), double starting_value,
                              double fx, const double precision, const int m) {
  ;
  double fprimx = derivative(f, starting_value, precision);
  // printf("x=%lf; f(x)=%lf; f'(x)=%lf;\n", starting_value, fx, fprimx);
  return starting_value - m * (fx / fprimx);
}

double newton_raphson_mul(double (*f)(double), double starting_value,
                          const int iterations, const double precision,
                          const int m, const bool check_contraction) {
  int count = 0;
  double fx;
  if (check_contraction &&
      !is_newton_raphson_contraction(f, starting_value, precision)) {
    printf("Is not contraction around %lf\n", starting_value);
    return NAN;
  }
  fx = (*f)(starting_value);
  while (fabs(fx) > precision && count < iterations) {
    starting_value =
        newton_raphson_phi_opt(f, starting_value, fx, precision, m);
    ++count;
    fx = (*f)(starting_value);
  }
  // printf("count: %d\n", count+1);
  return starting_value;
}

double newton_raphson(double (*f)(double), double starting_value,
                      const int iterations, const double precision,
                      const bool check_contraction) {
  return newton_raphson_mul(f, starting_value, iterations, precision, 1,
                            check_contraction);
}

double secant_phi(double xk_1, double xk_2, double fxk_1,
                  double fxk_2, const int m) {
  return xk_1 - m * ((fxk_1 * (xk_1 - xk_2)) / ((fxk_1 - fxk_2)));
}

double secant(double (*f)(double), double xk, int iterations,
              const double precision) {
  double xk_2 = xk;
  double fxk_2 = (*f)(xk_2);
  double xk_1 = xk_2 - fxk_2 / derivative(f, xk_2, precision);
  double fxk_1 = (*f)(xk_1);
  double fx;
  do {
    xk = secant_phi(xk_1, xk_2, fxk_1, fxk_2, 1);
    fx = (*f)(xk);
    iterations--;
    xk_2 = xk_1;
    fxk_2 = fxk_1;
    xk_1 = xk;
    fxk_1 = fx;
  } while (fabs(fx) > precision && iterations > 0);
  return xk;
}

int main() {
  fesetround(FE_TONEAREST);
  /*
  double res = find_std_fixed_point(&cosSq, 0.5, 200000, P);
  printf("0 %lf\n", res);
  double res1 = find_std_fixed_point(&alternate, 0.5, 200000, P);
  printf("1 %lf\n", res1);
  double res2 = find_std_fixed_point(&cube, 0.5, 20000, P);
  printf("2 %lf\n", res2);


  double res3 = newton_raphson(&alternate, 0.9, 200000, P, true);
  printf("3 %lf\n", res3);
  double res4 = newton_raphson(&cube, 0.5, 200000, P, true);
  printf("4 %lf\n", res4);
  double res5 = newton_raphson(&cosSq, 1.5, 200000, P, true);
  printf("5 %lf\n", res5);

  double ex1a = newton_raphson(&cust1, 1, 2, P, false);
  printf("6 %lf\n", ex1a);
  double ex2a = newton_raphson(&cust2, 1, 2, P, false);
  printf("7 %lf\n", ex2a);
  double ex3a = newton_raphson(&cust3, 1, 2, P, false);
  printf("8 %lf\n", ex3a);

  double ex2ex3_1 = newton_raphson(&ex2_3a, 1, POW, P, false);
  printf("9 %lf\n", ex2ex3_1);

  double ex2ex3_2 = newton_raphson(&ex2_3a, 20, POW, P, false);
  printf("10 %lf\n", ex2ex3_2);

  double ex2ex3_3 = newton_raphson(&ex2_3a, -20, POW, P, false);
  printf("11 %lf\n", ex2ex3_3);


  printf("res: %lf\n", ex2_3a(ex2ex3_1));
  printf("res: %lf\n", ex2_3a(ex2ex3_2));
  printf("res: %lf\n", ex2_3a(ex2ex3_3));
  */
  for (double i = -2 ; i < 2.00001; i += 0.01) {
    double xN = newton_raphson(&ex2_3b, i, POW * 10, P, false);
    double xS = secant(&ex2_3b, i, POW*10, P);
    double resN = ex2_3b(xN);
    double resS = ex2_3b(xS);
    printf("N i: %lf, x: %lf, res: %lf\n", i, xN, resN);
    printf("S i: %lf, x: %lf, res: %lf\n", i, xS, resS);
  }
}
