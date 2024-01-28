#include <fenv.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/_types/_sigaltstack.h>
#include "../utils.h"

#define P 0.0000001
#define POW 1000000

typedef struct properties_t {
  bool is_contraction;
  bool is_monotonic;
  bool is_self_mapping;
} properties_t;

properties_t* find_properties(double (*f)(double), double min, double max,
                    const double precision) {
  double max_derivative;
  double res;
  properties_t * props = malloc(sizeof(properties_t));
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
    props->is_monotonic = props->is_monotonic ? !(res < 0) : props->is_monotonic;
    props->is_self_mapping = res > min && res < max && props->is_self_mapping;  
    res = fabs(res);
    max_derivative = res > max_derivative ? res : max_derivative;
  }
  props->is_contraction = max_derivative < 1; 
  return props;
}

double find_std_fixed_point(double (*f)(double), double starting_value, const int iterations, const double precision){
  
  double min = starting_value - (2*precision);
  double max =  starting_value + (2*precision); 
  properties_t * props = find_properties(f, min, max, precision); 
  if (!props->is_contraction){
    printf("Not a contraction within region\n");
    return NAN; 
  }
  int count = 0;
  double prev = starting_value;
  starting_value = (*f)(starting_value); 
  while(fabs(starting_value - prev) > precision && count < iterations-1){
    prev = starting_value; 
    starting_value = (*f)(starting_value);
    ++count; 
  }
  free(props); 
  printf("count: %d \n", count+1); 
  return starting_value; 
}

extern double derivative(double (*)(double), double, double);

bool is_newton_raphson_contraction(double (*f)(double), double starting_value, const double precision){
  double min = starting_value - (POW*precision);
  double max =  starting_value + (POW*precision);
  double x0 = (min - (*f)(min)/derivative(f,min,precision));
  min += precision; 
  double x1 = (min - (*f)(min)/derivative(f,min,precision));
  double max_derivative = fabs((x1-x0)/precision);
  double res; 
  while (max_derivative < 1 && (min += precision) < max){
    x0 = x1; 
    x1 = (min - (*f)(min)/derivative(f,min,precision));
    res = fabs((x1-x0)/precision);
    max_derivative = max_derivative > res ? max_derivative : res;
  }
  return !(max_derivative > 1); 
}


double newton_raphson_phi(double (*f)(double), double starting_value, const double precision){
  double fx = (*f)(starting_value);
  double fprimx = derivative(f, starting_value, precision);
  printf("x=%lf; f(x)=%lf; f'(x)=%lf;\n", starting_value, fx, fprimx); 
  return starting_value - (fx/fprimx);
}

double newton_raphson(double (*f)(double), double starting_value, const int iterations, const double precision, const bool check_contraction){
  double x1;
  int count = 0;
  double prev;
  if (check_contraction && !is_newton_raphson_contraction(f, starting_value, precision)){
    printf("Is not contraction around %lf\n", starting_value);
    return NAN; 
  }
  x1 = newton_raphson_phi(f, starting_value, precision); 
  prev = starting_value;
  starting_value = x1;
  while(fabs(starting_value - prev) > precision && count < iterations-1){
    prev = starting_value; 
    starting_value = newton_raphson_phi(f, starting_value, precision); 
    ++count; 
  }
  printf("count: %d\n", count+1); 
  return starting_value; 
}

double square(double x){
  return x*x; 
}

double cube(double x){
  return x*x*x; 
}

double cosSq(double x){
  double cosx = cos(x);
  return cosx*cosx; 
}

double cust1(double x){
  double x2 = x*x; 
  return x2*x + x2 - 1; 
}

double cust2(double x){
  return x*x + (1/(x+1)) - 3*x; 
}

double cust3(double x){
  return 5*x - 10;
}

double alternate(double x) {
  return x/sqrt(pow(x,2) + 1); 
}

int main(){
  fesetround(FE_TONEAREST);

  double res = find_std_fixed_point(&cosSq, 0.5, 200000, P);
  printf("0 %lf\n", res);
  double res1 = find_std_fixed_point(&alternate, 0.9, 200000, P);
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
}

