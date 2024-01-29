#include "../range_t.h"
#include "../test_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define P 6
#define MAXCOUNT 1000
extern int min(int,int); 

range_t *bisect(double (*f)(double), double needed_res, double a, double b) {
  range_t *res = malloc(sizeof(range_t));
  int count = 0;
  double fa = needed_res - (*f)(a);
  double fb = needed_res - (*f)(b);
  
  if (!(fa < 0 || fb < 0)) {
    res->min = 0.0 / 0.0;
    res->max = 0.0 / 0.0;
    return res;
  }
  int max_iteration = min((int)(log2((b-a)*pow(10, P)) + 1), MAXCOUNT) + 1;
  double x = (b + a) / 2;
  double fx;

  do {
    fx = needed_res - (*f)(x);
    if (fa * fx < 0) { // branch instead of multiply ?
      b = x;
      fb = needed_res - (*f)(b);
    } else {
      a = x;
      fa = needed_res - (*f)(a);
    }
    x = (b + a) / 2;
    ++count;
    // printf("%d: %lf, %lf \n",count, a, b);
  } while (count < max_iteration && fx != 0);

  if (fx != 0) {
    res->min = a;
    res->max = b;
  } else {
    res->min = res->max = x;
  }
  res->count = count;

  return res;
}


int main() {
  printf("RES:\n");
  double x = 3.0;
  for (int i = 0; i < 20; ++i) {
    range_t *res = bisect(&cube, x, -(x + 1.0), x + 1.0);
    printf("%d, %lf: %lf, %lf; %d \n", i, x, res->min, res->max, res->count);
    free(res);
    x *= 3;
  }
  x = -3.0;
  for (int i = 0; i < 20; ++i) {
    range_t *res = bisect(&cube, x, x - 1.0, -(x - 1.0));
    printf("%d, %lf: %lf, %lf; %d \n", i, x, res->min, res->max, res->count);
    free(res);
    x *= 3;
  }
}