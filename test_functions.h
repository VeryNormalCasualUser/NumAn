#include <math.h>

double ex2_3a(double x){
  double x2 = x*x;
  return x2*x - 6*x2 + 4*x + 12;  
}

double ex2_3b(double x){
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  double x6 = x4*x2;
  return exp(pow(sin(x), 3)) + x6 - 2*x4 - x3 - 1;
}

double ex2_3b_unopt(double x){
  return exp(pow(sin(x), 3)) + pow(x,6) - 2*pow(x,4) - pow(x,3) - 1;
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

double ex3_e(double x){
  return pow(x,3) - x;
}

double ex4(double x){
  return exp(x) + sin(x) - 4;
}
