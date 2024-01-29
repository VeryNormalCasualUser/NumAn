#include <math.h>

double ex2_3(double x){
  double x2 = x*x;
  return x2*x - 6*x2 + 4*x + 12;  
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