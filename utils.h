inline int min(int a, int b) { return a < b ? a : b; }

double derivative(double (*f)(double), double value, double precision){
  double x0 = (*f)(value - precision);
  double x1 = (*f)(value + precision);
  return (x1 - x0)/(2*precision); 
}