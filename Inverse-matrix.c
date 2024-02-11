#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void minor(unsigned int colMatrix, unsigned int size, double minorMatrix[size][size],
           double newMinorMatrix[size][size]) {
  unsigned int col, row, row2 = 0, col2 = 0;
  for (row = 1; row < size; row++) {
    for (col = 0; col < size; col++) {
      if (col == colMatrix) {
        continue;
      }
      newMinorMatrix[row2][col2] = minorMatrix[row][col];
      col2++;
      if (col2 == (size - 1)) {
        row2++;
        col2 = 0;
      }
    }
  }
  return;
}

double determinant(unsigned int size, double minorMatrix[size][size]) {
  unsigned int col;
  double sum = 0, newMinorMatrix[size][size];
  if (size == 1) {
    return minorMatrix[0][0];
  } else if (size == 2) {
    return (minorMatrix[0][0] * minorMatrix[1][1] -
            minorMatrix[0][1] * minorMatrix[1][0]);
  } else {
    for (col = 0; col < size; col++) {
      minor(col, size, minorMatrix, newMinorMatrix); // function
      sum +=
          (double)(minorMatrix[0][col] * pow(-1, col) *
                   determinant(size - 1, newMinorMatrix)); // function
    }
  }
  return sum;
}

void transpose(unsigned int size, double cofactorMatrix[size][size],
               double determinte, double coutMatrix[size][size],
               double transposeMatrix[size][size]) {
  unsigned int row, col;
  for (row = 0; row < size; row++) {
    for (col = 0; col < size; col++) {
      transposeMatrix[row][col] = cofactorMatrix[col][row];
      coutMatrix[row][col] =
          cofactorMatrix[col][row] / determinte; // adjoint method
    }
  }
  return;
}

void cofactor(unsigned int size, double cinMatrix[size][size], double determinte,
              double coutMatrix[size][size], double transposeMatrix[size][size]) {
  double minorMatrix[size][size], cofactorMatrix[size][size];
  unsigned int col3, row3, row2, col2, row, col;
  for (row3 = 0; row3 < size; row3++) {
    for (col3 = 0; col3 < size; col3++) {
      row2 = 0;
      col2 = 0;
      for (row = 0; row < size; row++) {
        for (col = 0; col < size; col++) {
          if (row != row3 && col != col3) {
            minorMatrix[row2][col2] = cinMatrix[row][col];
            if (col2 < (size - 2)) {
              col2++;
            } else {
              col2 = 0;
              row2++;
            }
          }
        }
      }
      cofactorMatrix[row3][col3] =
          pow(-1, (row3 + col3)) * determinant((size - 1), minorMatrix);
    }
  }
  transpose(size, cofactorMatrix, determinte, coutMatrix,
            transposeMatrix); // function
  return;
}

void inverse(unsigned int size, double cinMatrix[size][size], double determinte,
             double coutMatrix[size][size],
             double transposeMatrix[size][size]) {
  assert(determinte); 
  if (size == 1) {
    coutMatrix[0][0] = 1;
  } else {
    cofactor(size, cinMatrix, determinte, coutMatrix,
             transposeMatrix); // function
  }
  return;
}