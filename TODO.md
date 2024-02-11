Better check for contractions in newton raphson method: |(xk - r)| < |(xk+1 - r)| for example
Bigger Test Suite

Implement Regula Falsi

Implement backward error check:
    1. backward error (relative and absolute): Slides Lecture 9 page 16
Implement condition number calculation
Implemet different accuracy at different stages of calculation based on condition number
    Condition number of A := kp(A)
    E.g.: kp(A) = 10^n then n digits will be lost in solving Ax = b
    1 < kp(A) < 10^8 - well conditioned
    10^8 <= kp(A) < 10^16 - ill conditioned
    10^16 <= kp(A)  - completely unreliable
    kp(a) = ||A||p||A^-1||p - expensive to calculate

Implement Gauss and LU factorization by pivoting (use the custom pointer initilization matrix in Jonas Skeppstedts book )
    See L10 slides (page 8)

Implement checking for strict diagonal dominance (for matrixes).
    if a matrix is stricty diagonally dominant then Jacobi and Gauss-Seidel will always converge

implement Jacobi iteration method: (L10 p40)
    -maybe add pivoting

Implement Gauss-Seidel method (L10 p49)

Implement SOR (L10 p60)


General:
    Check for memory leaks using valgrind
    Move out testing funtions to seperate files
    remove prints

